import os
import sys
import tables
from astropy.io import fits
from astropy.io.ascii import SExtractor
import h5py
import numpy as np
from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.constants import d_sun, L_sun, M_sun
from astropy import units as u

from config import config, save_dir, cat_version, home


# M/L calc
## Filter
def calc_MsunFilter(filtercurve):
    '''
        Calculates the solar absolute Magnitude on filter.
    '''
    # The solar spectra is located on the data directory which comes with the pystarlight distribution:
    data_dir = '%s/workspace/pystarlight/data/solar_spectrum/' % home  # os.path.dirname(inspect.getfile(pystarlight)) + '/../../data/solar_spectrum/'
    sun_spec, header_sun_spec = fits.getdata(data_dir + 'sun_reference_stis_002.fits',
                                             header=True)  # Read Solar Spectrum
    Fsun_corr = np.interp(filtercurve['lambda'], sun_spec['WAVELENGTH'], sun_spec['FLUX'])  # Correct to filter lambdas
    # FsunFilter = np.trapz(Fsun_corr * filtercurve['transm'],
    #                       filtercurve[
    #                           'lambda'])  # / np.trapz(filtercurve['transm'], filtercurve['lambda']) # Calc solar luminosity (ergs/s/cm2)
    MsunFilter = -2.5 * np.log10((d_sun ** 2 / ((10 * u.pc).to(u.cm).value ** 2)) * np.trapz(
        filtercurve['lambda'] * filtercurve['transm'] * Fsun_corr, filtercurve['lambda']) / (
                                     np.trapz(filtercurve['transm'] / filtercurve['lambda'],
                                              filtercurve['lambda']))) - 2.41  # Convert: erg/s/cm2 --> L_sun

    return MsunFilter


filter_m2l = readfilterfile(config['filter_m2l'], norm=False).data
aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
filter_m2l_new['lambda'] = aux_l
filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
filter_m2l = filter_m2l_new
MsunFilter = calc_MsunFilter(filter_m2l)

# out_fname = 'alhambra_fit_%s_fl%i.hdf5' % (sys.argv[1].lower(), config['f_l'])
out_fname = 'kk_alhambra_fit_fl%i.hdf5' % config['f_l']
if os.path.exists(out_fname):
    print 'Output file %s already exists' % out_fname
    sys.exit()


f_out = h5py.File(out_fname, 'w')  # FIXME:

f_bpz_fit = h5py.File(config['fit_bpz_outfile'], 'r')
model_parameters = np.copy(f_bpz_fit.get('model_parameters').value)

f_bpz_fit_post = h5py.File(config['fit_bpz_post'], 'r')
likelihood = f_bpz_fit_post.get('likelihood')
model_parameters_post = np.ma.masked_array(np.copy(f_bpz_fit_post.get('parameters').value))

f_cat = SExtractor()
f_cat = f_cat.read('/Users/william/data/alhambra_gold_feb2016/alhambragold_added.cat')
# try:
# f_bpz = h5py.File('/Volumes/unsafe RAID 0/new_pdfs/compressed_alhambra.%s.ColorProBPZ.hdf5' % (sys.argv[1]), 'r')
# except IOError:
#     f_bpz = h5py.File('/Volumes/unsafe RAID 0/new_pdfs/compressed_%s_colorproext_%s_ISO_phz_eB11.hdf5' % (
#         sys.argv[1][:6].lower(), sys.argv[1][-1]), 'r')



# pzt = f_bpz.get('FullProbability')[...likelihood.shape[0], :]
h5file = tables.File('/Users/william/data/alhambra_gold_feb2016/alhambragold_added_%s_1e-4_B13v6_eB11.h5' % cat_version)

pzt = np.zeros((len(h5file.root.Posterior), len(h5file.root.z), len(h5file.root.xt)), "float")
# pzt = np.zeros((n_galaxies, len(h5file.root.z), len(h5file.root.xt)), "float")

n_galaxies = pzt.shape[0]
# n_galaxies = 10

for j, x in enumerate(h5file.root.Posterior[:pzt.shape[0]]):
    gz = h5file.root.goodz[j]
    gt = h5file.root.goodt[j]
    if x.sum() > 0:
        pzt[j][np.outer(gz, gt)] += (x / x.sum())

pzt = pzt[:, :likelihood.shape[0], :]
pzt_reduced = np.zeros((n_galaxies, pzt.shape[1], pzt.shape[2]))
galaxies_reduced = np.zeros(n_galaxies, dtype=int)

i_gal = 0
for i in xrange(len(h5file.root.Posterior)):
    if np.sum(pzt[i]) > .98:
        if int(h5file.root.bpz[i]["ID"]) != f_cat[i]["ID"]:
            print 'error'
            sys.exit(1)
        pzt_reduced[i_gal] = pzt[i]
        galaxies_reduced[i_gal] = i
        print 'i_gal', i_gal, i, f_cat[i]["ID"]
        i_gal += 1
        if i_gal == n_galaxies:
            break

pzt_reduced = pzt_reduced[:i_gal-1]
galaxies_reduced = galaxies_reduced[:i_gal-1]

# pzt = np.ma.masked_equal(pzt[mask], 0)

# Save which galaxies are the ones with prob > 98% within the redshift range of our bpz_fit file
f_out.create_dataset('/gal_alhambra_seq_id', data=galaxies_reduced)

# f_cat = f_cat[mask]
f_out.create_dataset('/bpz_catalog', data=f_cat[galaxies_reduced])  # Only the galaxies with p > 98%

# Convert model_post from ndarray to a common array:
model_post_names = np.sort(model_parameters_post.dtype.names)
model_parameters_post_n = np.zeros((model_parameters_post.shape[0], model_parameters_post.shape[1],
                                    model_parameters_post.shape[2], len(model_post_names)))
for i in range(len(model_post_names)):
    model_parameters_post_n[..., :, i] = model_parameters_post[model_post_names[i]]
model_parameters_post = model_parameters_post_n

# likelihood_flat = likelihood.value.ravel()
# gal_parameters_likelihood = np.zeros((n_galaxies, len(likelihood_flat)))
# gal_parameters_flat = np.zeros(len(likelihood_flat), model_parameters.shape[-1])
# gal_parameters_likelihood = np.zeros((n_galaxies, ))

gal_mass_bins = np.empty((pzt_reduced.shape[0], pzt.shape[1], pzt.shape[2]))

# f_out.create_dataset('/gal_parameters_names', data=model_post_names)
# f_out.create_dataset('/gal_parameters_names_post', data=model_post_names)
# gal_parameters = np.zeros((n_galaxies, model_parameters.shape[3]))
# gal_parameters_post = np.zeros((n_galaxies, model_parameters_post.shape[3]))
gal_parameters_names = f_out.create_dataset('/gal_parameters_names',
                                            data=np.array(
                                                (['a_v', 'at_flux', 'm2l', 'metallicity', 'luminosity', 'mass'])))
n_params = len(gal_parameters_names)
param_resolution = 1000  # FIXME: Defines PDF resolution. 200 points for each parameter.
# gal_parameters_likelihood = np.zeros((n_galaxies, param_resolution-1, n_params))
gal_parameters_likelihood = f_out.create_dataset('/gal_parameters_likelihood',
                                                 shape=(pzt_reduced.shape[0], param_resolution - 1, n_params))
gal_parameters_bins = np.zeros((param_resolution, n_params))
# f_out.create_dataset('/gal_parameters_bins', shape=(param_resolution, n_params))    # To store the bins points...

# To store the bins points...
gal_parameters_bins_ds = f_out.create_dataset('/gal_parameters_bins', shape=(param_resolution - 1, n_params))

# Bins for stellar population parameters
model_parameters_names_ipar = [int(np.argwhere(gal_parameters_names.value[i_par] == model_post_names)) for i_par in
                               range(n_params - 2)]

for i_par in range(n_params - 2):
    gal_parameters_bins[:, i_par] = np.linspace(model_parameters_post[..., :, model_parameters_names_ipar[i_par]].min(),
                                                model_parameters_post[..., :, model_parameters_names_ipar[i_par]].max(),
                                                param_resolution)
    # Store the distribution center:
    gal_parameters_bins_ds[:, i_par] = gal_parameters_bins[:, i_par][1:] - (gal_parameters_bins[:, i_par][
                                                                            1:] - gal_parameters_bins[:,
                                                                                  i_par][:-1]) / 2
# For the galaxy luminosity L
i_par += 1
gal_luminosity = 10 ** (-.4 * (
    h5file.root.Absolute_Magnitude_zT_for_m0eq20[:likelihood.shape[0]] - 20 + \
    h5file.root.bpz[galaxies_reduced]["m0"][:, np.newaxis, np.newaxis] - MsunFilter))
    # f_cat[galaxies_reduced]["F814W"][:, np.newaxis, np.newaxis] - MsunFilter))
gal_parameters_bins[:, i_par] = 10 ** np.linspace(np.log10(np.min(gal_luminosity)), np.log10(np.max(gal_luminosity)),
                                                  param_resolution)
gal_parameters_bins_ds[:, i_par] = gal_parameters_bins[:, i_par][1:] - (gal_parameters_bins[:, i_par][1:] -
                                                                        gal_parameters_bins[:, i_par][:-1]) / 2
# For the galaxy mass
i_par += 1
gal_parameters_bins[:, i_par] = 10 ** np.linspace(3, 13, param_resolution)
gal_parameters_bins_ds[:, i_par] = gal_parameters_bins[:, i_par][1:] - (gal_parameters_bins[:, i_par][1:] -
                                                                        gal_parameters_bins[:, i_par][:-1]) / 2

# plt.clf()

for i_gal in range(pzt_reduced.shape[0]):
    print 'i_gal: ', i_gal
    aux_likelihood = (pzt_reduced[i_gal][..., np.newaxis] * likelihood.value ** config['f_l'])
    aux_likelihood /= aux_likelihood.sum()
    for i_par in range(n_params - 2):
        h, b = np.histogram(model_parameters_post[..., model_parameters_names_ipar[i_par]], weights=aux_likelihood,
                            bins=gal_parameters_bins[:, i_par])
        h /= h.sum()
        gal_parameters_likelihood[i_gal, :, i_par] = h
    # Galaxy luminosities
    i_par += 1
    h, b = np.histogram(gal_luminosity[i_gal], weights=pzt_reduced[i_gal], bins=gal_parameters_bins[:, i_par])
    h /= h.sum()
    gal_parameters_likelihood[i_gal, :, i_par] = h
    # Galaxy mass
    i_par += 1
    aux_mass = model_parameters_post[..., int(np.argwhere(model_post_names == 'm2l'))] * gal_luminosity[
        i_gal, ..., np.newaxis]
    h, b = np.histogram(aux_mass, weights=aux_likelihood, bins=gal_parameters_bins[:, i_par])
    h /= h.sum()
    gal_parameters_likelihood[i_gal, :, i_par] = h

# priors = f_out.create_dataset('/gal_parameters_priors', shape=(param_resolution - 1, n_params))
#
# for i_par in range(n_params - 2):
#     h, b = np.histogram(model_parameters_post[..., :, i_par], bins=gal_parameters_bins[:, i_par])
#     h = np.array(h, dtype='f')
#     h /= h.sum()
#     priors[:, i_par] = h
#     # if i_par == 3:
#     #     plt.plot(b[:-1], h, lw=2, alpha=.5)

print 'done galaxies!'

f_out.close()

print 'done writing to disk! 1'
f_bpz_fit_post.close()
print 'done writing to disk! 2'
f_bpz_fit.close()
print 'done writing to disk! 3'
# f_bpz.close()
