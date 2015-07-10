import os
import sys

from astropy.io.ascii import SExtractor

__author__ = 'william'
import h5py

import numpy as np


config = {'magerr': .1,
          'f_l': 1}  # This is a likelihood cooking-factor. The final likelihood is l = l ** f_l


out_fname = 'alhambra_fit_%s_fl%i.hdf5' % (sys.argv[1].lower(), config['f_l'])
if os.path.exists(out_fname):
    print 'Output file already exists'
    sys.exit()


f_out = h5py.File(out_fname, 'w')  # FIXME:

f_bpz_fit = h5py.File('bpz_fit_magerr_%3.2f.hdf5' % config['magerr'], 'r')
model_parameters = np.copy(f_bpz_fit.get('model_parameters').value)

f_bpz_fit_post = h5py.File('bpz_fit_magerr_%3.2f_post.hdf5' % config['magerr'], 'r')
likelihood = f_bpz_fit_post.get('likelihood')
model_parameters_post = np.ma.masked_array(np.copy(f_bpz_fit_post.get('parameters').value))

f_cat = SExtractor()
f_cat = f_cat.read(
    '/Users/william/doutorado/Alhambra/catalogs/latest/alhambra.%s.ColorProBPZ.cat' % sys.argv[1].upper())
try:
    f_bpz = h5py.File('/Volumes/unsafe RAID 0/new_pdfs/compressed_%s_colorproext_%s_ISO.hdf5' % (
        sys.argv[1][:6].lower(), sys.argv[1][-1]), 'r')
except IOError:
    f_bpz = h5py.File('/Volumes/unsafe RAID 0/new_pdfs/compressed_%s_colorproext_%s_ISO_phz_eB11.hdf5' % (
        sys.argv[1][:6].lower(), sys.argv[1][-1]), 'r')

pzt = f_bpz.get('FullProbability')[:, :likelihood.shape[0], :]
# mask = np.sum(pzt, axis=(1, 2)) > .95
# pzt = np.ma.masked_equal(pzt[mask], 0)

# f_cat = f_cat[mask]
f_out.create_dataset('/bpz_catalog', data=f_cat)

# Convert model_post from ndarray to a common array:
model_post_names = np.sort(model_parameters_post.dtype.names)
model_parameters_post_n = np.zeros((model_parameters_post.shape[0], model_parameters_post.shape[1],
                                    model_parameters_post.shape[2], len(model_post_names)))
for i in range(len(model_post_names)):
    model_parameters_post_n[:, :, :, i] = model_parameters_post[model_post_names[i]]
model_parameters_post = model_parameters_post_n

# likelihood_flat = likelihood.value.ravel()
# gal_parameters_likelihood = np.zeros((pzt.shape[0], len(likelihood_flat)))
# gal_parameters_flat = np.zeros(len(likelihood_flat), model_parameters.shape[-1])
# gal_parameters_likelihood = np.zeros((pzt.shape[0], ))


# f_out.create_dataset('/gal_parameters_names', data=model_post_names)
# f_out.create_dataset('/gal_parameters_names_post', data=model_post_names)
# gal_parameters = np.zeros((pzt.shape[0], model_parameters.shape[3]))
# gal_parameters_post = np.zeros((pzt.shape[0], model_parameters_post.shape[3]))
gal_parameters_names = f_out.create_dataset('/gal_parameters_names', data=np.array(np.sort(['a_v', 'age', 'm2l', 'metallicity'])))
n_params = len(gal_parameters_names)
param_resolution = 1000  #FIXME: Defines PDF resolution. 200 points for each parameter.
# gal_parameters_likelihood = np.zeros((pzt.shape[0], param_resolution-1, n_params))
gal_parameters_likelihood = f_out.create_dataset('/gal_parameters_likelihood', shape=(pzt.shape[0], param_resolution-1, n_params))
gal_parameters_bins = np.zeros((param_resolution, n_params)) #f_out.create_dataset('/gal_parameters_bins', shape=(param_resolution, n_params))    # To store the bins points...

# To store the bins points...
aux = f_out.create_dataset('/gal_parameters_bins', shape=(param_resolution-1, n_params))

for i_par in range(n_params):
    gal_parameters_bins[:, i_par] = np.linspace(model_parameters_post[:, :, :, i_par].min(),
                                                model_parameters_post[:, :, :, i_par].max(), param_resolution)
    # Store the distribution center:
    aux[:, i_par] = gal_parameters_bins[:, i_par][1:] - (gal_parameters_bins[:, i_par][1:] - gal_parameters_bins[:, i_par][:-1])/2


# plt.clf()
for i_gal in range(pzt.shape[0]):
    print 'i_gal: ', i_gal
    aux = (pzt[i_gal][:, :, np.newaxis] * likelihood.value ** config['f_l']) / np.sum(pzt[i_gal])
    for i_par in range(n_params):
        h, b = np.histogram(model_parameters_post[:, :, :, i_par], weights=aux, bins=gal_parameters_bins[:, i_par])
        gal_parameters_likelihood[i_gal, :, i_par] = h
        h /= h.sum()
        # if i_par == 3:
        #     plt.plot(b[:-1], h)


priors = f_out.create_dataset('/gal_parameters_priors', shape=(param_resolution-1, n_params))

for i_par in range(n_params):
    h, b = np.histogram(model_parameters_post[:, :, :, i_par], bins=gal_parameters_bins[:, i_par])
    h = np.array(h, dtype='f')
    h /= h.sum()
    priors[:, i_par] = h
    # if i_par == 3:
    #     plt.plot(b[:-1], h, lw=2, alpha=.5)


print 'done galaxies!'

f_out.close()

print 'done writing to disk! 1'
f_bpz_fit_post.close()
print 'done writing to disk! 2'
f_bpz_fit.close()
print 'done writing to disk! 3'
f_bpz.close()