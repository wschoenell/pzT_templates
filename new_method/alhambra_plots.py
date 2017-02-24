import inspect
import os
from astropy.io import fits
import h5py
from magal.plots.mosaic import get_mosaic
import matplotlib.pyplot as plt
import numpy as np
import pystarlight
from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.constants import d_sun
from pystarlight.util.constants import L_sun
from astropy.io.ascii import SExtractor
from config import config

__author__ = 'william'

# config = {'magerr': .1,
#           'filter_m2l': '/Users/william/doutorado/photo_filters/B_Johnson.res',
#           'f_l': 1}  # This is a likelihood cooking-factor. The final likelihood is l = l ** f_l


def plot_hist(x, axis, color='black', bins=50, label=None, max=None):
    aux_hst, aux_bins = np.histogram(x, bins=bins)
    if max:
        aux_hst = np.array(aux_hst, dtype=np.float)
        aux_hst /= aux_hst.max()
        aux_hst *= max
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    axis.plot(pts, aux_hst, linewidth=2, color=color, label=label)


## M2L-related stuff
# Calculate L_sun on filter
def calc_LsunFilter(filtercurve):
    '''
        Calculates the solar luminosity on self.filter in the same way as Bruzual does in Galaxev.
                NVA helped on this!
    '''
    # The solar spectra is located on the data directory which comes with the pystarlight distribution:
    # data_dir = os.path.dirname(inspect.getfile(pystarlight)) + '/../../data/solar_spectrum/'
    data_dir = '/Users/william/workspace/pystarlight/data/solar_spectrum/'
    sun_spec, header_sun_spec = fits.getdata(data_dir + 'sun_reference_stis_002.fits',
                                             header=True)  # Read Solar Spectrum
    Lsun_corr = np.interp(filtercurve['lambda'], sun_spec['WAVELENGTH'], sun_spec['FLUX'])  # Correct to filter lambdas
    LsunFilter = np.trapz(Lsun_corr * filtercurve['transm'],
                          filtercurve['lambda'])  # Calc solar luminosity (ergs/s/cm2)
    LsunFilter = (LsunFilter * 4 * np.pi * d_sun ** 2) / L_sun  # Convert: erg/s/cm2 --> L_sun

    return LsunFilter


# from pystarlight, eval l2m on filter
filter_m2l = readfilterfile(config['filter_m2l']).data
aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
filter_m2l_new['lambda'] = aux_l
filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
filter_m2l = filter_m2l_new
LsunFilter = calc_LsunFilter(filter_m2l)

# Calulate "pivotal?" 1/lambda_*
aux_l = np.trapz(filter_m2l['transm'] / filter_m2l['lambda'], filter_m2l['lambda']) / \
        np.trapz(filter_m2l['transm'] * filter_m2l['lambda'], filter_m2l['lambda'])  # Eq 5. from tmp05.pdf Cids

## Pre-processing: t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity
## Post-processing: a_v, age, m2l, metallicity

alhambra_fields = {'f02': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]],  # field, pointing, ccd
                   'f03': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]],
                   'f04': [[1, [1, 2, 3, 4]]],
                   'f05': [[1, [1, 2, 3, 4]]],
                   'f06': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]],
                   'f08': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]]}
# 'f07': [[3, [1, 2, 3, 4]], [4, [1, 2, 3, 4]]],

n_galaxies = 10000

f_color = {'f02': 'red', 'f03': 'black', 'f04': 'green', 'f05': 'blue', 'f06': 'brown', 'f07': 'purple',
           'f08': 'yellow'}

avg_age = None
avg_m2l = None
m_abs = None

fig = get_mosaic(2, 2, i_fig=1)

st = SExtractor()

aux_i = 0
redshift_bins = [[0, 0.3], [0.3, 0.6], [0.6, 0.9], [0.9, 1.2]]
# for field in alhambra_fields.keys():
f_avg_age = None
f_avg_m2l = None
f_m_abs = None
f_z_ml = None
f_mask = None
# for aux in alhambra_fields[field]:
#     pointing = aux[0]
#     for ccd in aux[1]:
out_fname = 'alhambra_fit_fl%i.hdf5' % config['f_l']
f_out = h5py.File(out_fname)
bpz_catalog1 = st.read('/Users/william/data/alhambra_gold_feb2016/alhambragold_added.cat')[
               :n_galaxies]  # '/Users/william/doutorado/Alhambra/catalogs/latest/alhambra.%sP0%iC0%i.ColorProBPZ.cat' % (field.upper(), pointing, ccd))
bpz_catalog2 = h5py.File('/Users/william/data/alhambra_gold_feb2016/alhambragold_added_4f254712ff_1e-4_B13v6_eB11.h5')[
                   'bpz'][:n_galaxies]

gal_parameters_likelihood = np.copy(f_out['gal_parameters_likelihood'])
gal_parameters_bins = np.copy(f_out['gal_parameters_bins'])

mask = np.bitwise_and(bpz_catalog1['stell'] < .55, bpz_catalog2['odds'] > .3)
mask = np.bitwise_and(mask, bpz_catalog1['F814W'] > 23.5)

tmp_avg_age = [
    np.ma.average(gal_parameters_bins[:, 1],
                  weights=np.ma.masked_invalid(gal_parameters_likelihood[i_gal, :, 1]) / np.ma.masked_invalid(
                      gal_parameters_likelihood[i_gal, :, 1]).sum()) for
    i_gal in range(gal_parameters_likelihood.shape[0])]

tmp_avg_m2l = [
    np.ma.average(gal_parameters_bins[:, 2],
                  weights=np.ma.masked_invalid(gal_parameters_likelihood[i_gal, :, 2]) / np.ma.masked_invalid(
                      gal_parameters_likelihood[i_gal, :, 2]).sum()) for
    i_gal in range(gal_parameters_likelihood.shape[0])]

if avg_age is None:
    avg_age = np.array(tmp_avg_age)
    avg_m2l = np.array(tmp_avg_m2l)
    m_abs = np.array(bpz_catalog2['Mabs'].copy())
    z_ml = np.array(bpz_catalog2['zml'].copy())
else:
    avg_age = np.append(avg_age, np.array(tmp_avg_age))
    avg_m2l = np.append(avg_m2l, np.array(tmp_avg_m2l))
    m_abs = np.append(m_abs, np.array(bpz_catalog2['Mabs'].copy()))
    z_ml = np.append(z_ml, np.array(bpz_catalog2['zml'].copy()))

if f_avg_age is None:
    f_avg_age = np.array(tmp_avg_age)
    f_avg_m2l = np.array(tmp_avg_m2l)
    f_m_abs = np.array(bpz_catalog2['Mabs'].copy())
    f_z_ml = np.array(bpz_catalog2['zml'].copy())
    f_mask = mask
else:
    f_avg_age = np.append(f_avg_age, np.array(tmp_avg_age))
    f_avg_m2l = np.append(f_avg_m2l, np.array(tmp_avg_m2l))
    f_m_abs = np.append(f_m_abs, bpz_catalog2['Mabs'])
    f_z_ml = np.append(f_z_ml, bpz_catalog2['zml'])
    f_mask = np.append(f_mask, mask)

i_ax = 0
for aux_z in redshift_bins:
    mask_z = np.bitwise_and(f_z_ml >= aux_z[0], f_z_ml < aux_z[1])
    mask_z = np.bitwise_and(f_mask, mask_z)
    ax = fig.axes[i_ax]
    plot_hist(f_avg_age[mask_z], ax, bins=np.arange(7, 10, .05))  # , color=f_color[field], label=field)
    # ax.hist(f_avg_age, bins=np.arange(7, 10, .05), alpha=.8)
    i_ax += 1

i_ax = 0
for aux_z in redshift_bins:
    ax = fig.axes[i_ax]
    if i_ax == 1:
        ax.legend()
    ax.text(7.3, ax.get_ylim()[0] + .9 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            '%3.2f <= z <  %3.2f' % (aux_z[0], aux_z[1]))
    ax.axes.get_yaxis().set_visible(False)
    i_ax += 1

i_ax = 0
for aux_z in redshift_bins:
    mask_z = np.bitwise_and(z_ml >= aux_z[0], z_ml < aux_z[1])
    ax = fig.axes[i_ax]
    plot_hist(avg_age[mask_z], ax, bins=np.arange(7, 10, .05), color='orange', max=ax.get_ylim()[1])
    i_ax += 1

plt.title('avg_age')
plt.savefig('fig1.png')

raw_input('Next plot...')

# plt.clf()
# plt.hist(avg_age, bins=np.arange(7, 10, .05), alpha=.8)
# plt.xlabel('mean age')
# plt.savefig('fig1.png')
# raw_input('Next plot...')
#

# f_gal = 10 ** (-.4 * (m_abs + 2.41)) / LsunFilter  # Galaxies fluxes in F_\odot
#
# mass = np.log10(10 ** avg_m2l * f_gal)
#
# plt.clf()
# plt.hexbin(mass, avg_age, bins='log')
# plt.xlabel('mass')
# plt.ylabel('mean age')
# plt.savefig('fig2.png')
# raw_input('Next...')


f_gal = 10 ** (-.4 * (m_abs + 2.41)) / LsunFilter  # Galaxies fluxes in F_\odot
mass = np.log10(10 ** avg_m2l * f_gal)

i_ax = 0
fig = get_mosaic(2, 2, i_fig=2, y_shareaxis=True)
x_min, x_max, y_min, y_max = np.nanmin(mass), np.nanmax(mass), avg_age.min(), avg_age.max()
for aux_z in redshift_bins:
    mask_z = np.bitwise_and(z_ml >= aux_z[0], z_ml < aux_z[1])
    ax = fig.axes[i_ax]
    ax.hexbin(mass[mask_z], avg_age[mask_z], bins='log', extent=[x_min, x_max, y_min, y_max])
    ax.text(ax.get_xlim()[0] + .08 * (ax.get_xlim()[1] - ax.get_xlim()[0]), 9.5,
            '%3.2f <= z <  %3.2f' % (aux_z[0], aux_z[1]), color='white')
    i_ax += 1

# plt.xlabel('mass')
plt.title('mean age')
plt.savefig('fig2.png')
raw_input('Next...')


# plt.clf()
# plt.plot((avg_m2l*10**(-.4 * bpz_catalog['M_ABS_1'])/c.L_sun.cgs.value)[mask], np.array(avg_age)[mask], '.')
#
# plt.clf()
# plt.plot(bpz_catalog['Stell_Mass_1'][mask], np.array(avg_age)[mask], '.')

# plt.plot(f['bpz_catalog']['zb_max_1'], f['gal_parameters_post'][:, 1], '.', color='blue')
# plt.plot(f['bpz_catalog']['zb_max_1'][mask], f['gal_parameters_post'][:, 1][mask], '.', color='red')
#
#
# plt.figure(2)
# plt.clf()
# plt.plot(f['gal_parameters_post'][:,2], f['gal_parameters_post'][:,1], '.')
