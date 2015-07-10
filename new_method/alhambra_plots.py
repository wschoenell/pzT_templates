import inspect
import os

from astropy.io import fits
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pystarlight
from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.constants import d_sun
from pystarlight.util.constants import L_sun
import sys

__author__ = 'william'

config = {'magerr': .1,
          'filter_m2l': '/Users/william/doutorado/photo_filters/B_Johnson.res',
          'f_l': 1}  # This is a likelihood cooking-factor. The final likelihood is l = l ** f_l

## M2L-related stuff
# Calculate L_sun on filter
def calc_LsunFilter(filtercurve):
        '''
            Calculates the solar luminosity on self.filter in the same way as Bruzual does in Galaxev.
                    NVA helped on this!
        '''
        #The solar spectra is located on the data directory which comes with the pystarlight distribution:
        data_dir = os.path.dirname(inspect.getfile(pystarlight))+'/../../data/solar_spectrum/'
        sun_spec, header_sun_spec = fits.getdata(data_dir+'sun_reference_stis_002.fits', header=True) #Read Solar Spectrum
        Lsun_corr = np.interp(filtercurve['lambda'], sun_spec['WAVELENGTH'], sun_spec['FLUX']) # Correct to filter lambdas
        LsunFilter = np.trapz(Lsun_corr * filtercurve['transm'], filtercurve['lambda']) # Calc solar luminosity (ergs/s/cm2)
        LsunFilter = (LsunFilter * 4*np.pi * d_sun**2 ) / L_sun # Convert: erg/s/cm2 --> L_sun

        return LsunFilter

#from pystarlight, eval l2m on filter
filter_m2l = readfilterfile(config['filter_m2l']).data
aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
filter_m2l_new['lambda'] = aux_l
filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
filter_m2l = filter_m2l_new
LsunFilter = calc_LsunFilter(filter_m2l)

# Calulate "pivotal?" 1/lambda_*
aux_l = np.trapz(filter_m2l['transm']/filter_m2l['lambda'], filter_m2l['lambda']) /\
        np.trapz(filter_m2l['transm']*filter_m2l['lambda'], filter_m2l['lambda'])       # Eq 5. from tmp05.pdf Cids


## Pre-processing: t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity
## Post-processing: a_v, age, m2l, metallicity

alhambra_fields = np.loadtxt(sys.argv[1], dtype='|S20')

avg_age = None
avg_m2l = None
m_abs = None

for field in alhambra_fields:
    f_out = h5py.File('alhambra_fit_%s_fl%i.hdf5' % (field.lower(), config['f_l']), 'r')

    mask = np.bitwise_and(f_out['bpz_catalog']['Stellar_Flag'] < .1, f_out['bpz_catalog']['Odds_1'] > .9)

    tmp_avg_age = [np.average(f_out[u'gal_parameters_bins'][:, 1], weights=f_out[u'gal_parameters_likelihood'][i_gal, :, 1]) for i_gal
               in range(f_out[u'gal_parameters_likelihood'].shape[0])]
    tmp_avg_m2l = [np.average(f_out[u'gal_parameters_bins'][:, 2], weights=f_out[u'gal_parameters_likelihood'][i_gal, :, 2]) for i_gal
               in range(f_out[u'gal_parameters_likelihood'].shape[0])]

    if avg_age is None:
        avg_age = np.array(tmp_avg_age)
        avg_m2l = np.array(tmp_avg_m2l)
        m_abs = f_out['bpz_catalog']['M_ABS_1']
    else:
        avg_age = np.append(avg_age, np.array(tmp_avg_age))
        avg_m2l = np.append(avg_m2l, np.array(tmp_avg_m2l))
        m_abs = np.append(m_abs, f_out['bpz_catalog']['M_ABS_1'])

plt.clf()
plt.hist(avg_age, bins=50, alpha=.8)
plt.xlabel('mean age')
plt.savefig('fig1.png')
raw_input('Next plot...')


f_gal = 10 ** (-.4 * (m_abs + 2.41)) / LsunFilter # Galaxies fluxes in F_\odot

mass = np.log10(10**avg_m2l * f_gal)

plt.clf()
plt.hexbin(mass, avg_age, bins='log')
plt.xlabel('mass')
plt.ylabel('mean age')
plt.savefig('fig2.png')
raw_input('Next...')


# plt.clf()
# plt.plot((avg_m2l*10**(-.4 * f_out['bpz_catalog']['M_ABS_1'])/c.L_sun.cgs.value)[mask], np.array(avg_age)[mask], '.')
#
# plt.clf()
# plt.plot(f_out['bpz_catalog']['Stell_Mass_1'][mask], np.array(avg_age)[mask], '.')

# plt.plot(f['bpz_catalog']['zb_max_1'], f['gal_parameters_post'][:, 1], '.', color='blue')
# plt.plot(f['bpz_catalog']['zb_max_1'][mask], f['gal_parameters_post'][:, 1][mask], '.', color='red')
#
#
# plt.figure(2)
# plt.clf()
# plt.plot(f['gal_parameters_post'][:,2], f['gal_parameters_post'][:,1], '.')