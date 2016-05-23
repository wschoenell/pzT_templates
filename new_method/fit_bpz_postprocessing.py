import inspect
import os

from astropy.io import fits
from magal.util.stellarpop import n_component
from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.base import StarlightBase, pystarlight
from pystarlight.util.constants import d_sun, L_sun


__author__ = 'william'

import h5py
import numpy as np

# TODO: Import config from hdf5 output.
## 0 - Configuration
from config import config, save_dir


def calc_LsunFilter(filtercurve):
    '''
        Calculates the solar luminosity on self.filter in the same way as Bruzual does in Galaxev.
                NVA helped on this!
    '''
    # The solar spectra is located on the data directory which comes with the pystarlight distribution:
    data_dir = '/Users/william/workspace/pystarlight/data/solar_spectrum/' #os.path.dirname(inspect.getfile(pystarlight)) + '/../../data/solar_spectrum/'
    sun_spec, header_sun_spec = fits.getdata(data_dir + 'sun_reference_stis_002.fits',
                                             header=True)  # Read Solar Spectrum
    Lsun_corr = np.interp(filtercurve['lambda'], sun_spec['WAVELENGTH'], sun_spec['FLUX'])  # Correct to filter lambdas
    LsunFilter = np.trapz(Lsun_corr * filtercurve['transm'],
                          filtercurve['lambda'])  # Calc solar luminosity (ergs/s/cm2)
    LsunFilter = (LsunFilter * 4 * np.pi * d_sun ** 2) / L_sun  # Convert: erg/s/cm2 --> L_sun

    return LsunFilter


class PostParameters(object):
    def __init__(self, bt, l2m_norm_ssp, m2l_ssp=None):
        self.bt = bt
        self.l2m_norm_ssp = l2m_norm_ssp
        ## Metallicity
        aux = np.log10(self.bt.metBase)
        aux2 = (aux[1:] - aux[0:-1]) / 2
        self.met_min = aux - np.append(aux2[0], aux2)
        self.met_min[0] += -0.00001
        self.met_max = aux + np.append(aux2, aux2[-1])
        self.met_max[-1] += 0.00001
        self.m2l_ssp = m2l_ssp
        # self.met_low = min(self.met_min)
        # self.met_upp = max(self.met_max)

    def get_parameters(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity):
        i_met = np.argwhere(np.bitwise_and(metallicity >= self.met_min, metallicity <= self.met_max))[0, 0]
        aux_p = {'a_v': a_v, 'metallicity': self.bt.metBase[i_met]}
        csp_model = n_component(self.bt.ageBase)
        csp_model.add_exp(t0_young, tau_young, frac_young)
        csp_model.add_exp(t0_old, tau_old, 1 - frac_young)

        sfh = csp_model.get_sfh()
        aux_p['age'] = np.sum(sfh * self.l2m_norm_ssp[i_met] * np.log10(self.bt.ageBase)) / np.sum(
            sfh * self.l2m_norm_ssp[i_met])

        aux_p['m2l'] = np.average(self.m2l_ssp[i_met], weights=sfh)

        return aux_p

# file_fit = h5py.File('/Users/william/Downloads/bpz_fit_full_newmask_nointerp.hdf5', 'r')
# file_fit = h5py.File('/Users/william/Downloads/bpz_fit_nointerp_newmask_chi2tx.hdf5', 'r')
# file_fit = h5py.File('/Users/william/tmp_out/out_bpzfit/bpz_fit_nointerp_newmask_chi2tx_CCM_test.hdf5', 'r')
file_fit = h5py.File('%s/%s' % (save_dir, config['fit_bpz_outfile']), 'r')

# Doing this copy due bug: https://github.com/h5py/h5py/issues/480
params = np.copy(file_fit.get('/model_parameters').value)

n_z, n_t, n_models, n_parameters = params.shape

file_fit_processed = h5py.File('%s/%s' % (save_dir, config['fit_bpz_post']), 'w')
post_parameters_data = np.empty(
    dtype=np.dtype([('age', np.float), ('metallicity', np.float), ('a_v', np.float), ('m2l', np.float)]),
    shape=(n_z, n_t, n_models))

bt = StarlightBase(config['base_file'], config['base_path'])
bt.ageBase[bt.ageBase == 0] = 100  # Avoid NaN on log10(age)
i_norm = np.argmin((bt.l_ssp - config['lambda_norm']) ** 2)
# i_met = int(np.argwhere(bt.metBase == metallicity))
l2m_norm_ssp = np.empty((len(bt.metBase), len(bt.ageBase)))
for i_met in range(len(bt.metBase)):
    for i_age in range(len(bt.ageBase)):
        l2m_norm_ssp[i_met, i_age] = bt.f_ssp[i_met, i_age, i_norm]

# from pystarlight, eval l2m on filter
filter_m2l = readfilterfile(config['filter_m2l']).data
aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
filter_m2l_new['lambda'] = aux_l
filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
filter_m2l = filter_m2l_new
LsunFilter = calc_LsunFilter(filter_m2l)

m2l_ssp = np.empty((len(bt.metBase), len(bt.ageBase)))
for i_met in range(len(bt.metBase)):
    for i_age in range(len(bt.ageBase)):
        aux = np.interp(filter_m2l['lambda'], bt.l_ssp, bt.f_ssp[i_met, i_age])
        m2l_ssp[i_met, i_age] = bt.Mstars[i_met, i_age] * LsunFilter / np.trapz(aux * filter_m2l['transm'],
                                                                                filter_m2l['lambda'])

post_parameters = PostParameters(bt, l2m_norm_ssp, m2l_ssp)

for i_z in range(n_z):
    print 'i_z', i_z
    for i_t in range(n_t):
        print 'i_t', i_t
        for i_model in range(n_models):
            for k, v in post_parameters.get_parameters(*params[i_z, i_t, i_model]).iteritems():
                post_parameters_data[i_z, i_t, i_model][k] = v

file_fit_processed.create_dataset('/parameters', data=post_parameters_data, compression='gzip')
l = np.exp(file_fit.get('/model_lnprob') - np.max(file_fit.get('/model_lnprob'), axis=2)[:, :, np.newaxis])
l /= l.sum(axis=2)[:, :, np.newaxis]
file_fit_processed.create_dataset('/likelihood', data=l)

#
file_fit_processed.close()
file_fit.close()
