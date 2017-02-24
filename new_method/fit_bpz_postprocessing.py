import os

from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.base import StarlightBase
from pystarlight.util.redenninglaws import Cardelli_RedLaw

from common import calc_LsunFilter
from model import Model

__author__ = 'william'

import h5py
import numpy as np
# TODO: Import config from hdf5 output.
## 0 - Configuration
from config import config, save_dir


class PostParameters(object):
    def __init__(self, bt, l2m_norm_ssp, l2m_ssp=None, l2m_av_sampling=None):
        self.bt = bt
        self.l2m_norm_ssp = l2m_norm_ssp
        ## Metallicity
        aux = np.log10(self.bt.metBase)
        aux2 = (aux[1:] - aux[0:-1]) / 2
        self.met_min = aux - np.append(aux2[0], aux2)
        self.met_min[0] += -0.00001
        self.met_max = aux + np.append(aux2, aux2[-1])
        self.met_max[-1] += 0.00001
        self.l2m_ssp = l2m_ssp
        self.av_sampling = l2m_av_sampling
        # self.l2m_ssp_spec = l2m_ssp_spec
        # self.met_low = min(self.met_min)
        # self.met_upp = max(self.met_max)

    def get_parameters(self, t0_young, tauT_young, t0_old, tauT_old, lg_frac_young, a_v, metallicity):
        i_met = np.argwhere(np.bitwise_and(metallicity >= self.met_min, metallicity <= self.met_max))[0, 0]
        aux_p = {'a_v': a_v, 'metallicity': self.bt.metBase[i_met]}

        sfh = Model.get_sfh(bt, 10 ** t0_young, 10 ** (tauT_young + t0_young), 10 ** t0_old, 10 ** (tauT_old + t0_old),
                            10 ** lg_frac_young)

        aux_p['at_flux'] = np.sum(sfh * self.l2m_norm_ssp[i_met] * np.log10(self.bt.ageBase)) / np.sum(
            sfh * self.l2m_norm_ssp[i_met])

        aux_p['at_mass'] = np.sum(sfh * np.log10(self.bt.ageBase)) / np.sum(sfh)

        i_av = np.searchsorted(self.av_sampling, a_v)
        if i_av >= len(self.av_sampling):
            print a_v, i_av, l2m_ssp.shape
            i_av -= 1

        aux_p['m2l'] = 1 / np.sum(self.l2m_ssp[i_met, :, i_av] * sfh)

        # aux_p['m2l_spec'] = [1 / np.average(self.l2m_ssp_spec[i_met, :, i], weights=sfh) for i in
        #                      range(len(config["m2l_lambdas"]))]

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
    dtype=np.dtype([('at_mass', np.float), ('at_flux', np.float), ('metallicity', np.float), ('a_v', np.float),
                    ('m2l', np.float)]),
    shape=(n_z, n_t, n_models))

# m2l_spec = np.zeros((n_z, n_t, n_models, len(config["m2l_lambdas"])))

bt = StarlightBase(config['base_file'], config['base_path'])
bt.ageBase[bt.ageBase == 0] = 1e5  # Avoid NaN on log10(age)
i_norm = np.argmin((bt.l_ssp - config['lambda_norm']) ** 2)
# i_met = int(np.argwhere(bt.metBase == metallicity))
l2m_norm_ssp = np.empty((len(bt.metBase), len(bt.ageBase)))
for i_met in range(len(bt.metBase)):
    for i_age in range(len(bt.ageBase)):
        l2m_norm_ssp[i_met, i_age] = bt.f_ssp[i_met, i_age, i_norm]

# from pystarlight, eval l2m on filter
filter_m2l = readfilterfile(config['filter_m2l'], ).data
aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
filter_m2l_new['lambda'] = aux_l
filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
filter_m2l = filter_m2l_new
LsunFilter = calc_LsunFilter(filter_m2l)

aux_m2l_lambdas_i = [int(np.argwhere(bt.l_ssp == i)) for i in config['m2l_lambdas']]

av_step = 1e-3
av_sampling = np.arange(config['AV_min'], config['AV_max'], av_step)
av_sampling = np.append(av_sampling, config['AV_max']) # to deal with boundaries
q_lambda = Cardelli_RedLaw(filter_m2l['lambda'])

# l2m_ssp_spec = np.zeros((len(bt.metBase), len(bt.ageBase), len(config["m2l_lambdas"])))
cache_fname = "cache_M2L_%s_%s_avmin%3.2f_avmax%3.2f_avstep%6.4f.npz" % (
os.path.basename(config['base_file']), os.path.basename(config['filter_m2l']), config['AV_min'], config['AV_max'], 1e-3)
if os.path.exists(cache_fname):
    l2m_ssp = np.load(cache_fname)["arr_0"]
else:
    l2m_ssp = np.empty((len(bt.metBase), len(bt.ageBase), len(av_sampling)))
    for i_met in range(len(bt.metBase)):
        for i_age in range(len(bt.ageBase)):
            print i_met, i_age
            for i_av in range(len(av_sampling)):
                gamma = np.interp(filter_m2l['lambda'], bt.l_ssp, bt.f_ssp[i_met, i_age])
                l2m_ssp[i_met, i_age, i_av] = np.trapz(gamma * 10 ** (-0.4 * q_lambda) * filter_m2l['transm'],
                                                       filter_m2l['lambda']) / (bt.Mstars[i_met, i_age] * LsunFilter)
                # l2m_ssp_spec[i_met, i_age] = bt.f_ssp[
                #     i_met, i_age, aux_m2l_lambdas_i]  # bt.Mstars[i_met, i_age] / bt.f_ssp[i_met, i_age, aux_m2l_lambdas_i]
    np.savez(cache_fname, l2m_ssp)

post_parameters = PostParameters(bt, l2m_norm_ssp, l2m_ssp, av_sampling)

for i_z in range(n_z):
    print 'i_z', i_z
    for i_t in range(n_t):
        print 'i_t', i_t
        for i_model in range(n_models):
            for k, v in post_parameters.get_parameters(*params[i_z, i_t, i_model]).iteritems():
                # if k == "m2l_spec":
                #     m2l_spec[i_z, i_t, i_model] = v
                # else:
                post_parameters_data[i_z, i_t, i_model][k] = v

file_fit_processed.create_dataset('/parameters', data=post_parameters_data, compression='gzip')
l = np.exp(file_fit.get('/model_lnprob') - np.max(file_fit.get('/model_lnprob'), axis=2)[:, :, np.newaxis])
l /= l.sum(axis=2)[:, :, np.newaxis]
file_fit_processed.create_dataset('/likelihood', data=l)
# file_fit_processed.create_dataset('/m2l_spec', data=m2l_spec)
config['filters'] = ','.join(sorted(config['filters'].keys()))
config['m2l_lambdas'] = str(config['m2l_lambdas'])
file_fit_processed.attrs.update(config)

#
file_fit_processed.close()
file_fit.close()
