from magal.util.stellarpop import n_component
from pystarlight.util.base import StarlightBase

__author__ = 'william'

import h5py
import numpy as np

# TODO: Import config from hdf5 output.
## 0 - Configuration
config = {'bpz_library': '../templates/eB11.list', 'n_interpolations': 0, 'bpz_library_dir': '../no_elines/',
          'z_ini': 1e-4, 'z_fin': 7.0, 'z_delta': 0.001,
          'base_file': '/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'base_path': 'Base.bc03.Padova1994.chab.All',
          'AV_min': 0, 'AV_max': 2,
          'taylor_file': '/Users/william/tmp_out/m_a.hdf5',
          'filters_dir': '/Users/william/doutorado/photo_filters/Alhambra_Filters',
          'filters': {'F_365': 'F_365_1.res',
                      'F_396': 'F_396_1.res',
                      'F_427': 'F_427_1.res',
                      'F_458': 'F_458_1.res',
                      'F_489': 'F_489_1.res',
                      'F_520': 'F_520_1.res',
                      'F_551': 'F_551_1.res',
                      'F_582': 'F_582_1.res',
                      'F_613': 'F_613_1.res',
                      'F_644': 'F_644_1.res',
                      'F_675': 'F_675_1.res',
                      'F_706': 'F_706_1.res',
                      'F_737': 'F_737_1.res',
                      'F_768': 'F_768_1.res',
                      'F_799': 'F_799_1.res',
                      'F_830': 'F_830_1.res',
                      'F_861': 'F_861_1.res',
                      'F_892': 'F_892_1.res',
                      'F_923': 'F_923_1.res',
                      'F_954': 'F_954_1.res'}}  #,
                      # 'F_H': 'F_H_1.res',
                      # 'F_J': 'F_J_1.res',
                      # 'F_KS': 'F_KS_1.res'}}

config['lambda_norm'] = 4020.

class PostParameters(object):

    def __init__(self, bt, l2m_norm_ssp):
        self.bt = bt
        self.l2m_norm_ssp = l2m_norm_ssp

    def get_parameters(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v):
        aux_p = {'a_v': a_v}
        csp_model = n_component(self.bt.ageBase)
        csp_model.add_exp(t0_young, tau_young, frac_young)
        csp_model.add_exp(t0_old, tau_old, 1 - frac_young)

        self.sfh = csp_model.get_sfh()
        aux_p['age'] = np.sum(self.sfh * self.l2m_norm_ssp * np.log10(self.bt.ageBase)) / np.sum(
            self.sfh * self.l2m_norm_ssp)

        aux_p['m2l'] = np.log10(1/np.average(self.l2m_norm_ssp, weights=self.sfh))

        return aux_p



# file_fit = h5py.File('/Users/william/Downloads/bpz_fit_full_newmask_nointerp.hdf5', 'r')
# file_fit = h5py.File('/Users/william/Downloads/bpz_fit_nointerp_newmask_chi2tx.hdf5', 'r')
file_fit = h5py.File('/Users/william/Downloads/bpz_fit_nointerp_newmask_chi2tx_CCM_test.hdf5', 'r')
params = file_fit.get('/model_parameters')
n_z, n_t, n_models, n_parameters = params.shape

file_fit_processed = h5py.File('test.hdf5', 'w')
post_parameters_data = np.empty(
    dtype=np.dtype([('age', np.float), ('metallicity', np.float), ('a_v', np.float), ('m2l', np.float)]),
    shape=(n_z, n_t, n_models))


metallicity = 0.02  # FIXME:

bt = StarlightBase(config['base_file'], config['base_path'])
bt.ageBase[bt.ageBase == 0] = 100  # Avoid NaN on log10(age)
i_norm = np.argmin((bt.l_ssp - config['lambda_norm']) ** 2)
i_met = int(np.argwhere(bt.metBase == metallicity))
l2m_norm_ssp = [bt.f_ssp[i_met, i_age, i_norm] for i_age in range(len(bt.ageBase))]

post_parameters = PostParameters(bt, l2m_norm_ssp)


for i_z in [0]:  # FIXME
    for i_t in range(n_t):
        for i_model in range(n_models):
            for k, v in post_parameters.get_parameters(*params[i_z, i_t, i_model]).iteritems():
                post_parameters_data[i_z, i_t, i_model][k] = v

file_fit_processed.create_dataset('/parameters', data=post_parameters_data, compression='gzip')
l = np.exp(file_fit.get('/model_lnprob') - np.min(file_fit.get('/model_lnprob'), axis=2)[:, :, np.newaxis])
l /= l.sum(axis=2)[:,:,np.newaxis]
file_fit_processed.create_dataset('/likelihood', data=l)

#
file_fit_processed.close()
file_fit.close()