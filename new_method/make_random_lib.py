from pystarlight.util.base import StarlightBase

from model import Model
from config import config
import numpy as np

lib_size = 500
lib_dir = './random_library_single_burst_nolimit/'

m = Model(0, 0, 0, config['base_file'], config['base_path'])
m.bt.ageBase[0] = 1e5

params = []
for i in range(lib_size):
    while True:
        p0 = m.get_p0()
        lg_metallicity = p0[6]
        metallicity = m.bt.metBase[
            int(np.argwhere(np.bitwise_and(lg_metallicity >= m.met_min, lg_metallicity < m.met_max)))]
        # get_spec: t0_young, tau_young, t0_old, tau_old, frac_young, a_v, met
        # get_p0: t0_young, tauT_young, t0_old, tauT_old, lg_frac_young, a_v, lg_metallicity
        # lnprior: t0_young, tauT_young, t0_old, tauT_old, lg_frac_young, a_v, lg_metallicity
        if np.isfinite(m.lnprior(*p0)):
            p0[1] = 10 ** (p0[0] + p0[1])  # tau_young
            p0[0] = 10 ** p0[0]  # t0_young
            p0[3] = 10 ** (p0[2] + p0[3])  # tau_old
            p0[2] = 10 ** p0[2]  # t0_old
            p0[4] = 10 ** p0[4]  # frac_young
            p0[6] = metallicity

            # check if at_mass > 8.5
            # sfh = m.get_sfh(m.bt, *p0[0:5])
            # at_mass = np.sum(sfh * np.log10(m.bt.ageBase)) / np.sum(sfh)
            # if at_mass > 8.5:
            #     spec = np.transpose([m.bt.l_ssp, m.get_spec(*p0)])
            #     np.savetxt(lib_dir + '/random_%i.txt' % i, spec)
            #     p0.insert(0, 'random_%i.txt' % i)
            #     params.append(p0)
            #     print 'Added: ', p0
            #     break
            # else:
            #     print 'at_mass < 8.5. Skipping...'
        else:
            print 'Bad params: ', i, str(p0)
            print m.met_min, m.met_max
            # raw_input()

np.savetxt(lib_dir + '/list_random_single_nolimit.txt', params, fmt='%s')
