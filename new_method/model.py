## 3 - Model
import sys
from magal.util.stellarpop import n_component
from pystarlight.util.base import StarlightBase
from astropy.cosmology import WMAP9 as cosmo
import numpy as np
from pystarlight.util.redenninglaws import Cardelli_RedLaw
from config import config
from common import load_taylor, av_taylor_coeff, mag_in_z, find_nearest_idx


class Model(object):
    def __init__(self, template_magnitudes, filters, z, base_file, base_path, taylor_file=None, magerr=1e-14,
                 is_single=False):

        # Template magnitudes
        self.template_magnitudes = template_magnitudes
        # self.chi2_constant = np.average(template_magnitudes)  # This will be summed to chi2, avoiding to get them ~ 0.
        self.filters = filters
        self.z = z
        self.z_age_limit = cosmo.age(z).to('yr').value

        # Base
        self.bt = StarlightBase(base_file, base_path)
        ## Metallicity
        aux = np.log10(self.bt.metBase)
        aux2 = (aux[1:] - aux[0:-1]) / 2
        self.met_min = aux - np.append(aux2[0], aux2)
        self.met_max = aux + np.append(aux2, aux2[-1])
        self.met_low = min(self.met_min)
        self.met_upp = max(self.met_max)

        # Extinction Law
        self.q = Cardelli_RedLaw(self.bt.l_ssp)

        # If Taylor expansion:
        if taylor_file:
            self.m, self.a, self.tz, self.int_f = load_taylor(config, filters)

        # Magnitude error
        self.magerr = magerr

    @staticmethod
    def get_sfh(bt, t0_young, tau_young, t0_old, tau_old, frac_young):
        # 1 - Eval the SFH
        csp_model = n_component(bt.ageBase)
        # If frac_young == 0, then only eval the old component and vice-versa
        if frac_young > 0:
            csp_model.add_exp(t0_young, tau_young, frac_young)
        if frac_young < 1:
            csp_model.add_exp(t0_old, tau_old, 1 - frac_young)
        if csp_model._total_M_fraction != 1.0:
            print 'ValueError: ', csp_model._total_M_fraction
            raise ValueError
        return csp_model.get_sfh()

    def get_spec(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity=0.02):
        csp_model = self.get_sfh(self.bt, t0_young, tau_young, t0_old, tau_old, frac_young)

        # 2 - Eval the correspondent spectrum
        try:
            i_met = int(np.argwhere(self.bt.metBase == metallicity))
        except TypeError:
            print 'TypeError:', self.bt.metBase, metallicity, np.argwhere(self.bt.metBase == metallicity)
            sys.exit(1)
        spec = self.bt.f_ssp[i_met] * csp_model[:, np.newaxis]  # Base spectra [??units??]
        spec *= 10 ** (-0.4 * (self.q * a_v))

        return spec.sum(axis=0)

    def mag_in_z_taylor(self, sfh, av, i_z, i_met, m, a):
        n = m.shape[4]
        return -2.5 * np.log10(
            np.sum(av_taylor_coeff(n, av, a[i_z, i_met]) * sfh[:, np.newaxis, np.newaxis] * m[i_z, i_met],
                   axis=(0, 2)) / self.int_f) - 2.41

    def get_mags(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity):
        # sfh, av, i_z, i_met, m, a, filters)
        i_z = find_nearest_idx(self.tz, self.z)
        try:
            i_met = int(np.argwhere(
                np.bitwise_and(np.log10(metallicity) >= self.met_min, np.log10(metallicity) < self.met_max)))
        except:

            print 'get_mags error:', np.argwhere(np.bitwise_and(np.log10(metallicity) >= self.met_min, np.log10(
                metallicity) < self.met_max)), self.met_min, np.log10(metallicity), self.met_max
            sys.exit(1)

        if config["chi2_spectrum"]:  # Use the pedestrian mode, calculating spectra and everything
            spec = self.get_spec(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, self.bt.metBase[i_met])
            return mag_in_z(self.bt.l_ssp, spec, self.z, self.filters)
        else:  # Shortcut by using taylor expansion on AV
            return self.mag_in_z_taylor(self.get_sfh(self.bt, t0_young, tau_young, t0_old, tau_old, frac_young), a_v,
                                        i_z, i_met, self.m, self.a)

    def lnprob(self, x):
        t0_young, tauT_young, t0_old, tauT_old, lg_frac_young, a_v, lg_metallicity = x
        if np.isfinite(self.lnprior(t0_young, tauT_young, t0_old, tauT_old, lg_frac_young, a_v, lg_metallicity)):
            ## For single burst emulation:
            if config['lg_frac_young_min'] == config['lg_frac_young_max'] == -1e30:
                lg_frac_young = -np.inf

            # # Txitxo chi^2:
            # aux_s = self.get_mags(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity)
            # return -0.5 * (np.sum(aux_s**2) - np.sum(aux_s * self.template_magnitudes)**2 / np.sum(self.template_magnitudes**2)) * (1/self.magerr)**2
            # William's magnitude chi^2:
            # \chi^2 = \sum O - M + a  ## a is the scaling factor.
            aux_s = self.template_magnitudes - self.get_mags(10 ** t0_young, 10 ** (tauT_young + t0_young),
                                                             10 ** t0_old, 10 ** (tauT_old + t0_old),
                                                             10 ** lg_frac_young, a_v, 10 ** lg_metallicity)
            return -0.5 * np.sum((aux_s - np.mean(aux_s)) ** 2) * self.magerr ** -2
        else:
            return -np.inf

    def lnprior(self, t0_young, tauT_young, t0_old, tauT_old, lg_frac_young, a_v, lg_metallicity):

        # metallicity
        if np.float64(lg_metallicity) < self.met_low or np.float64(lg_metallicity) >= self.met_upp:
            # print '1'
            return -np.inf

        # old
        if np.log10(config['t0_old_min']) > t0_old or t0_old > np.log10(self.z_age_limit):
            # print '2'
            return -np.inf

        # tau_old = tauT_old * 10 ** t0_old
        if tauT_old < np.log10(config['tau_t0_old_min']) or tauT_old > np.log10(config['tau_t0_young_max']):
            # print '3'
            return -np.inf

        # extinction
        if config['AV_min'] > a_v or a_v > config['AV_max']:
            # print '8'
            return -np.inf

        # Items below are not taken in account on single burst case!
        if config['lg_frac_young_min'] == config['lg_frac_young_max'] == -1e30:
            return 1

        # young
        if t0_young > t0_old:
            # print '4'
            return np.inf
        if np.log10(config['t0_young_min']) >= t0_young or t0_young > np.log10(config['t0_young_max']):
            # print '5'
            return -np.inf
        if tauT_young < np.log10(config['tau_t0_young_min']) or tauT_young > np.log10(config['tau_t0_young_max']):
            # print '6'
            return -np.inf


        # frac_young
        ## single burst
        if lg_frac_young > 0:
            # print '7'
            return -np.inf
        elif config['lg_frac_young_min'] > lg_frac_young or lg_frac_young > config['lg_frac_young_max']:
            return -np.inf
        # if frac_young != 0:
        #     return -np.inf

        return 1

    def get_p0(self):
        lg_t0_young = np.log10(config['t0_young_min']) + np.random.rand() * (
            np.log10(config['t0_young_max']) - np.log10(config['t0_young_min']))
        lg_tauT_young = np.log10(config['tau_t0_young_min']) + np.random.rand() * (
            np.log10(config['tau_t0_young_max']) - np.log10(config['tau_t0_young_min']))
        lg_t0_old = np.log10(config['t0_old_min']) + np.random.rand() * (
            np.log10(self.z_age_limit) - np.log10(config['t0_old_min']))
        lg_tauT_old = np.log10(config['tau_t0_old_min']) + np.random.rand() * (
            np.log10(config['tau_t0_old_max']) - np.log10(config['tau_t0_old_min']))
        if config['lg_frac_young_min'] == config['lg_frac_young_max']:
            lg_frac_young = config['lg_frac_young_max']
        else:
            lg_frac_young = config['lg_frac_young_max'] + np.log10(np.random.rand())
        a_v = config['AV_min'] + np.random.random() * (config['AV_max'] - config['AV_min'])
        lg_metallicity = self.met_low + np.random.random() * (self.met_upp - self.met_low)

        return [lg_t0_young, lg_tauT_young, lg_t0_old, lg_tauT_old, lg_frac_young, a_v, lg_metallicity]
