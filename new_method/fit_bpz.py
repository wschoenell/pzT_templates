import multiprocessing
import h5py
from magal.util.stellarpop import n_component
from pystarlight.util.base import StarlightBase
from pystarlight.util.redenninglaws import Cardelli_RedLaw
import time

__author__ = 'william'

import emcee
from astropy.cosmology import WMAP9 as cosmo
import numpy as np
import matplotlib.pyplot as plt

from common import params, mag_in_z, template_in_z, mag_in_z_taylor, load_filters, av_taylor_coeff

## 0 - Configuration
config = {'bpz_library': '../templates/eB11.list', 'n_interpolations': 7, 'bpz_library_dir': '../no_elines/',
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

pool = multiprocessing.Pool()
map_function = pool.map
plot = False


#####
z = np.arange(config['z_ini'], config['z_fin'] + config['z_delta'], config['z_delta'])

#####

## 1 - Load filter files.
filters, filter_lambdas = load_filters(config)


## 2 - Calculate bpz magnitudes over different redshifts.

#### 2.1 - Interpolate templates
templates = np.loadtxt(config['bpz_library'], np.str)
templates_data = [np.loadtxt('%s/%s' % (config['bpz_library_dir'], t), dtype=np.dtype([('lambda', np.float), ('flux', np.float)])) for t in templates]
interpolations = np.linspace(0., len(templates)-1, len(templates)+(len(templates)-1)*config['n_interpolations'])
templates_data_interpolated = []
if plot: plt.clf(); ymin = 1e10; ymax = -1e10; xlim = 3000, 9000
for interpolation in interpolations:
    if interpolation != round(interpolation):
        i_0 = int(interpolation)
        i_1 = int(interpolation + 1)
        print interpolation, i_0, i_1, interpolation - i_0, i_1 - interpolation
        aux = np.copy(templates_data[i_0])
        aux['flux'] *= interpolation - i_0
        aux['flux'] += np.interp(aux['lambda'], templates_data[i_1]['lambda'], templates_data[i_1]['flux']) * (i_1 - interpolation)
        templates_data_interpolated.append(aux)
    else:
        print interpolation
        templates_data_interpolated.append(templates_data[int(interpolation)])

    if plot:
        mask_lambda = np.bitwise_and(xlim[0] < templates_data_interpolated[-1]['lambda'], templates_data_interpolated[-1]['lambda'] < xlim[1])
        ymin = np.min(np.concatenate(([ymin], templates_data_interpolated[-1]['flux'][mask_lambda])))
        ymax = np.max(np.concatenate(([ymax], templates_data_interpolated[-1]['flux'][mask_lambda])))
        plt.plot(templates_data_interpolated[-1]['lambda'], templates_data_interpolated[-1]['flux'])

if plot:
    plt.xlim(xlim)
    plt.ylim(ymin, ymax)

#### 2.2 - Calculate magnitudes
z = np.arange(config['z_ini'], config['z_fin'] + config['z_delta'], config['z_delta'])
z_len, templates_len, magnitudes_len = len(z), len(templates_data_interpolated), len(config['filters'])
templates_magnitudes = np.empty((z_len, templates_len, magnitudes_len), dtype=np.float)  # N_z, N_t, N_mag
p = params(z_len, templates_len, z, templates_data_interpolated, filters)
for result in map_function(template_in_z, p):
    i_z, i_t, mags = result
    templates_magnitudes[i_z, i_t] = mags

pool.close()

# f = h5py.File('bpz_templates_mags.hdf5', 'w')
# f.create_dataset('/filter_lambdas', data=np.array(filter_lambdas))
# f.create_dataset('/bpz_template_magnitudes', data=templates_magnitudes)
# f.close()
# #
# sys.exit()

## 3 - Model
class Model(object):

    def __init__(self, template_magnitudes, filters, z, base_file, base_path, taylor_file=None):

        #Template magnitudes
        self.template_magnitudes = template_magnitudes
        self.chi2_constant = np.average(template_magnitudes)  # This will be summed to chi2, avoiding to get them ~ 0.
        self.filters = filters
        self.z = z
        self.z_age_limit = cosmo.age(z).to('yr').value

        # Base
        self.bt = StarlightBase(base_file, base_path)

        # Extinction Law
        self.q = Cardelli_RedLaw(self.bt.l_ssp)

        # If Taylor expansion:
        if taylor_file:
            tf = h5py.File(taylor_file)
            self.m = tf['m']
            self.a = tf['a']
            self.tz = tf['redshift']
            self.int_f = np.array([np.trapz(filters[fid]['R'] * filters[fid]['lambda']**-1) for fid in np.sort(filters.keys())])

    def get_sfh(self, t0_young, tau_young, t0_old, tau_old, frac_young):
        # 1 - Eval the SFH
        csp_model = n_component(self.bt.ageBase)
        csp_model.add_exp(t0_young, tau_young, frac_young)
        csp_model.add_exp(t0_old, tau_old, 1 - frac_young)
        return csp_model.get_sfh()

    def get_spec(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity=0.02):
        csp_model = self.get_sfh(frac_young, t0_old, t0_young, tau_old, tau_young)

        # 2 - Eval the correspondent spectrum
        i_met = int(np.argwhere(self.bt.metBase == metallicity))
        spec = self.bt.f_ssp[i_met] * csp_model[:, np.newaxis] / self.bt.Mstars[i_met][:, np.newaxis]  # Base spectra [??units??]
        spec *= 10 ** (-0.4 * (self.q * a_v))

        return spec.sum(axis=0)

    def mag_in_z_taylor(self, sfh, av, i_z, i_met, m, a):
        n = m.shape[4]
        return -2.5 * np.log10(np.sum(av_taylor_coeff(n, av, a[i_z, i_met]) * sfh[:, np.newaxis, np.newaxis] * m[i_z, i_met], axis=(0,2)) / self.int_f) - 2.41

    def get_mags(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity=0.02):
        #sfh, av, i_z, i_met, m, a, filters)
        i_z = int(np.argwhere(self.z == self.tz))
        i_met = int(np.argwhere(self.bt.metBase == metallicity))
        return self.mag_in_z_taylor(self.get_sfh(t0_young, tau_young, t0_old, tau_old, frac_young), a_v, i_z, i_met, self.m, self.a)
        #return mag_in_z(self.bt.l_ssp, self.get_spec(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity), self.z, self.filters)

    def lnprob(self, x):
        t0_young, tau_young, t0_old, tau_old, frac_young, a_v = x
        if np.isfinite(self.lnprior(t0_young, tau_young, t0_old, tau_old, frac_young, a_v)):
            # Txitxo chi^2:
            aux_s = self.get_mags(t0_young, tau_young, t0_old, tau_old, frac_young, a_v)
            return -0.5 * (np.sum(aux_s**2) - np.sum(aux_s * self.template_magnitudes)**2 / np.sum(self.template_magnitudes**2))
        else:
            return -np.inf


    def lnprior(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v):

        # old
        if 1e9 > t0_old or t0_old > self.z_age_limit:
            return -np.inf
        if tau_old < 0:
            return -np.inf

        # young
        if t0_young > self.z_age_limit:
            return -np.inf
        if 6e6 > t0_young or t0_young > 5e9:
            return -np.inf
        if 0 > frac_young or frac_young > 1:
            return -np.inf
        if tau_young < 0:
            return -np.inf

        # extinction
        if config['AV_min'] > a_v or a_v > config['AV_max']:
            return -np.inf

        return 1

    def get_p0(self):
        t0_young = np.log10(6e6) + np.random.rand() * (np.log10(5e9) - np.log10(6e6))
        tau_young = np.log10(.001) + np.random.rand() * (np.log10(.001) - np.log10(1000))
        t0_old = np.log10(1e9) + np.random.rand() * (np.log10(self.z_age_limit) - np.log10(1e9))
        tau_old = np.log10(.001) + np.random.rand() * (np.log10(.001) - np.log10(1000))
        frac_young = np.random.rand()
        a_v = 2 * np.random.rand()

        return [10**t0_young, 10**tau_young, 10**t0_old, 10**tau_old, frac_young, a_v]


# Do the modeling
print 'starting model...'

i_z = 0

test = False

if test:
    ndim, nwalkers = 6, 14
    n_samples = 200
    n_steps = 100
    n_burnin = int(n_steps*.3)
    out_fname = 'bpz_fit_full_newmask_kk_test.hdf5'
    models_fit = [0]
    print 'Running for %i templates' % len(templates_data_interpolated)
    models_fit = range(len(templates_data_interpolated))
else:
    ndim, nwalkers = 6, 100
    n_samples = 200
    n_steps = 1000
    n_burnin = int(n_steps*.3)
    out_fname = 'bpz_fit_nointerp_newmask_chi2tx_CCM_test.hdf5'
    models_fit = range(len(templates_data_interpolated))

f_fit = h5py.File(out_fname, 'w')

model_lnprob = f_fit.create_dataset('/model_lnprob', shape=(z_len, templates_len, n_samples), compression='gzip')
model_parameters = f_fit.create_dataset('/model_parameters', shape=(z_len, templates_len, n_samples, ndim), compression='gzip')
model_magnitudes = f_fit.create_dataset('/model_magnitudes', shape=(z_len, templates_len, n_samples, magnitudes_len), compression='gzip')
f_fit.create_dataset('/template_magnitudes', data=templates_magnitudes, compression='gzip')

for i_t in models_fit:  #np.array(range(0, len(templates_data_interpolated), 7))[::-1]:
    model = Model(templates_magnitudes[i_z, i_t], filters, z[i_z], config['base_file'], config['base_path'], config['taylor_file'])
    p0 = [model.get_p0() for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, model.lnprob)
    t0 = time.time()
    sampler.run_mcmc(p0, n_steps)
    print 'Total time: %3.2f secs' % (time.time() - t0)
    samples = sampler.chain[:, n_burnin:, :].reshape((-1, ndim))
    samples_lnprob = sampler.lnprobability[:, n_burnin:].reshape((-1))

    plt.clf()
    random_i = np.random.randint(len(samples), size=n_samples)
    model_lnprob[i_z, i_t] = samples_lnprob[random_i]
    model_parameters[i_z, i_t] = samples[random_i]
    for i_sample, (t0_young, tau_young, t0_old, tau_old, frac_young, a_v) in enumerate(model_parameters[i_z, i_t]):
        aux_mags = model.get_mags(t0_young, tau_young, t0_old, tau_old, frac_young, a_v)
        model_magnitudes[i_z, i_t, i_sample] = aux_mags
        plt.plot(filter_lambdas, aux_mags + np.mean(model.template_magnitudes - aux_mags), color="k", alpha=0.1)
    plt.plot(filter_lambdas, templates_magnitudes[i_z, i_t], color="r", lw=2, alpha=0.6)
    # plt.ylim(21, 25)
    # if not test:
    plt.savefig('fit_template_%i_noIR_chi2tx.png' % i_t)
    np.savez('fit_template_%i_noIR_chi2tx.npz' % i_t, samples, samples_lnprob, templates_magnitudes[i_z, i_t])
    print 'saving template %i' % i_t
    # plt.errorbar(x, y, yerr=yerr, fmt=".k")

f_fit.close()

## 99 - End