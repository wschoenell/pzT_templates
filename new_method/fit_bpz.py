import multiprocessing
import h5py
from magal.util.stellarpop import n_component
from pystarlight.util.base import StarlightBase
from pystarlight.util.redenninglaws import Cardelli_RedLaw
import time
import sys
import emcee
from astropy.cosmology import WMAP9 as cosmo
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import params, template_in_z, load_filters, av_taylor_coeff
## 0 - Configuration
from config import config

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
templates_data = [
    np.loadtxt('%s/%s' % (config['bpz_library_dir'], t), dtype=np.dtype([('lambda', np.float), ('flux', np.float)])) for
    t in templates]
interpolations = np.linspace(0., len(templates) - 1, len(templates) + (len(templates) - 1) * config['n_interpolations'])
templates_data_interpolated = []
if plot: plt.clf(); ymin = 1e10; ymax = -1e10; xlim = 3000, 9000
for interpolation in interpolations:
    if interpolation != round(interpolation):
        i_0 = int(interpolation)
        i_1 = int(interpolation + 1)
        print interpolation, i_0, i_1, interpolation - i_0, i_1 - interpolation
        aux = np.copy(templates_data[i_0])
        aux['flux'] *= interpolation - i_0
        aux['flux'] += np.interp(aux['lambda'], templates_data[i_1]['lambda'], templates_data[i_1]['flux']) * (
            i_1 - interpolation)
        templates_data_interpolated.append(aux)
    else:
        print interpolation
        templates_data_interpolated.append(templates_data[int(interpolation)])

    if plot:
        mask_lambda = np.bitwise_and(xlim[0] < templates_data_interpolated[-1]['lambda'],
                                     templates_data_interpolated[-1]['lambda'] < xlim[1])
        ymin = np.min(np.concatenate(([ymin], templates_data_interpolated[-1]['flux'][mask_lambda])))
        ymax = np.max(np.concatenate(([ymax], templates_data_interpolated[-1]['flux'][mask_lambda])))
        if interpolation == round(interpolation):
            plt.plot(templates_data_interpolated[-1]['lambda'], templates_data_interpolated[-1]['flux'], lw=2)
        else:
            plt.plot(templates_data_interpolated[-1]['lambda'], templates_data_interpolated[-1]['flux'])

if plot:
    plt.xlim(xlim)
    plt.ylim(ymin, ymax)
    plt.savefig('spec.png' % interpolation)
    sys.exit()

#### 2.2 - Calculate magnitudes
z = np.arange(config['z_ini'], config['z_fin'] + config['z_delta'], config['z_delta'])
z_len, templates_len, magnitudes_len = len(z), len(templates_data_interpolated), len(config['filters'])
templates_magnitudes = np.empty((z_len, templates_len, magnitudes_len), dtype=np.float)  # N_z, N_t, N_mag
p = params(z_len, templates_len, z, templates_data_interpolated, filters)
for result in map_function(template_in_z, p):
    i_z, i_t, mags = result
    templates_magnitudes[i_z, i_t] = mags

pool.close()


## 3 - Model
class Model(object):
    def __init__(self, template_magnitudes, filters, z, base_file, base_path, taylor_file=None, magerr=1e-14):

        # Template magnitudes
        self.template_magnitudes = template_magnitudes
        self.chi2_constant = np.average(template_magnitudes)  # This will be summed to chi2, avoiding to get them ~ 0.
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
            tf = h5py.File(taylor_file, 'r')
            self.m = np.copy(tf['m'].value)
            self.a = np.copy(tf['a'].value)
            self.tz = np.copy(tf['redshift'].value)
            self.int_f = np.array(
                [np.trapz(filters[fid]['R'] * filters[fid]['lambda'] ** -1) for fid in np.sort(filters.keys())])
            tf.close()

        # Magnitude error
        self.magerr = magerr

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
        i_z = int(np.argwhere(self.z == self.tz))
        try:
            i_met = int(np.argwhere(np.bitwise_and(metallicity >= self.met_min, metallicity < self.met_max)))
        except:
            print np.argwhere(np.bitwise_and(metallicity >= self.met_min,
                                             metallicity < self.met_max)), self.met_min, metallicity, self.met_max
            sys.exit(1)
        return self.mag_in_z_taylor(self.get_sfh(t0_young, tau_young, t0_old, tau_old, frac_young), a_v, i_z, i_met,
                                    self.m, self.a)
        # return mag_in_z(self.bt.l_ssp, self.get_spec(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity), self.z, self.filters)

    def lnprob(self, x):
        t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity = x
        if np.isfinite(self.lnprior(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity)):
            # # Txitxo chi^2:
            # aux_s = self.get_mags(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity)
            # return -0.5 * (np.sum(aux_s**2) - np.sum(aux_s * self.template_magnitudes)**2 / np.sum(self.template_magnitudes**2)) * (1/self.magerr)**2
            # William's magnitude chi^2:
            # \chi^2 = \sum O - M + a  ## a is the scaling factor.
            aux_s = self.template_magnitudes - self.get_mags(t0_young, tau_young, t0_old, tau_old, frac_young, a_v,
                                                             metallicity)
            return -0.5 * np.sum((aux_s - np.mean(aux_s)) ** 2) * self.magerr ** -2
        else:
            return -np.inf

    def lnprior(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity):

        # metallicity
        if np.float64(metallicity) < self.met_low or np.float64(metallicity) >= self.met_upp:
            return -np.inf

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
        a_v = config['AV_min'] + np.random.random() * (config['AV_max'] - config['AV_min'])  # 2 * np.random.rand()
        metallicity = self.met_low + np.random.random() * (self.met_upp - self.met_low)

        return [10 ** t0_young, 10 ** tau_young, 10 ** t0_old, 10 ** tau_old, frac_young, a_v, metallicity]


# Do the modeling
print 'starting model...'

test = False

if test:
    ndim, nwalkers = 7, 14
    n_samples = 200
    n_steps = 100
    n_burnin = int(n_steps * .3)
    out_fname = 'bpz_fit_full_newmask_kk_test.hdf5'
    models_fit = [0]
    z_fit = [0, 1, 2]
    print 'Running for %i templates' % len(templates_data_interpolated)
    models_fit = range(len(templates_data_interpolated))
else:
    ndim, nwalkers = 7, config['n_walkers']
    n_samples = config['n_samples']
    n_steps = config['n_steps']
    n_burnin = int(n_steps * .3)
    out_fname = config['fit_bpz_outfile']
    models_fit = range(len(templates_data_interpolated))
    z_fit = range(z_len)

f_fit = h5py.File(out_fname, 'w')

model_lnprob = f_fit.create_dataset('/model_lnprob', shape=(z_len, templates_len, n_samples), compression='gzip')
model_parameters = f_fit.create_dataset('/model_parameters', shape=(z_len, templates_len, n_samples, ndim),
                                        compression='gzip')
f_fit.create_dataset('/template_magnitudes', data=templates_magnitudes, compression='gzip')


def params_fit(models_fit, z_fit):
    for i_t in models_fit:
        for i_z in z_fit:
            yield i_t, i_z


def run_fit(p):
    i_t, i_z = p
    model = Model(templates_magnitudes[i_z, i_t], filters, z[i_z], config['base_file'], config['base_path'],
                  config['taylor_file'], config['magerr'])
    p0 = [model.get_p0() for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, model.lnprob)
    t0 = time.time()
    sampler.run_mcmc(p0, n_steps)
    print 'i_t, i_z, time: %i, %i, %3.2f secs' % (i_t, i_z, time.time() - t0)
    samples = sampler.chain[:, n_burnin:, :].reshape((-1, ndim))
    samples_lnprob = sampler.lnprobability[:, n_burnin:].reshape((-1))

    # plt.clf()
    random_i = np.random.randint(len(samples), size=n_samples)

    return i_t, i_z, samples_lnprob[random_i], samples[random_i]


pool = multiprocessing.Pool()
plot = False

print 'N_models, N_z = ', len(models_fit), len(z_fit)
fit_pars = params_fit(models_fit, z_fit)
for i_t, i_z, lnprob, parameters in pool.map(run_fit, fit_pars):
    model_lnprob[i_z, i_t] = lnprob
    model_parameters[i_z, i_t] = parameters

# Save the magnitudes.
model_magnitudes = np.empty((z_len, templates_len, n_samples, magnitudes_len))
for i_z in z_fit:
    for i_t in range(templates_len):
        model = Model(templates_magnitudes[i_z, i_t], filters, z[i_z], config['base_file'], config['base_path'],
                      config['taylor_file'], config['magerr'])
        for i_sample, (t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity) in enumerate(
                model_parameters[i_z, i_t]):
            print i_sample, i_z, i_t
            aux_mags = model.get_mags(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity)
            model_magnitudes[i_z, i_t, i_sample] = aux_mags

f_fit.create_dataset('/model_magnitudes', data=model_magnitudes, compression='gzip')
# # End of saving magnitudes #


f_fit.close()

pool.close()

## 99 - End
