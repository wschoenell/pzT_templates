import multiprocessing
import random
import sys
import time
import emcee
import h5py
import matplotlib
import numpy as np
from model import Model

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import params, template_in_z, load_filters
## 0 - Configuration
from config import config, save_dir

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
templates = np.loadtxt(config['bpz_library'], np.str, usecols=(0,))
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
z_len, templates_len, magnitudes_len = len(z), len(templates_data_interpolated), len(config['filters'])
templates_magnitudes = np.empty((z_len, templates_len, magnitudes_len), dtype=np.float)  # N_z, N_t, N_mag
p = params(z_len, templates_len, z, templates_data_interpolated, filters)
for result in map_function(template_in_z, p):
    i_z, i_t, mags = result
    if config['mag_noise'] == 0:
        templates_magnitudes[i_z, i_t] = mags
    else:
        templates_magnitudes[i_z, i_t] = mags + np.random.normal(scale=config['mag_noise'], size=len(mags))

pool.close()

# Do the modeling
print 'starting model...'

test = False

if test:
    ndim, nwalkers = 7, 14
    n_samples = 10
    n_steps = 100
    n_burnin = int(n_steps * .3)
    out_fname = 'bpz_fit_full_newmask_kk_test.hdf5'
    models_fit = [0]
    z_fit = range(z_len)
    print 'Running for %i templates' % len(templates_data_interpolated)
    models_fit = range(len(templates_data_interpolated))
else:
    ndim, nwalkers = 7, config['n_walkers']
    n_samples = config['n_samples']
    n_steps = config['n_steps']
    n_burnin = int(n_steps * .3)
    out_fname = '%s/%s' % (save_dir, config['fit_bpz_outfile'])
    models_fit = range(len(templates_data_interpolated))
    z_fit = range(z_len)

f_fit = h5py.File(out_fname, 'w')
# Store config used on the run as attributes
aux_config = config.copy()
aux_config['filters'] = ','.join(sorted(aux_config['filters'].keys()))
f_fit.attrs.update(aux_config)


model_lnprob = f_fit.create_dataset('/model_lnprob', shape=(z_len, templates_len, n_samples), compression='gzip')
model_parameters = f_fit.create_dataset('/model_parameters', shape=(z_len, templates_len, n_samples, ndim),
                                        compression='gzip')
acceptance_fraction = f_fit.create_dataset('/acceptance_fraction', shape=(z_len, templates_len, nwalkers),
                                           compression='gzip')
full_chain = f_fit.create_dataset('/full_chain', shape=(z_len, templates_len, config['save_chains'], n_steps, 7),
                                  compression='gzip')
full_lnprob = f_fit.create_dataset('/full_lnprob', shape=(z_len, templates_len, config['save_chains'], n_steps),
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
    print 'i_t, i_z, acceptance_fraction, time: %i, %i, %.2f, %3.2f secs' % (
        i_t, i_z, np.average(sampler.acceptance_fraction), time.time() - t0)
    samples = sampler.chain[:, n_burnin:, :].reshape((-1, ndim))
    samples_lnprob = sampler.lnprobability[:, n_burnin:].reshape((-1))

    # plt.clf()
    # random_i = np.random.randint(len(samples), size=n_samples)
    random_i = list(np.sort(random.sample(range(len(samples)), n_samples)))
    random_chains = random.sample(range(nwalkers), config['save_chains'])

    return i_t, i_z, samples_lnprob[random_i], samples[random_i], sampler.acceptance_fraction, \
           sampler.chain[random_chains], sampler.lnprobability[random_chains]
    # np.random.random_integers(0, nwalkers - 1, config['save_chains'])]


pool = multiprocessing.Pool()
plot = False

print 'N_models, N_z = ', len(models_fit), len(z_fit)
fit_pars = params_fit(models_fit, z_fit)
for i_t, i_z, lnprob, parameters, acceptance_fraction[i_z, i_t], full_chain[i_z, i_t], full_lnprob[
    i_z, i_t] in pool.map(run_fit, fit_pars):
    model_lnprob[i_z, i_t] = lnprob
    model_parameters[i_z, i_t] = parameters

# Save the magnitudes. #
model_magnitudes = np.empty((z_len, templates_len, n_samples, magnitudes_len))
for i_z in z_fit:
    print 'Evaluating magnitudes for z = %.2f' % z_fit[i_z]
    for i_t in range(templates_len):
        model = Model(templates_magnitudes[i_z, i_t], filters, z[i_z], config['base_file'], config['base_path'],
                      config['taylor_file'], config['magerr'])
        for i_sample, (
        lg_t0_young, lg_tauT_young, lg_t0_old, lg_tauT_old, lg_frac_young, a_v, lg_metallicity) in enumerate(
                model_parameters[i_z, i_t]):
            # print i_sample, i_z, i_t
            aux_mags = model.get_mags(10 ** lg_t0_young, 10 ** (lg_tauT_young + lg_t0_young),
                                      10 ** lg_t0_old, 10 ** (lg_tauT_old + lg_t0_old),
                                      10 ** lg_frac_young, a_v, 10 ** lg_metallicity)
            model_magnitudes[i_z, i_t, i_sample] = aux_mags

f_fit.create_dataset('/model_magnitudes', data=model_magnitudes, compression='gzip')
# End of saving magnitudes #


f_fit.close()

pool.close()

## 99 - End
