import multiprocessing
import os
from astropy import units
import emcee
import h5py
from magal.util.stellarpop import n_component
from pystarlight.util.base import StarlightBase
from pystarlight.util.redenninglaws import Cardelli_RedLaw
import time


__author__ = 'william'
import matplotlib.pyplot as plt
import numpy as np

from astropy.cosmology import WMAP9 as cosmo

u_Lsun = units.Lsun.to('erg/s')

def taylor_params(bt, n, z_ini, z, filters):
    for i_z in range(len(z)):
        for i_met in range(bt.nMet):
            for i_age in range(bt.nAges):
                yield i_z+z_ini, i_met, i_age, n, z[i_z], bt.l_ssp, bt.f_ssp[i_met, i_age], filters

def m_factor(p):
    i_z, i_met, i_age, n, z, l_ssp, f_ssp, filters = p
    d_L = cosmo.luminosity_distance(z).to('cm')
    l = np.copy(l_ssp)
    q_lambda = Cardelli_RedLaw(l)  # Extinction occurs at z = 0.
    l *= (1 + z)
    f = np.copy(f_ssp)
    f *= np.copy(u_Lsun) / (4 * np.pi * np.power(d_L, 2) * (1 + z))
    out_m = []
    out_a = []
    for filter_id in np.sort(filters.keys()):
        aux_lambda = filters[filter_id]['lambda']
        aux_q = np.interp(aux_lambda, l, q_lambda)
        a = np.interp(np.average(aux_lambda), l, q_lambda)
        out_m.append(
            [np.trapz(np.interp(aux_lambda, l, f) * ((aux_q - a) ** i_taylor) * aux_lambda * filters[filter_id]['R'])
             for i_taylor in range(n)])
        out_a.append(a)
    return i_z, i_met, i_age, out_m, out_a
    # f = np.interp(filt['lambda'], l, )
    # fid = 'F_892'
    # return np.trapz(f * (Cardelli_RedLaw(filters[fid]['lambda']) - a)**i_taylor * filters[fid]['lambda'] * filters[fid]['R'])

def taylor_matrix_ssp(outfile, bt, n, z, filters, multiproc=True):
    if multiproc:
        pool = multiprocessing.Pool()
        map_function = pool.map
    else:
        map_function = map

    out = h5py.File(outfile, 'w')

    Nz = len(z)
    m = out.create_dataset('m', shape=(Nz, bt.nMet, bt.nAges, len(filters), n), dtype=np.float32)#, compression='gzip', compression_opts=4)
    a = out.create_dataset('a', shape=(Nz, bt.nMet, bt.nAges, len(filters)), dtype=np.float32)#, compression='gzip', compression_opts=4)
    out.create_dataset('redshift', data=z)
    nz_max = 10
    n_split = Nz/nz_max  # Limit number of redshifts per loop
    z_ini = 0
    t0 = time.time()
    for aux_z in np.array_split(z, n_split):
        print 't = ', time.time() - t0
        print 'Working on z_min, z_max = ', np.min(aux_z), np.max(aux_z)
        p = taylor_params(bt, n, z_ini, aux_z, filters)
        z_ini += len(aux_z)
        for result in map_function(m_factor, p):
            i_z, i_met, i_age, out_m, out_a = result
            m[i_z, i_met, i_age] = out_m
            a[i_z, i_met, i_age] = out_a

    if multiproc:
        pool.close()

    return m, a, out

def av_taylor_coeff(n, av, a):
    '''
    Returns the coefficients to the Taylor expansion on a of the reddening term:
    :math:`r_\lambda = 10^{-0.4 A_V q_\lambda}`

    The Taylor expansion is on the form:
    :math:`r_\lambda = \sum_{n=0}^\infty \frac{(-0.4 A_V \ln 10)^n}{n!} 10^{-0.4 A_V a} (q_\lambda - A_V)^n`
    '''

    return np.rollaxis(np.array([((-0.92103403719761845 * av)**i)/np.math.factorial(i) * 10**(-0.4 * av * a) for i in range(n)]), 0, 3)


def mag_in_z_taylor(sfh, av, i_z, i_met, m, a):
    n = m.shape[4]
    int_f = np.array([np.trapz(filters[fid]['R'] * filters[fid]['lambda']**-1) for fid in np.sort(filters.keys())])
    return -2.5 * np.log10(np.sum(av_taylor_coeff(n, av, a[i_z, i_met]) * sfh[:, np.newaxis, np.newaxis] * m[i_z, i_met], axis=(0,2)) / int_f) - 2.41

def mag_in_z(l, f, z, filters):
    l, f = np.copy(l), np.copy(f)
    d_L = cosmo.luminosity_distance(z).to('cm')

    l *= 1 + z
    f /= (1 + z) / (units.Lsun.to('erg/s') / (4 * np.pi * np.power(d_L, 2)))

    return np.array([-2.5 * np.log10(
                    np.trapz(np.interp(filters[fid]['lambda'], l, f) * filters[fid]['R'] *
                            filters[fid]['lambda'], filters[fid]['lambda']) /
                    np.trapz(filters[fid]['R'] / filters[fid]['lambda'], filters[fid]['lambda']) )
    for fid in np.sort(filters.keys())], dtype=np.float) - 2.41

def params(z_len, templates_len, z_range, template_library, filters):
    for i_z in range(z_len):
        for i_t in range(templates_len):
            print i_z, i_t
            yield i_z, i_t, z_range[i_z], template_library[i_t], filters

def template_in_z(p):
    i_z, i_t, z, spec, filters = p
    return i_z, i_t, mag_in_z(spec['lambda'], spec['flux'], z, filters)

## 0 - Configuration
config = {'bpz_library': '../templates/eB11.list', 'n_interpolations': 0, 'bpz_library_dir': '../no_elines/',
          'z_ini': 1e-4, 'z_fin': 7.0, 'z_delta': 0.001,
          'base_file': '/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'base_path': 'Base.bc03.Padova1994.chab.All',
          'AV_min': 0, 'AV_max': 2,
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

# pool = multiprocessing.Pool()
map_function = map
plot = False


#####
z = np.arange(config['z_ini'], config['z_fin'] + config['z_delta'], config['z_delta'])

#####

## 1 - Load filter files.
filters = {}
for fid in np.sort(config['filters'].keys()):
    filters[fid] = np.loadtxt('%s/%s' % (config['filters_dir'], config['filters'][fid]), dtype=np.dtype([('lambda', np.float), ('R', np.float)]))

# filter_lambdas = {}
# for fid in np.sort(config['filters'].keys()):
#     filter_lambdas[fid] = np.average(filters[fid]['lambda'], weights=filters[fid]['R'])
#
filter_lambdas = [np.average(filters[fid]['lambda'], weights=filters[fid]['R']) for fid in np.sort(config['filters'].keys())]

## 2 - Calculate bpz magnitudes over different redshifts.

def do_not_run():
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

    # pool.close()

    # f = h5py.File('bpz_templates_mags.hdf5', 'w')
    # f.create_dataset('/filter_lambdas', data=np.array(filter_lambdas))
    # f.create_dataset('/bpz_template_magnitudes', data=templates_magnitudes)
    # f.close()
    # #
    # sys.exit()

    ## 3 - Model
    class Model(object):

        def __init__(self, template_magnitudes, filters, z, base_file, base_path):

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

        def get_spec(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity=0.02):
            # 1 - Eval the SFH
            csp_model = n_component(self.bt.ageBase)
            csp_model.add_exp(t0_young, tau_young, frac_young)
            csp_model.add_exp(t0_old, tau_old, 1 - frac_young)

            # 2 - Eval the correspondent spectrum
            i_met = int(np.argwhere(self.bt.metBase == metallicity))
            spec = self.bt.f_ssp[i_met] * csp_model.get_sfh()[:, np.newaxis] / self.bt.Mstars[i_met][:, np.newaxis]  # Base spectra [??units??]
            spec *= 10 ** (-0.4 * (self.q * a_v))

            return spec.sum(axis=0)

        def get_mags(self, t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity=0.02):
            return mag_in_z(self.bt.l_ssp, self.get_spec(t0_young, tau_young, t0_old, tau_old, frac_young, a_v, metallicity), self.z, self.filters)

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

    test = True

    if test:
        ndim, nwalkers = 6, 14
        n_samples = 200
        n_steps = 10
        n_burnin = int(n_steps*.3)
        out_fname = 'bpz_fit_full_newmask_kk.hdf5'
        models_fit = [0]
    else:
        ndim, nwalkers = 6, 100
        n_samples = 200
        n_steps = 1000
        n_burnin = int(n_steps*.3)
        out_fname = 'bpz_fit_nointerp_newmask_chi2tx_CCM.hdf5'
        models_fit = range(len(templates_data_interpolated))

    f_fit = h5py.File(out_fname, 'w')

    model_lnprob = f_fit.create_dataset('/model_lnprob', shape=(z_len, templates_len, n_samples))
    model_parameters = f_fit.create_dataset('/model_parameters', shape=(z_len, templates_len, n_samples, ndim))
    model_magnitudes = f_fit.create_dataset('/model_magnitudes', shape=(z_len, templates_len, n_samples, magnitudes_len))

    for i_t in models_fit:  #np.array(range(0, len(templates_data_interpolated), 7))[::-1]:
        model = Model(templates_magnitudes[i_z, i_t], filters, z[0], config['base_file'], config['base_path'])
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
        plt.ylim(21, 25)
        # if not test:
        plt.savefig('fit_template_%i_noIR_chi2tx.png' % i_t)
        np.savez('fit_template_%i_noIR_chi2tx.npz' % i_t, samples, samples_lnprob)
        print 'saving template %i' % i_t
        # plt.errorbar(x, y, yerr=yerr, fmt=".k")

    f_fit.close()

    ## 99 - End


    #
    # samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

    upp_av, low_av = config['AV_max'], config['AV_min']
    x0 = (upp_av + low_av)/2


def test_magnitudes():
    bt = StarlightBase('/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'Base.bc03.Padova1994.chab.All')
    i_met = 0
    i_z = 0

    z = np.array([0.5, 1.0, 1.5])

    # Generate random parameters to test:
    z_age_limit = cosmo.age(0.01).to('yr').value
    t0_young = np.log10(6e6) + np.random.rand() * (np.log10(5e9) - np.log10(6e6))
    tau_young = np.log10(.001) + np.random.rand() * (np.log10(.001) - np.log10(1000))
    t0_old = np.log10(1e9) + np.random.rand() * (np.log10(z_age_limit) - np.log10(1e9))
    tau_old = np.log10(.001) + np.random.rand() * (np.log10(.001) - np.log10(1000))
    frac_young = np.random.rand()
    a_v = 2 * np.random.rand()

    print 'Parameters: ', z_age_limit, t0_young, tau_young, t0_old, tau_old, frac_young, a_v

    # Create the CSP with it
    csp_model = n_component(bt.ageBase)
    csp_model.add_exp(t0_young, tau_young, frac_young)
    csp_model.add_exp(t0_old, tau_old, 1 - frac_young)

    # Calculate the analytic integral
    spec = np.sum(bt.f_ssp[i_met] * csp_model.get_sfh()[:, np.newaxis] / bt.Mstars[i_met][:, np.newaxis], axis=0)  # Base spectra [??units??]
    spec *= 10 ** (-0.4 * (Cardelli_RedLaw(bt.l_ssp) * a_v))
    mags_analytic = mag_in_z(bt.l_ssp, spec, z[i_z], filters)

    # Calculate the Taylor Approximation
    for n_taylor in range(1,10):
        t0 = time.time()
        m, a, f = taylor_matrix_ssp(bt, n_taylor, z, filters)
        mags_taylor = mag_in_z_taylor(csp_model.get_sfh() / bt.Mstars[i_met], a_v, i_z, i_met, m, a)
        print 'n_taylor, t, max |delta_m| = ', n_taylor, time.time() - t0, np.max(np.abs(mags_taylor - mags_analytic))

def test_matrix_time():
    print 'Test time to calculate all the matrix over a range of redshifts of (0, 7, 0.001)'
    bt = StarlightBase('/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'Base.bc03.Padova1994.chab.All')
    z = np.arange(0.001, 1, 0.001)
    for n_taylor in range(1,10):
        t0 = time.time()
        aux_fname = '/Users/william/tmp_out/m_av.hdf5'
        os.unlink(aux_fname)
        m, a, f = taylor_matrix_ssp(aux_fname, bt, n_taylor, z, filters)
        print 'n_taylor, t:', n_taylor, time.time() - t0
        f.close()

# test_matrix_time()

def eval_matrix(n_taylor=6):
    print 'Calculating Taylor coefficients matrix for n_taylor = 5.'
    bt = StarlightBase('/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'Base.bc03.Padova1994.chab.All')
    z = np.arange(config['z_ini'], config['z_fin']+config['z_delta'], config['z_delta'])
    aux_fname = '/Users/william/tmp_out/m_a.hdf5'
    m, a, f = taylor_matrix_ssp(aux_fname, bt, n_taylor, z, filters)
    f.close()

eval_matrix()


# # Indefinite integral:
# def integ_analitic(low, up):
#     return -1.08574 * np.exp(-.921034 * up) +1.08574 * np.exp(-.921034 * low)
#
# def d_func(x, n):
#     # d/dx a^(-k x) = a^(-k x) * ln(a) * -k
#     return 10**(-.4 * x) * (-0.4 * 2.3025850929940459)**n
#
# def aux_taylor(x0, n, x):
#     print [d_func(x0, i)/np.math.factorial(i) for i in range(n)], i
#     return np.sum([(d_func(x0, i)/np.math.factorial(i)) * (x - x0)**i for i in range(n)])
#
# def integ_taylor(x0, n, low, up):
#     return aux_taylor(x0, n, up) - aux_taylor(x0, n, low)

