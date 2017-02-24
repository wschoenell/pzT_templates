import multiprocessing
import os
import time

import h5py
from astropy import units
from astropy.io import fits
from magal.util.stellarpop import n_component
from pystarlight.util.base import StarlightBase
from pystarlight.util.constants import d_sun, L_sun
from pystarlight.util.redenninglaws import Cardelli_RedLaw

from config import home

__author__ = 'william'
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

import math

u_Lsun = units.Lsun.to('erg/s')


def taylor_params(bt, n, z_ini, z, filters):
    for i_z in range(len(z)):
        for i_met in range(bt.nMet):
            for i_age in range(bt.nAges):
                yield i_z + z_ini, i_met, i_age, n, z[i_z], bt.l_ssp, bt.f_ssp[i_met, i_age], filters


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
    m = out.create_dataset('m', shape=(Nz, bt.nMet, bt.nAges, len(filters), n),
                           dtype=np.float32)  # , compression='gzip', compression_opts=4)
    a = out.create_dataset('a', shape=(Nz, bt.nMet, bt.nAges, len(filters)),
                           dtype=np.float32)  # , compression='gzip', compression_opts=4)
    out.create_dataset('redshift', data=z)
    nz_max = 10
    n_split = Nz / nz_max  # Limit number of redshifts per loop
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

    return np.rollaxis(
        np.array([((-0.92103403719761845 * av) ** i) / np.math.factorial(i) * 10 ** (-0.4 * av * a) for i in range(n)]),
        0, 3)


def mag_in_z_taylor(sfh, av, i_z, i_met, m, a, filters):
    n = m.shape[4]
    int_f = np.array([np.trapz(filters[fid]['R'] * filters[fid]['lambda'] ** -1) for fid in np.sort(filters.keys())])
    return -2.5 * np.log10(
        np.sum(av_taylor_coeff(n, av, a[i_z, i_met]) * sfh[:, np.newaxis, np.newaxis] * m[i_z, i_met],
               axis=(0, 2)) / int_f) - 2.41


def mag_in_z(l, f, z, filters):
    l, f = np.copy(l), np.copy(f)
    if z == 0:
        d_L = 3.0856775814671913e+19  # == 10 * units.parsec.to("cm")
    else:
        d_L = cosmo.luminosity_distance(z).to('cm')

    l *= 1 + z
    f /= (1 + z) / (units.Lsun.to('erg/s') / (4 * np.pi * np.power(d_L, 2)))

    return np.array([-2.5 * np.log10(
        np.trapz(np.interp(filters[fid]['lambda'], l, f) * filters[fid]['R'] *
                 filters[fid]['lambda'], filters[fid]['lambda']) /
        np.trapz(filters[fid]['R'] / filters[fid]['lambda'], filters[fid]['lambda']))
                     for fid in np.sort(filters.keys())], dtype=np.float) - 2.41


def params(z_len, templates_len, z_range, template_library, filters):
    for i_z in range(z_len):
        for i_t in range(templates_len):
            print i_z, i_t
            yield i_z, i_t, z_range[i_z], template_library[i_t], filters


def template_in_z(p):
    i_z, i_t, z, spec, filters = p
    return i_z, i_t, mag_in_z(spec['lambda'], spec['flux'], z, filters)


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
    spec = np.sum(bt.f_ssp[i_met] * csp_model.get_sfh()[:, np.newaxis], axis=0)  # Base spectra [??units??]
    spec *= 10 ** (-0.4 * (Cardelli_RedLaw(bt.l_ssp) * a_v))
    mags_analytic = mag_in_z(bt.l_ssp, spec, z[i_z], filters)

    # Calculate the Taylor Approximation
    for n_taylor in range(1, 10):
        t0 = time.time()
        m, a, f = taylor_matrix_ssp(bt, n_taylor, z, filters)
        mags_taylor = mag_in_z_taylor(csp_model.get_sfh(), a_v, i_z, i_met, m, a)
        print 'n_taylor, t, max |delta_m| = ', n_taylor, time.time() - t0, np.max(np.abs(mags_taylor - mags_analytic))


def test_matrix_time():
    print 'Test time to calculate all the matrix over a range of redshifts of (0, 7, 0.001)'
    bt = StarlightBase('/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'Base.bc03.Padova1994.chab.All')
    z = np.arange(0.001, 1, 0.001)
    for n_taylor in range(1, 10):
        t0 = time.time()
        aux_fname = '/Users/william/tmp_out/m_av.hdf5'
        os.unlink(aux_fname)
        m, a, f = taylor_matrix_ssp(aux_fname, bt, n_taylor, z, filters)
        print 'n_taylor, t:', n_taylor, time.time() - t0
        f.close()


# test_matrix_time()

def load_filters(config):
    filters = {}
    for fid in np.sort(config['filters'].keys()):
        filters[fid] = np.loadtxt('%s/%s' % (config['filters_dir'], config['filters'][fid]),
                                  dtype=np.dtype([('lambda', np.float), ('R', np.float)]))
        filters[fid] = filters[fid][filters[fid]['R'] > 0]

    filter_lambdas = [np.average(filters[fid]['lambda'], weights=filters[fid]['R']) for fid in
                      np.sort(config['filters'].keys())]

    return filters, filter_lambdas


def eval_matrix(config, n_taylor=10):
    filters, filter_lambdas = load_filters(config)
    print 'Calculating Taylor coefficients matrix for n_taylor = %i.' % n_taylor
    bt = StarlightBase(config['base_file'], config['base_path'])
    z = np.arange(config['z_ini'], config['z_fin'] + config['z_delta'], config['z_delta'])
    m, a, f = taylor_matrix_ssp(config['taylor_file'], bt, n_taylor, z, filters)
    config['filters'] = ','.join(sorted(config['filters'].keys()))
    f.attrs.update(config)
    f.close()


def load_taylor(config, filters):
    tf = h5py.File(config['taylor_file'], 'r')
    m = np.copy(tf['m'].value)
    a = np.copy(tf['a'].value)
    tz = np.copy(tf['redshift'].value)
    int_f = np.array(
        [np.trapz(filters[fid]['R'] * filters[fid]['lambda'] ** -1) for fid in np.sort(filters.keys())])

    # Sanity checks on Taylor file:
    try:
        for kw in ['z_ini', 'z_delta', 'base_path']:  # , 'z_fin' is not checked!
            if config[kw] != tf.attrs[kw]:
                print 'Bad taylor_file. kw = %s. val = %s' % (kw, tf.attrs[kw])
                # raise Exception
                # raw_input('Press ENTER to ignore or CTRL+C to abort.')
        for kw in ['base_file']:
            if os.path.basename(config[kw]) != os.path.basename(tf.attrs[kw]):
                print 'Bad taylor_file. kw = %s. val = %s' % (kw, tf.attrs[kw])
                raise Exception
        if tf.attrs['filters'] != ','.join(sorted(config['filters'].keys())):
            print 'Bad taylor_file. (filters)'
            raise Exception
    except KeyError:
        print('KeyError: %s' % kw)
        raise KeyError

    tf.close()
    return m, a, tz, int_f


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


def calc_LsunFilter(filtercurve):
    '''
        Calculates the solar luminosity on self.filter in the same way as Bruzual does in Galaxev.
                NVA helped on this!
    '''
    # The solar spectra is located on the data directory which comes with the pystarlight distribution:
    data_dir = '%s/workspace/pystarlight/data/solar_spectrum/' % home  # os.path.dirname(inspect.getfile(pystarlight)) + '/../../data/solar_spectrum/'
    sun_spec, header_sun_spec = fits.getdata(data_dir + 'sun_reference_stis_002.fits',
                                             header=True)  # Read Solar Spectrum
    Lsun_corr = np.interp(filtercurve['lambda'], sun_spec['WAVELENGTH'], sun_spec['FLUX'])  # Correct to filter lambdas
    LsunFilter = np.trapz(Lsun_corr * filtercurve['transm'],filtercurve['lambda'])
    # / np.trapz(filtercurve['transm'], filtercurve['lambda']) # Calc solar luminosity (ergs/s/cm2)
    LsunFilter = (LsunFilter * 4 * np.pi * d_sun ** 2) / L_sun  # Convert: erg/s/cm2 --> L_sun

    return LsunFilter


def find_nearest_idx(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
        return idx - 1
    else:
        return idx
