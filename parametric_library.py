import atpy
from pystarlight.util.base import StarlightBase

__author__ = 'william'

import os
import numpy as np
import matplotlib.pyplot as plt
from magal.library import LibraryModel
from multiprocessing.pool import Pool
from astropy import units
from astropy.cosmology import WMAP9 as cosmo
import h5py
from magal.io.readfilterset import FilterSet
from magal.photometry.syntphot import spec2filterset
from magal.util.cosmo import zcor
import pystarlight.io

from magal.fit.stats import chi2 as chi2_magal

# 0 - Open bpz templates and interpolate them
def bpz_templates():
    templates_dir = '/Users/william/workspace/pzT_templates/no_elines/'
    tpl = np.loadtxt('templates/eB11.list', np.str)
    interpolations = np.linspace(0, len(tpl) - 1, (len(tpl) - 1) * 8 + 1)
    i_interpolations = np.array(interpolations, dtype=np.int)

    templates = []

    for i in range(len(interpolations) - 1):
        i_0 = i_interpolations[i]
        i_1 = i_interpolations[i] + 1
        i_frac = interpolations[i] - i_interpolations[i]
        print i_0, i_1, interpolations[i]

        s1 = np.loadtxt('%s/%s' % (templates_dir, tpl[i_0])).T
        s2 = np.loadtxt('%s/%s' % (templates_dir, tpl[i_1])).T
        s2_ = np.interp(s1[0], s2[0], s2[1])
        s_int = s1.copy()
        s_int[1] = s1[1] * (1 - i_frac) + s2_ * i_frac
        this_template = np.empty(len(s_int[0]), dtype=np.dtype([('wl', np.float), ('flux', np.float)]))
        this_template['wl'] = s_int[0]
        this_template['flux'] = s_int[1]
        this_template['flux'] /= this_template['flux'][int(np.argwhere(this_template['wl'] >= 4020)[0])]
        templates.append(this_template)
    this_template = np.loadtxt('%s/%s' % (templates_dir, tpl[i_1]), dtype=np.dtype([('wl', np.float), ('flux', np.float)]))
    this_template['flux'] /= this_template['flux'][int(np.argwhere(this_template['wl'] >= 4020)[0])]
    templates.append(this_template)
    # np.savetxt( 'template_int_%3.2f.txt' % (interpolations[i+1]), s2 )
    return templates


# LIBRARY RANGES
def lib_guess(min_max_update=None):

    min_max = {'t0_young': (np.log10(6e6), np.log10(5e9) + 0.1),
               'theta_young': np.log10([0.001, 1000]),
               't0_old': (np.log10(1e9), np.log10(14e9) + 0.1),
               'theta_old': np.log10([0.001, 1000]),
               'frac_young': [0, 1],
               'tau_v': [0, 2]}

    aux = {'Z': [0.02]}
    # aux = {'Z': [0.00010,
    #              0.00040,
    #              0.00400,
    #              0.00800,
    #              0.02000,
    #              0.05000]}

    if min_max_update:
        min_max.update(min_max_update)

    out = {}
    for key in min_max.keys():
        out[key] = min_max[key][0] + np.random.rand() * (min_max[key][1] - min_max[key][0])

    out['tau_young'] = out['t0_young'] * 10**out.pop('theta_young')
    out['tau_old'] = out['t0_old'] * 10**out.pop('theta_old')
    out['t0_young'] = 10 ** out['t0_young']
    out['t0_old'] = 10 ** out['t0_old']

    out['Z'] = np.random.choice(aux['Z'])

    return out



def mag_in_z(spec, z, f):
    d_L = cosmo.luminosity_distance(z).to('cm')
    k_cosmo = units.Lsun.to('erg/s') / (4 * np.pi * np.power(d_L, 2))

    model_spec = zcor(spec, z)
    model_spec['flux'] *= k_cosmo
    x = spec2filterset(f.filterset, model_spec)

    return x

def chi22(l, t):  # Txitxo's chi2
    # print l.shape, t.shape
    # raw_input('debug_stop_here')
    # print 'l', np.min(l), np.max(l), np.shape(l)
    # print 't', np.min(t), np.max(t), np.shape(t)
    l = 10**(-.4*l)  # Input is in magnitudes. Convert to Flux.
    t = 10**(-.4*t)
    print 'l', np.min(l), np.max(l), np.shape(l)
    print 't', np.min(t), np.max(t), np.shape(t)
    return np.sum(l ** 2, axis=1) - np.sum(l * t, axis=1) ** 2 / np.sum(t ** 2, axis=1)

def chi2(model, lib):
    return np.array([chi2_magal(lib[z], model[z])[2] for z in range(lib.shape[0])])

def aux_chi2(m_o, m_l):
    dm = np.array(m_o) - np.array(m_l)
    return np.sum((dm - np.mean(dm)) ** 2)

def chi2(model, lib):
    return [aux_chi2(lib[z], model[z]) for z in range(lib.shape[0])]

def template_in_z(p):
    i_z, i_t, tpl = p
    # if i_z % 100 == 0:
    # if i_t % 11 == 0:
    # print 'i_t, i_z =', i_t, i_z
    return i_z, i_t, mag_in_z(tpl[i_t], z_lib[i_z], f)['m_ab']


def params(n_z, n_t, tpl):
    for i_z in range(n_z):
        for i_t in range(n_t):
            yield i_z, i_t, tpl

if __name__ == '__main__':
    # #### Defs #####
    base_path = 'Base.bc03.Padova1994.chab.All'
    base_file = '%s/BasesDir/Base.bc03.Padova1994.chab.All.hdf5' % os.environ.get('HOME')
    # base_file = '%s/Base.BC03.N' % bases_dir

    # 1 - Define library type
    # dt = np.dtype([('t0_young', np.float), ('tau_young', np.float), ('t0_old', np.float), ('tau_old', np.float),
    #                ('frac_young', np.float), ('tau_v', np.float), ('Z', np.float)])
    # lib_param = np.zeros(1, dtype=dt)
    # lib = LibraryModel('two_exp', lib_param, base_path, base_file, 'light', 4020.)

    # plots = False
    # if plots:
    # plt.clf()
    #     for i_guess in range(100):
    #         guess = lib_guess()
    #         for key in guess.keys():
    #             lib_param[key] = guess[key]
    #         print lib_param
    #
    #         spec = lib.get_model_spectrum(0)[0]
    #         # plt.clf()
    #         cut = np.bitwise_and(spec['wl'] > 3000, spec['wl'] < 9000)
    #         plt.plot(spec['wl'][cut], spec['flux'][cut])
    #         # raw_input('next...')

    f = FilterSet('/Users/william/doutorado/photo_filters/Alhambra_24.hdf5')
    f.load('Alhambra_24', 1)
    n_filters = len(f.filter_wls)

    aux_base = StarlightBase(base_file, base_path)
    base_maxage = np.max(aux_base.ageBase)
    min_z_age = cosmo.age(0.001).to('yr').value
    if base_maxage > min_z_age:
        base_maxage = min_z_age


    allow_overwrite_bpz = True
    allow_overwrite_lib = True

    bpz_lib_file = '/Users/william/tmp_out/pzT_parametric_noelines.hdf5'
    parametric_lib = '/Users/william/tmp_out/zT_weights_noelines.hdf5'
    if allow_overwrite_bpz:
        try:
            os.unlink(bpz_lib_file)
        except OSError:
            pass
    if allow_overwrite_lib:
        try:
            os.unlink(parametric_lib)
        except OSError:
            pass

    if not os.path.exists(bpz_lib_file):
        t = bpz_templates()
        n_t = len(t)

        # z_lib = np.arange(0.0001, 3.001, 0.001)
        z_lib = np.arange(0.0001, 1.501, 0.001)
        n_z = len(z_lib)

        bpz_template_library = np.empty((n_z, n_t, n_filters))

        p = params(n_z, n_t, t)

        pool = Pool()
        map_function = pool.map

        for result in map_function(template_in_z, p):
            i_z, i_t, mags = result
            bpz_template_library[i_z, i_t] = mags

        pool.close()

        f_parametric = h5py.File(bpz_lib_file)
        f_parametric.create_dataset('/templates_bpz', data=bpz_template_library)
        # f_parametric.create_dataset('/template', data=???TEMPLATE???)
        f_parametric.create_dataset('/redshift', data=z_lib)
        f_parametric.close()
        del bpz_template_library

    #### PLOT missing redshifts ####
    plot = False
    if plot:
        plt.clf()
        for i_t in range(n_t):
            print 'i_t:', i_t
            a = [np.isinf(bpz_template_library[i_z, i_t, :]).sum() for i_z in range(7001)]
            plt.plot(z_lib, a)
        plt.xlabel('z')
        plt.ylabel('Number of filters w/ no data.')
    #### PLOT missing redshifts ####

    f_parametric = h5py.File(bpz_lib_file, 'r')
    bpz_template_library = f_parametric.get('/templates_bpz')
    z_lib = f_parametric.get('/redshift')
    n_z, n_t = bpz_template_library.shape[0:2]

    max_age = cosmo.age(z_lib).to('yr').value

    dt = np.dtype([('t0_young', np.float), ('tau_young', np.float), ('t0_old', np.float), ('tau_old', np.float),
                   ('frac_young', np.float), ('tau_v', np.float), ('Z', np.float)])
    lib_param = np.zeros(1, dtype=dt)
    lib = LibraryModel('two_exp', lib_param, base_path, base_file, 'light', 4020.)

    save_model = True

    n_models = 1000
    tpl_weight = np.zeros((n_models, n_z, n_t))
    dt = np.dtype([('t0_young', np.float), ('tau_young', np.float), ('t0_old', np.float), ('tau_old', np.float),
                   ('frac_young', np.float), ('tau_v', np.float), ('Z', np.float), ('mean_age', np.float)])
    tpl_params = np.empty(n_models, dtype=dt)
    if save_model:
        model_library = np.empty((n_z, n_models, n_filters))


    # raw_input('Little stop')


    # multiprocessing = False
    multiprocessing = True

    if multiprocessing:
        pool = Pool()
        map_function = pool.map
    else:
        map_function = map

    for i_model in range(n_models):
        guess = lib_guess(min_max_update={'t0_old': (np.log10(1e9), np.log10(base_maxage - 1))})
        for k, v in guess.iteritems():
            lib_param[k] = v
            tpl_params[i_model][k] = v
        mask_age = np.bitwise_or(guess['t0_old'] > max_age, guess['t0_young'] > max_age)
        n_z_age = n_z - mask_age.sum()

        model_mags = np.empty((n_z, 1, n_filters))

        csp_tpl = [lib.get_model_spectrum(0).copy()]
        csp_tpl[0]['flux'][0] /= csp_tpl[0]['flux'][0][int(np.argwhere(csp_tpl[0][0]['wl'] >= 4020)[0])]
        a = lib._model
        print 'debug_sum> ', a.get_sfh().sum()
        mean_age = lib._csp_at_flux
        print mean_age
        # if mean_age < 10 or mean_age == np.nan:
        #     a.plot()
        #     raw_input('Error')
        print 'Running model %i' % i_model
        p = params(n_z_age, 1, csp_tpl)  # LIMIT THE TEMPLATES TO ONLY THE ONES with age < age_universe(z)
        for result in map_function(template_in_z, p):
            i_z, i_t, mags = result
            model_mags[i_z, i_t] = mags
            tpl_params[i_model]['mean_age'] = mean_age
        model_mags[mask_age] = np.nan

        for k in tpl_params.dtype.names:
            print '\t %s = %s' % (k, tpl_params[i_model][k])

        if save_model:
            model_library[:, i_model, :][:, np.newaxis, :] = model_mags

        for i_t in range(n_t):
            tpl_weight[i_model, :, i_t] = chi2(model_mags[:, 0, :], bpz_template_library[:, i_t, :])

    if multiprocessing:
        pool.close()

    f_weights = h5py.File(parametric_lib, 'w')
    f_weights.create_dataset('/zT_weights', data=np.exp(-0.5 * tpl_weight))
    f_weights.create_dataset('/tpl_params', data=tpl_params)
    f_weights.create_dataset('/redshift', data=z_lib)
    if save_model:
        f_weights.create_dataset('/tpl_magnitudes', data=model_library)
    f_weights.close()


def get_best_models(tpl_weight):
    for i in range(len(tpl_weight)):
        is_nan = np.isnan(tpl_weight[i])
        aux = tpl_weight[i][~is_nan].sum()