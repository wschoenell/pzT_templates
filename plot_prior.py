import atpy
from magal.util.stellarpop import n_component

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

# #### Defs #####
bases_dir = '%s/BasesDir/' % os.environ.get('HOME')
# base_file = '%s/Base.bc03.Padova1994.chab.All' % bases_dir
base_file = '%s/Base.BC03.N' % bases_dir

bt = atpy.Table(base_file, '/Users/william/BasesDir/', type='starlightv4_base', read_basedir=True)
base_ages = bt['age_base']
base_maxage = np.max(base_ages)
min_z_age = cosmo.age(0.001).to('yr').value
if base_maxage > min_z_age:
    base_maxage = min_z_age

from parametric_library import lib_guess

Z_model = 0.02  # Z_\odot

flag_Z = (bt['Z_base'] == bt['Z_base'][np.argmin((Z_model - bt['Z_base']) ** 2)])  # Get a mask to the desired metallicity


n_model = 10000

_i_norm = int(np.argwhere(bt.l_ssp[0] == 4020.))

at_flux = []
theta_young = []

for i_model in range(n_model):
    mod = n_component(base_ages[flag_Z])
    guess = lib_guess(min_max_update={'t0_old': (np.log10(1e9), np.log10(base_maxage - 1))})
    mod.add_exp(guess['t0_old'], guess['tau_old'], 1 - guess['frac_young'])
    mod.add_exp(guess['t0_young'], guess['tau_young'], guess['frac_young'])
    sfh = mod.get_sfh()

    _ages = mod.ages_end - (mod.ages_end - mod.ages_start) / 2
    l2m_norm_ssp = [bt[flag_Z][i]['f_ssp'][_i_norm] for i in range(len(_ages))]
    csp_at_flux = np.sum(sfh * l2m_norm_ssp * np.log10(_ages)) / np.sum(sfh * l2m_norm_ssp)
    at_flux.append(csp_at_flux)
    theta_young.append(guess['tau_young']/guess['t0_young'])

plt.clf()
plt.hist(theta_young, bins=50)
plt.clf()
plt.hist(at_flux, bins=50)