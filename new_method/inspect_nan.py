import numpy as np
import h5py

from model import Model
from config import config
from common import load_filters

file_fit = h5py.File('bpz_fit_magerr_0.0050_lib_list_random_single_z0.00_0.00_nwalk_70_nsamp_1000_nstep_3000.hdf5')

filters, filter_lambdas = load_filters(config)

isinf = np.isinf(file_fit['full_lnprob'])

bad_chains = list()
probs = list()

i = 1

# [0, chain[1]]


for chain in np.argwhere(isinf): #[:10]:
    model = Model(file_fit['template_magnitudes'][0, chain[1]], filters, config['z_ini'], config['base_file'],
                  config['base_path'], config['taylor_file'], config['magerr'])
    bad_chains.append(file_fit['full_chain'][chain[0], chain[1], chain[2], chain[3], :])
    probs.append(model.lnprob(bad_chains[-1]))
    print 'i ', i
    i += 1

bad_chains = np.array(bad_chains)
