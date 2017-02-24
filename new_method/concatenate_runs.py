import h5py

# from config import config

# def check_config(config1, config2):
#     for kw in

combine_datasets = {u'acceptance_fraction': 0,
                    u'full_chain': 0,
                    u'full_lnprob': 0,
                    u'model_lnprob': 0,
                    u'model_parameters': 0,
                    u'template_magnitudes': 0}

combine_files = ['bpz_fit_magerr_single_0.0500_lib_eB11_z0.00_0.10_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.10_0.20_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.20_0.30_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.30_0.40_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.40_0.50_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.50_0.60_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.60_0.70_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.70_0.80_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.80_0.90_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z0.90_1.00_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 'bpz_fit_magerr_single_0.0500_lib_eB11_z1.00_1.10_nwalk_70_nsamp_1000_nstep_3000.hdf5',
                 ]
f = h5py.File(combine_files[0], 'r')
f_combined = h5py.File('combined_bpz_fit_magerr_single_0.0500_lib_eB11_z0.00_1.10_nwalk_70_nsamp_1000_nstep_3000.hdf5', 'w')

n_z = f['acceptance_fraction'].shape[0]

for ds in combine_datasets:
    aux_shape = list(f[ds].shape)
    aux_shape[0] *= len(combine_files)   # combine_datasets[ds]
    aux_shape = tuple(aux_shape)
    f_combined.create_dataset(ds, shape=aux_shape, compression='gzip')
f.close()

i_file = 0
for fname in combine_files:
    print 'Combining file: %s' % fname
    f = h5py.File(fname)
    for ds in combine_datasets:
        f_combined[ds][n_z*i_file:n_z*(i_file+1)] = f[ds][:n_z]#.copy()
    i_file += 1
    f.close()

f_combined.close()