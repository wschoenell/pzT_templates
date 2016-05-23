import os

home = os.path.expanduser('~/')
save_dir = '%s/data/out/' % home
config = {
    'bpz_library': '%s/workspace/pzT_templates/templates/Base.BC03.N.list' % home, 'n_interpolations': 1,
    'bpz_library_dir': '/Users/william/data/BasesDir/',
    'z_ini': 1e-3, 'z_fin': 1e-3, 'z_delta': 0.05,
    # 'bpz_library': '%s/workspace/pzT_templates/templates/eB11.list' % home, 'n_interpolations': 3,
    # 'bpz_library_dir': '../no_elines/',
    # 'z_ini': 1e-3, 'z_fin': 1.5, 'z_delta': 0.001,
    'base_file': '%s/data/BasesDir/Base.bc03.Padova1994.chab.All.hdf5' % home,
    'base_path': 'Base.bc03.Padova1994.chab.All', 'AV_min': 0, 'AV_max': 2,
    'taylor_file': '%s/m_a_interp.hdf5' % save_dir,
    'filters_dir': '%s/doutorado/photo_filters/Alhambra_Filters' % home, 'filters': {'F_365': 'F_365_1.res',
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
                                                                                     'F_954': 'F_954_1.res'},
    'filter_m2l': '%s/doutorado/photo_filters/B_Johnson.res' % home, 'magerr': 0.05,
    'lambda_norm': 4020.,
    'f_l': 1,  # This is a likelihood cooking-factor. The final likelihood is l = l ** f_l
    'n_walkers': 70,
    'n_samples': 1000,
    'n_steps': 3000,
}

config['fit_bpz_outfile'] = 'bpz_fit_magerr_%3.2f_lib_%s_z%.2f_%.2f_nwalk_%i_nsamp_%i_nstep_%i.hdf5' % (
    config['magerr'], os.path.splitext(os.path.basename(config['bpz_library']))[0], config['z_ini'], config['z_fin'],
    config['n_walkers'], config['n_samples'], config['n_steps'])

prefix, suffix = os.path.splitext(config['fit_bpz_outfile'])
config['fit_bpz_post'] = prefix + '_post.hdf5'

cat_version = '4f254712ff'

# 'F_H': 'F_H_1.res',
# 'F_J': 'F_J_1.res',
# 'F_KS': 'F_KS_1.res'}}
