import os

home = os.path.expanduser('~/')
# save_dir = '%s/data/out/' % home
save_dir = './'
config = {
    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library/list_random.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library/list_random.txt' % home,
    # 'bpz_library_true_values_cols': (1, 7, 6),

    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_single/list_random_single.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_single/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_single/list_random_single.txt' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),


    # OLLLLDDD!!!! MOCK *DOUBLE* Burst
    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_double/list_random_double_new.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_double/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_double/list_random_double_new.txt' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),

    # OLLLLDDD!!!! MOCK Single Burst
    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_single_old/list_random_single_old.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_single_old/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_single_old/list_random_single_old.txt' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),

    # MOCK Single Burst
    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_single_burst/list_random_single.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_single_burst/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_single_burst/list_random_single.txt' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),

    # MOCK *DOUBLE* Burst
    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_double_burst/list_random_double.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_double_burst/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_double_burst/list_random_double.txt' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),

    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_plot_taylor/list_random_single.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_plot_taylor/' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_plot_taylor/list_random_single.txt' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),  # age, metallicity, extinction

    # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_plot_taylor_double/list_random_double_fit.txt' % home,
    # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_plot_taylor_double/list_random_double_fit.txt' % home,
    # # 'bpz_library': '%s/workspace/pzT_templates/new_method/random_library_plot_taylor_double/list_random_double.txt' % home,
    # # 'bpz_library_true_values': '%s/workspace/pzT_templates/new_method/random_library_plot_taylor_double/list_random_double.txt' % home,
    # 'n_interpolations': 0,
    # 'bpz_library_dir': '/%s/workspace/pzT_templates/new_method/random_library_plot_taylor_double/' % home,
    # 'bpz_library_true_values_cols': (3, 7, 6),  # age, metallicity, extinction

    # 'bpz_library': '%s/workspace/pzT_templates/templates/Base.BC03.N.list' % home, 'n_interpolations': 0,
    # 'bpz_library_dir': '%s/data/BasesDir/' % home,
    # 'bpz_library_true_values': '%s/mestrado/BasesDir/Base.BC03.N' % home,
    # 'bpz_library_true_values_cols': (1, 2, 4),

    'bpz_library': '%s/workspace/pzT_templates/templates/eB11.list' % home,
    'n_interpolations': 3,
    'bpz_library_dir': '../no_elines/',

    # 'combined': False,
    # 'z_ini': 0.00100, 'z_fin': 0.10000, 'z_delta': 0.00100,
    # 'z_ini': 0.30100, 'z_fin': 0.40000, 'z_delta': 0.00100,
    # 'z_ini': 0.50100, 'z_fin': 0.60000, 'z_delta': 0.00100,
    # 'z_ini': 0.60100, 'z_fin': 0.70000, 'z_delta': 0.00100,
    # 'z_ini': 0.40100, 'z_fin': 0.50000, 'z_delta': 0.00100,
    # 'z_ini': 0.80100, 'z_fin': 0.90000, 'z_delta': 0.00100,
    # 'z_ini': 1.10100, 'z_fin': 1.20000, 'z_delta': 0.00100,

    'combined': True,
    'z_ini': 0.00100, 'z_fin': 0.40000, 'z_delta': 0.00100,
    # 'z_ini': 0.00100, 'z_fin': 1.10000, 'z_delta': 0.00100,

    'taylor_file': '%s/data/out/m_a_interp_1e-3.hdf5' % home,

    # 'z_ini': 1e-3, 'z_fin': 1e-3, 'z_delta': 0.05,
    # 'taylor_file': '%s/data/out/m_a_interp_new.hdf5' % home,

    'base_file': '%s/data/BasesDir/Base.bc03.Padova1994.chab.All.no22and32.hdf5' % home,
    'base_path': 'Base.bc03.Padova1994.chab.All.no22and32',

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
    # 'filter_m2l': '%s/doutorado/photo_filters/B_Johnson.res' % home,
    # 'filter_m2l': '%s/doutorado/photo_filters/sdss/i.dat' % home,
    # 'filter_m2l': '/Users/william/Downloads/M2L_filter/V.dat',
    'filter_m2l': '%s/doutorado/photo_filters/Alhambra_Filters/HST_F814W.res' % home,
    'magerr': 0.05,
    'mag_noise': 0,  # .05,
    'm2l_lambdas': [4020, 5000, 7000],
    'lambda_norm': 5500.,
    'f_l': 1,  # This is a likelihood cooking-factor. The final likelihood is l = l ** f_l
    ### Standard
    'n_walkers': 70,
    'n_samples': 1000,
    'n_steps': 3000,
    'save_chains': 7,
    ### Extended
    # 'n_walkers': 70,
    # 'n_samples': 2000,
    # 'n_steps': 6000,
    # 'save_chains': 14,
    ### Test
    # 'n_walkers': 70,
    # 'n_samples': 700,
    # 'n_steps': 1000,
    # 'save_chains': 70,
    'chi2_spectrum': False,  # If True, evaluates the chi2 with the spectrum, not using the AV approximation
    # Prior configuration
    # Extinction
    'AV_min': -0.1, 'AV_max': 2,
    # Old component:
    # 't0_old_min': 5e9 + 1,  - Now defined if single/double burst
    # t0_old_max is the age of the Universe at given redshift
    'tau_t0_old_min': 0.1,
    'tau_t0_old_max': 10,
    # Young component:
    't0_young_min': 1e6,
    't0_young_max': 5e9,
    'tau_t0_young_min': 0.1,
    'tau_t0_young_max': 10,
    # Maximum light fraction on the Young component
    ## Single burst - Only "old" component
    't0_old_min': 1e6,
    'lg_frac_young_min': -1e30,
    'lg_frac_young_max': -1e30,
    ## Double burst
    # 't0_old_min': 5e9 + 1,
    # 'lg_frac_young_min': -np.inf,
    # 'lg_frac_young_max': 0
}

config['fit_bpz_outfile'] = 'bpz_fit_magerr_%s_%6.4f_lib_%s_z%.2f_%.2f_nwalk_%i_nsamp_%i_nstep_%i.hdf5' % (
    'single' if config['lg_frac_young_min'] == -1e30 else 'double',
    # config['fit_bpz_outfile'] = 'specmag_bpz_fit_magerr_%3.2f_lib_%s_z%.2f_%.2f_nwalk_%i_nsamp_%i_nstep_%i.hdf5' % (
    config['magerr'], os.path.splitext(os.path.basename(config['bpz_library']))[0], config['z_ini'], config['z_fin'],
    config['n_walkers'], config['n_samples'], config['n_steps'])

if config['combined']:
    config['fit_bpz_outfile'] = "combined_" + config['fit_bpz_outfile']

# Append magnitude noise on the end of output filename, if necessary.
if config['mag_noise'] != 0:
    f, e = os.path.splitext(config['fit_bpz_outfile'])
    config['fit_bpz_outfile'] = f + '_noise_%.3f' % config['mag_noise'] + e

if config["chi2_spectrum"]:
    config['fit_bpz_outfile'] = "specmag_" + config['fit_bpz_outfile']

# config['fit_bpz_outfile'] = 'test_'+config['fit_bpz_outfile']

prefix, suffix = os.path.splitext(config['fit_bpz_outfile'])
config['fit_bpz_post'] = prefix + '_post.hdf5'

cat_version = '4f254712ff'


# 'F_H': 'F_H_1.res',
# 'F_J': 'F_J_1.res',
# 'F_KS': 'F_KS_1.res'}}