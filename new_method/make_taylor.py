from common import eval_matrix

__author__ = 'william'

## 0 - Configuration
config = {'bpz_library': '../templates/eB11.list', 'n_interpolations': 7, 'bpz_library_dir': '../no_elines/',
          'z_ini': 1e-4, 'z_fin': 7.0, 'z_delta': 0.001,
          'base_file': '/Users/william/BasesDir/Base.bc03.Padova1994.chab.All.hdf5', 'base_path': 'Base.bc03.Padova1994.chab.All',
          'AV_min': 0, 'AV_max': 2,
          'taylor_file': '/Users/william/tmp_out/m_a_interp.hdf5',
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

eval_matrix(config)