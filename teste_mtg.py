

from magal.library import LibraryModel
from magal.io.filterset import FilterSet
from magal.photometry.syntphot import spec2filterset

FSet = FilterSet('JPAS_56.hdf5')
FSet.load('JPAS_56', 1)

lib_param = random_param()  ## 't0_young', 'tau_young', 't0_old', 'tau_old', 'frac_young', 'Z', 'tau_v'
lib = LibraryModel('two_exp', lib_param, bases_dir, base_file, 'light', 4020.)
mags = spec2filterset(f.filterset, lib.get_model_spectrum(0))


