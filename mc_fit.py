__author__ = 'william'

import numpy as np

from astropy import units
from astropy.cosmology import WMAP9 as cosmo

from magal.io.readfilterset import FilterSet
from magal.photometry.syntphot import spec2filterset
from magal.util.cosmo import zcor

from parametric_library import bpz_templates

def mag_in_z(spec, z, f):
    # d_L = cosmo.luminosity_distance(z).to('cm')
    # k_cosmo = units.Lsun.to('erg/s') / (4 * np.pi * np.power(d_L, 2))
    model_spec = zcor(spec, z)
    # model_spec['flux'] *= k_cosmo
    return spec2filterset(f.filterset, model_spec)

class FitBpz(object):

    def __init__(self, bpz_mag):
        self.bpz_mag = bpz_mag

    def lnprob(self, p):

        return


f = FilterSet('/Users/william/doutorado/photo_filters/Alhambra_24.hdf5')
f.load('Alhambra_24', 1)

bpz_spectra = bpz_templates()
bpz_mags = [mag_in_z(bpz_spectra[i_template], 0.001, f)['m_ab'] for i_template in range(len(bpz_spectra))]

#save bpz_spectra to a file
for i_template in range(len(bpz_spectra)):
    np.savetxt('templates_interp/bpz_template_%02i.txt' % i_template, bpz_spectra[i_template])