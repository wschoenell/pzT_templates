import h5py
import numpy as np
from astropy.io import ascii

alhambra_fields = {'f02': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]],  # field, pointing, ccd
                   'f03': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]],
                   'f04': [[1, [1, 2, 3, 4]]],
                   'f05': [[1, [1, 2, 3, 4]]],
                   'f06': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]],
                   'f08': [[1, [1, 2, 3, 4]], [2, [1, 2, 3, 4]]]
                   }


s = ascii.SExtractor()

for field in alhambra_fields.keys():
    for pointing in alhambra_fields[field]:
        for ccd in pointing[1]:
            field_id = '%sP0%iC0%i' % (field.upper(), pointing[0], ccd)
            a = s.read('/Users/william/doutorado/Alhambra/catalogs/latest/alhambra.%s.ColorProBPZ.cat' % field_id)
            b = h5py.File('/Volumes/unsafe RAID 0/new_pdfs/compressed_alhambra.%s.ColorProBPZ.hdf5' % field_id)['FullProbability']
            print field_id, len(a), len(b), len(a) - len(b)