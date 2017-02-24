import numpy as np
import h5py

from config import config
from common import av_taylor_coeff

tf = h5py.File(config['taylor_file'], 'r')
m = np.copy(tf['m'].value)
a = np.copy(tf['a'].value)

av = 1
i_z = 0
i_met
av_taylor_coeff(n, av, a[i_z, i_met])