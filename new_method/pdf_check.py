import tables
import numpy as np
from config import cat_version

h5file = tables.File('/Users/william/data/alhambra_gold_feb2016/alhambragold_added_%s_1e-4_B13v6_eB11.h5' % cat_version)

pzt = np.zeros((len(h5file.root.Posterior), len(h5file.root.z), len(h5file.root.xt)), "float")
for j, x in enumerate(h5file.root.Posterior):
    gz = h5file.root.goodz[j]
    gt = h5file.root.goodt[j]
    if x.sum() > 0:
        pzt[j][np.outer(gz, gt)] += (x / x.sum())

s = np.sum(pzt[:, :, 32:], axis=1)

# pzt = pzt #[:, :likelihood.shape[0], :]
