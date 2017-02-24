import h5py
import matplotlib.pyplot as plt
import tables
import numpy as np

from config import cat_version, config

f_color = {0: 'blue', 1: 'red', 2: 'green', 3: 'blue', 4: 'brown', 5: 'purple', 6: 'yellow', 7: 'grey', 8: 'red',
           9: 'green', 10: 'brown'}


def plot_hist(x, weights, axis, color='black', bins=50, label=None, max=None, linestyle='solid'):
    if isinstance(weights, np.ma.core.MaskedArray):
        x = x[~weights.mask]
        weights = weights[~weights.mask]
    aux_hst, aux_bins = np.histogram(x, weights=weights, bins=bins)
    if max:
        aux_hst = np.array(aux_hst, dtype=np.float)
        aux_hst /= aux_hst.max()
        aux_hst *= max
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    print '@> sum', np.sum(aux_hst)
    axis.plot(pts, aux_hst, linewidth=2, color=color, label=label, linestyle=linestyle)

# William fit catalog
out_fname = 'alhambra_fit_fl%i.hdf5' % config['f_l']
f_fit = h5py.File(out_fname)

# Alhambra Catalog:
h5file = tables.File('/Users/william/data/alhambra_gold_feb2016/alhambragold_added_%s_1e-4_B13v6_eB11.h5' % cat_version)
z = np.array(h5file.root.z)
pzt = np.zeros((len(h5file.root.Posterior), len(h5file.root.z), len(h5file.root.xt)), "float")
for j, x in enumerate(h5file.root.Posterior):
    gz = h5file.root.goodz[j]
    gt = h5file.root.goodt[j]
    if x.sum() > 0:
        pzt[j][np.outer(gz, gt)] += (x / x.sum())

# pzt = pzt[:, :likelihood.shape[0], :]

for i_gal in range(1, 10000, 900):
    # i_gal = 1
    mask = z < 1.0
    print '@@> P(z,T) (z > 1.0): %.2f' % np.sum(pzt[i_gal][~mask])
    weight = pzt[i_gal][mask]

    plt.figure(1) # BPZ plot
    plt.clf()
    ax = plt.subplot(111)

    xtemplate = np.array(h5file.root.xt, dtype=int)
    xtemplate_unique = np.unique(xtemplate)

    for i in xtemplate_unique:
        i_template = int(np.argwhere(xtemplate == i)[0])
        print i_template
        plot_hist(z[mask], np.sum(weight[:, xtemplate == i], axis=1), ax, color=f_color[i])

    plot_hist(z[mask], np.sum(weight, axis=1), ax, linestyle='dashed')
    plt.title('photo-z')
    plt.draw()

    plt.figure(2)  # M/L PLOT
    plt.clf()
    ax = plt.subplot(111)
    plot_hist(f_fit['gal_parameters_bins'][:, 2], np.ma.masked_invalid(f_fit['gal_parameters_likelihood'][i_gal, :, 2]), ax, linestyle='dashed')
    plt.title('M/L')
    plt.draw()

    plt.figure(2)  # AGE PLOT
    plt.clf()
    ax = plt.subplot(111)
    plot_hist(f_fit['gal_parameters_bins'][:, 2], np.ma.masked_invalid(f_fit['gal_parameters_likelihood'][i_gal, :, 2]), ax, linestyle='dashed')
    plt.title('M/L')
    plt.draw()
    raw_input('next')



# In [21]: f_fit['gal_parameters_likelihood'].shape
# Out[21]: (10000, 999, 4)
#
# In [22]: f_fit['gal_parameters_bins'].shape
# Out[22]: (999, 4)
#
# In [25]: print f_fit['gal_parameters_names'].value
# ['a_v' 'age' 'm2l' 'metallicity']




