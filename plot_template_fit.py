import h5py
import numpy as np
import matplotlib.pyplot as plt

__author__ = 'william'

f_weights = h5py.File('/Users/william/tmp_out/zT_weights_noelines.hdf5', 'r')
tpl_params = f_weights.get('/tpl_params')
zT_weights = np.ma.masked_invalid(f_weights.get('zT_weights'))

age_min = 8
age_max = 10.1
bins = np.arange(age_min, age_max, 0.05)

for i_t in range(0, 81, 8):
    plt.clf()
    aux_hst, aux_bins = np.histogram(tpl_params['mean_age'], bins=bins, normed=True)
    aux_hst /= np.max(aux_hst)
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    plt.plot(pts, aux_hst, linewidth=2, color='black')  # , label=label)
    aux_weight = zT_weights.sum(axis=1)[:, i_t]
    aux_weight /= aux_weight.sum()
    aux_hst, aux_bins = np.histogram(tpl_params['mean_age'], weights=aux_weight, bins=bins, normed=True)
    aux_hst /= np.max(aux_hst)
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    plt.plot(pts, aux_hst, linewidth=2, color='blue')  # , label=label)
    plt.savefig('pdf_%i.png' % i_t)

f_bpz_mags = h5py.File('/Users/william/tmp_out/pzT_parametric_noelines.hdf5', 'r')
bpz_mags = np.ma.masked_invalid(f_bpz_mags['/templates_bpz'])
tpl_mags = np.ma.masked_invalid(f_weights['tpl_magnitudes'])

i_z = 10
for i_t in range(0, 81, 8):
    plt.clf()
    l_ = range(bpz_mags.shape[2])
    plt.plot(l_, bpz_mags[i_z, i_t])
    model = np.ma.masked_array([tpl_mags[i_z][:,i_m] * zT_weights[:, i_z, i_t] / zT_weights[:, i_z, i_t].sum() for i_m in range(tpl_mags.shape[-1])]).sum(axis=1)
    plt.plot(l_, model)
    plt.savefig('model_%i.png' % i_t)