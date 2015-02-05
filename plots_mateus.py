import os
import atpy
import h5py
import numpy as np
import matplotlib.pyplot as plt

f_catalog = '/Users/william/doutorado/Alhambra/catalogs/catalog_v4f/alhambra.F06P01C01.ColorProBPZ.cat'
cat = atpy.Table(f_catalog, type='ascii')

f_weights = h5py.File('/Users/william/tmp_out/zT_weights_kk.hdf5', 'r')
tpl_params = f_weights.get('/tpl_params')
f_ps = h5py.File('test_jpmtg_kk.hdf5', 'r')  # From fit_alhambra
p_s = f_ps.get('/p_s').value  # ) #, axis=2)   #[:7000]  # FIXME:

n_zmax = f_weights.get('zT_weights').shape[1]
f_bpz = h5py.File('/Users/william/tmp_out/new_pdfs/compressed_f06p01_colorproext_1_ISO.hdf5', 'r')
mask_zmax = f_bpz.get('/FullProbability')[:, 1501:, :].sum(-1).sum(-1) < 0.005  # Mask out objects where p(zT| z>1.5)

# for i in range(len(p_s)):
#     p_s[i] = p_s[i] / p_s[i].sum()

# hist_a = np.array([tpl_params['mean_age'] * p_s[i] for i in range(len(p_s))])


# bins=50
# plt.hist(tpl_params['mean_age'], bins=bins, alpha=0.3, normed=True, label='prior')


# mask = (cat['tb_1'] <= 5)[:p_s.shape[0]]

# hist_wei = p_s[mask].sum(axis=0)
# hist_wei /= hist_wei.sum()

# plt.hist(tpl_params['mean_age'], bins=bins, weights=hist_wei, normed=True, label='posterior tb_1 <= 5', alpha=.3)

# mask = (cat['tb_1'] > 5)[:p_s.shape[0]]
# hist_wei = p_s[mask].sum(axis=0)
# plt.hist(tpl_params['mean_age'], bins=bins, weights=hist_wei, normed=True, label='posterior tb_1 > 5', alpha=.3)
#plt.hist(tpl_params['mean_age'], bins=bins, normed=True, alpha=0.3, weights=hist_a[mask].sum(axis=0), label='posterior')
# plt.hist(tpl_params['mean_age'], weights=p_s, bins=bins, alpha=0.3, normed=True, label='posterior')

def plot(p, bins, log=False):
    plt.clf()


    #### PRIOR ####
    if log:
        p_ = np.log10(tpl_params[p])
    else:
        p_ = tpl_params[p]
    aux_hst, aux_bins = np.histogram(p_, bins=bins, normed=True)
    aux_hst /= np.max(aux_hst)

    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    plt.plot(pts, aux_hst, linewidth=2, color='black', label='prior')

    #### POSTERIOR T <= 5 ####
    label = 'posterior t <= 4.5'
    mask = np.bitwise_and(cat['tb_1'][:p_s.shape[0]] <= 4.5, cat.data['MagPrior'][:p_s.shape[0]] < 21.5)
    mask = np.bitwise_and(mask, mask_zmax)
    hist_wei = p_s[mask].sum(axis=(0,2))
    print hist_wei, mask
    aux_hst, aux_bins = np.histogram(p_, weights=hist_wei, bins=bins)# , normed=True)
    aux_hst /= np.max(aux_hst)

    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    plt.plot(pts, aux_hst, linewidth=2, color='red', label=label)

    print '%s: %i objects' % (label, mask.sum())

    #### POSTERIOR T > 5 ####
    label = 'posterior t > 6.5'
    mask = np.bitwise_and(cat['tb_1'][:p_s.shape[0]] > 6.5, cat.data['MagPrior'][:p_s.shape[0]] < 21.5)
    mask = np.bitwise_and(mask, mask_zmax)
    hist_wei = p_s[mask].sum(axis=(0,2))
    aux_hst, aux_bins = np.histogram(p_, weights=hist_wei, bins=bins) #, normed=True)
    aux_hst /= np.max(aux_hst)

    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    plt.plot(pts, aux_hst, linewidth=2, color='blue', label=label)

    print '%s: %i objects' % (label, mask.sum())


    # plt.hist(tpl_params['mean_age'][mask], weights=p_s[mask].sum(axis=1), bins=bins, alpha=0.3, normed=True, label='posterior t >= 5 ')
    xlim = (bins.min(), bins.max())
    plt.xlim(xlim)



    # plt.legend()
    fname = '%s.png' % p
    if os.path.exists(fname):
        print('deleting %s' % fname)
        os.unlink(fname)
    plt.savefig(fname)


age_min = 8
age_max = 10.1

for p_ in ['t0_young', 't0_old']:
    plt.clf()
    bins = np.arange(age_min, age_max, 0.1)
    plot(p_, bins, log=True)

bins = np.arange(age_min, age_max, 0.1)
plot('mean_age', bins)

bins = np.arange(0, 1, 0.05)
plot('frac_young', bins)

tau_v_low = 0
tau_v_upp = 2.0

plt.clf()
bins = np.arange(tau_v_low, tau_v_upp, 0.05)

plot('tau_v', bins)

# plt.clf()
# mask = np.bitwise_and(cat['tb_1'] >= 4.5, cat.data['MagPrior'] < 21.5)
# mask = np.bitwise_and(mask, mask_zmax)
# p_s = np.sum(f_ps.get('/p_s').value[mask], axis=0)
# plt.plot(f_weights['redshift'], [((tpl_params['mean_age'] * p_s[:,i])/p_s[:,i].sum()).sum() for i in range(p_s.shape[1])])
#
# mask = np.bitwise_and(cat['tb_1'] > 6.5, cat.data['MagPrior'] < 21.5)
# mask = np.bitwise_and(mask, mask_zmax)
# p_s = np.sum(f_ps.get('/p_s').value[mask], axis=0)
# plt.plot(f_weights['redshift'], [((tpl_params['mean_age'] * p_s[:,i])/p_s[:,i].sum()).sum() for i in range(p_s.shape[1])])
#
#
# plt.savefig('age_redshift.png')


#tpl_params['mean_age']

#f_ps.close()
#f_weights.close()