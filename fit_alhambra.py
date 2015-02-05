import multiprocessing

__author__ = 'william'

import h5py
import numpy as np

f_wei = h5py.File('/Users/william/tmp_out/zT_weights_kk.hdf5', 'r')
tpl_params = f_wei.get('/tpl_params')
zT_weights = f_wei.get('/zT_weights') #.value[0:1000]
zT_weights = np.ma.masked_array(zT_weights, mask=np.isnan(zT_weights))
# zT_weights.mask += zT_weights >= 0.01

# zT_weights /= np.ma.sum(zT_weights)  # Normalizing...

f_bpz = h5py.File('/Users/william/tmp_out/new_pdfs/compressed_f06p01_colorproext_1_ISO.hdf5', 'r')
bpz_pdf = f_bpz.get('/FullProbability')[:, :zT_weights.shape[1], :]  # TODO: check if z(zT_weights) == z(bpz_pdf)
# bpz_pdf = np.ma.masked_array(bpz_pdf, mask=zT_weights.mask)
# p_Sk = np.dot(zT_weights, np.swapaxes(bpz_pdf,2,1))  # N_S x N_z

# bpz_pdf.shape = (n_obj, n_z, n_t)     (a, b, c)
# zT_weights.shape = (n_S, n_z, n_t)    (d, b, c)
# p_Sk.shape = (n_S, n_z)               (d, b)

(n_obj, n_z, n_t) = bpz_pdf.shape
n_s = zT_weights.shape[0]

# f_out = h5py.File('test_jpmtg_d.hdf5', 'w')
f_out = h5py.File('test_jpmtg_kk.hdf5', 'w')
print 'f_out'
p_s = f_out.create_dataset('/p_s', shape=(n_obj, n_s, n_t)) # , n_t

print '(n_obj, n_s, n_z, n_t): ', n_obj, n_s, n_z, n_t

def aux_p_s(i):
    if i % 10 == 0:
        print 'i_obj>', i
    return i, np.sum(bpz_pdf[i] * zT_weights, axis=1)


def params(i_start, i_end):
    for i in range(i_start, i_end):
        yield i


pool = multiprocessing.Pool()
map_function = pool.map

t_last = 0
# for i_obj in range(n_obj):
#     if i_obj % 10 == 0:
#         print 'i_obj>', i_obj, time.time() - t_last
#         t_last = time.time()

aux_step = 200
for iterac in range(0, n_obj, aux_step):
    if iterac+aux_step-1 > n_obj:
        aux_end = n_obj
    else:
        aux_end = iterac+aux_step
    print 'caca'
    print range(iterac, aux_end)
    aux_i = params(iterac, aux_end)

    for result in map_function(aux_p_s, aux_i):
        i_obj, ps = result
        p_s[i_obj] = ps #/ ps.sum()

# p_s = p_s.sum(axis=-1).sum(axis=-1)

f_bpz.close()
f_wei.close()
f_out.close()