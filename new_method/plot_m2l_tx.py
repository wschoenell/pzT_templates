import h5py
import matplotlib.pyplot as plt
import numpy as np

from common import load_filters, mag_in_z
from config import config, home

plt.ion()

config.update({'filters_dir': '%s/doutorado/photo_filters/sdss' % home, 'filters': {'g': 'g.dat',
                                                                                    'i': 'i.dat'}, })

filters, _ = load_filters(config)

templates = np.loadtxt(config['bpz_library'], np.str, usecols=(0,))
templates_data = [
    np.loadtxt('%s/%s' % (config['bpz_library_dir'], t), dtype=np.dtype([('lambda', np.float), ('flux', np.float)])) for
    t in templates]

interpolations = np.linspace(0., len(templates) - 1, len(templates) + (len(templates) - 1) * config['n_interpolations'])

templates_data_interpolated = list()
for interpolation in interpolations:
    if interpolation != round(interpolation):
        i_0 = int(interpolation)
        i_1 = int(interpolation + 1)
        print interpolation, i_0, i_1, interpolation - i_0, i_1 - interpolation
        aux = np.copy(templates_data[i_0])
        aux['flux'] *= interpolation - i_0
        aux['flux'] += np.interp(aux['lambda'], templates_data[i_1]['lambda'], templates_data[i_1]['flux']) * (
            i_1 - interpolation)
        templates_data_interpolated.append(aux)
    else:
        print interpolation
        templates_data_interpolated.append(templates_data[int(interpolation)])

m2l_taylor = list()
gi = list()
for i_template in range(len(templates_data_interpolated)):
    x = mag_in_z(templates_data_interpolated[i_template]['lambda'], templates_data_interpolated[i_template]['flux'], 0,
                 filters)
    gi.append(x[0] - x[1])
    # m2l_taylor.append(1.15 + 0.7 * gi[-1])
    m2l_taylor.append(-0.68 + 0.7 * gi[-1])

plt.figure(1)
plt.clf()
plt.scatter(np.arange(len(templates_data_interpolated)), gi)
plt.xlabel("template")
plt.ylabel("g - i")



fit_bpz_post = h5py.File("bpz_fit_magerr_single_0.0500_lib_eB11_z0.00_0.10_nwalk_70_nsamp_1000_nstep_3000_post.hdf5")
m2l = fit_bpz_post['parameters']["m2l"]
plt.figure(2)
plt.clf()
for i_template in range(m2l.shape[1]):
    plt.scatter(m2l.shape[2] * [i_template], np.log10(m2l[0, i_template]))

plt.figure(5)
plt.clf()
i_z = 74
for i_template in range(m2l.shape[1]):
    # FIXME: tirar esse for!
    plt.scatter(m2l.shape[2] * [gi[i_template]], np.log10(m2l[0, i_template]), alpha=.2, c=np.log10(fit_bpz_post['likelihood'][i_z, i_template])/np.log10(fit_bpz_post['likelihood'][i_z, i_template]).max())
# plt.figure(2)
# for i_z in range(10):
#     plt.scatter(range(m2l.shape[1]), np.average(m2l[i_z], weights=fit_bpz_post['likelihood'][i_z], axis=-1))
plt.colorbar()

plt.figure(2)
plt.scatter(np.arange(len(templates_data_interpolated)), m2l_taylor, color='black')
plt.scatter(range(m2l.shape[1]), np.average(np.log10(m2l[i_z]), weights=fit_bpz_post['likelihood'][i_z], axis=-1), color='magenta')
plt.xlabel("template")
plt.ylabel("log M/L")

plt.figure(5)
plt.scatter(gi, m2l_taylor, color='black')
i_z = 0
plt.scatter(gi, np.average(np.log10(m2l[i_z]), weights=fit_bpz_post['likelihood'][i_z], axis=-1), color='magenta')
plt.xlabel("(g - i)")
plt.ylabel("log M/L")

plt.figure(3)
plt.clf()
plt.scatter(gi, m2l_taylor)
plt.plot([-.2, 1.6], np.array([-.2, 1.6])*.7 - .68)
plt.xlabel("g - i")
plt.ylabel("M/L i-band")