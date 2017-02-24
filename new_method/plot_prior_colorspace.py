import os

import asciitable
import atpy
import numpy as np
import matplotlib.pyplot as plt
from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.base import StarlightBase
from pystarlight.util.redenninglaws import Cardelli_RedLaw
from config import config
from common import load_filters, mag_in_z, calc_LsunFilter
from model import Model

filters_plot = [['F_458', 'F_768'], ['F_365', 'F_644']]

i_z = 0

fnumbers = np.sort(config['filters'].keys())
fid = [[int(np.argwhere(fnumbers == j)) for j in i] for i in filters_plot]

# For the A_V Arrow:
e_c1 = Cardelli_RedLaw([4580])[0] - Cardelli_RedLaw([7680])[0]
e_c2 = Cardelli_RedLaw([3650])[0] - Cardelli_RedLaw([6440])[0]

# M/L calculations
## Filter
filter_m2l = readfilterfile(config['filter_m2l'], norm=False).data
aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
filter_m2l_new['lambda'] = aux_l
filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
filter_m2l = filter_m2l_new
LsunFilter = calc_LsunFilter(filter_m2l)

## Base
bt = StarlightBase(config['base_file'], config['base_path'])
m2l_ssp = np.empty((len(bt.metBase), len(bt.ageBase)))
i_norm = np.argmin((bt.l_ssp - config['lambda_norm']) ** 2)
l2m_norm_ssp = np.empty((len(bt.metBase), len(bt.ageBase)))
for i_met in range(len(bt.metBase)):
    for i_age in range(len(bt.ageBase)):
        aux = np.interp(filter_m2l['lambda'], bt.l_ssp, bt.f_ssp[i_met, i_age])
        m2l_ssp[i_met, i_age] = bt.Mstars[i_met, i_age] * LsunFilter / np.trapz(aux * filter_m2l['transm'],
                                                                                filter_m2l['lambda'])
        l2m_norm_ssp[i_met, i_age] = bt.f_ssp[i_met, i_age, i_norm]

# Open BPZ magnitudes
templates = np.loadtxt(config['bpz_library'], np.str)
aux = templates[:, 0] if len(templates.shape) > 1 else templates
templates_data = [
    np.loadtxt('%s/%s' % (config['bpz_library_dir'], t), dtype=np.dtype([('lambda', np.float), ('flux', np.float)]))
    for
    t in aux]

filters, _ = load_filters(config)

templates_mags = list()
m2l = list()
at_mass = list()
at_flux = list()
m2l_filter = list()
m2l_csp_spec = list()

# FIXME:
bt.ageBase[0] = 1e5

q = Cardelli_RedLaw(bt.l_ssp)

for i_template in range(len(templates_data)):
    templates_mags.append(
        mag_in_z(templates_data[i_template]['lambda'], templates_data[i_template]['flux'], config['z_ini'], filters))
    # m2l.append(1 / templates_data[i_template]['flux'][templates_data[i_template]['lambda'] == 5500])

    t0_young = float(templates[i_template][1])
    tau_young = float(templates[i_template][2])
    t0_old = float(templates[i_template][3])
    tau_old = float(templates[i_template][4])
    frac_young = float(templates[i_template][5])
    sfh = Model.get_sfh(bt, t0_young, tau_young, t0_old, tau_old, frac_young)
    i_met = int(np.argwhere(bt.metBase == float(templates[i_template][7])))

    at_flux.append(np.sum(sfh * l2m_norm_ssp[i_met] * np.log10(bt.ageBase) / np.sum(sfh * l2m_norm_ssp[i_met])))
    at_mass.append(np.sum(sfh * np.log10(bt.ageBase)))
    m2l_filter.append(1 / np.sum(sfh / m2l_ssp[i_met]))

    spec = (bt.f_ssp[i_met] * sfh[:, np.newaxis]).sum(axis=0)  # bt.Mstars[i_met][:, np.newaxis]
    a_v = float(templates[i_template][6])
    # spec *= 10 ** (-0.4 * (q * a_v))
    aux = np.interp(filter_m2l['lambda'], bt.l_ssp, spec)
    m2l_csp_spec.append(LsunFilter / np.trapz(aux * filter_m2l['transm'], filter_m2l['lambda']))
    m2l.append(1 / spec[bt.l_ssp == 5500])

    # plt.figure(10)
    # plt.clf()
    # plt.plot(bt.l_ssp, spec)
    # plt.plot(templates_data[i_template]['lambda'], templates_data[i_template]['flux'])
    # plt.xlim(3000, 9000)
    # plt.draw()
    # raw_input('next...')

at_mass = np.array(at_mass)
templates_mags = np.array(templates_mags)
m2l = np.array(m2l)
m2l_csp_spec = np.array(m2l_csp_spec)

color1 = np.ravel(templates_mags[:, fid[0][0]] - templates_mags[:, fid[0][1]])
color2 = np.ravel(templates_mags[:, fid[1][0]] - templates_mags[:, fid[1][1]])

# Template file variables:
#  t0_young, tauT_young, t0_old, tauT_old, frac_young, a_v, lg_metallicity

plt.clf()

ax = plt.subplot(2, 3, 1)
ax.scatter(color1, color2, c=np.array(at_flux, dtype=float), s=5, lw=0)
ax.set_title('at_flux @ 4020 AA')

ax = plt.subplot(2, 3, 2)
ax.scatter(color1, color2, c=np.array(templates[:, 6], dtype=float), s=5, lw=0)
ax.set_title('A_V')

ax = plt.subplot(2, 3, 4)
ax.scatter(color1, color2, c=np.array(templates[:, 7], dtype=float), s=5, lw=0)
ax.set_title('Z')
ax.set_xlabel("%s - %s" % (filters_plot[0][0], filters_plot[0][1]))
ax.set_ylabel("%s - %s" % (filters_plot[1][0], filters_plot[1][1]))

ax = plt.subplot(2, 3, 5)
ax.scatter(color1, color2, c=np.log10(m2l), s=5, lw=0)
ax.set_title('log M/L @ 4020 AA')

ax = plt.subplot(2, 3, 6)
ax.scatter(color1, color2, c=np.log10(m2l_filter), s=5, lw=0)
ax.set_title('log M/L_{%s}' % os.path.splitext(os.path.basename(config['filter_m2l']))[0])

ax = plt.subplot(2, 3, 3)

# Opens bc03 file
bc03color = np.loadtxt('bc2003_hr_m62_chab_ssp.4color', usecols=(0, 4, 5),
                       dtype=np.dtype([('log_age', float), ('BC03_M2L_B', float), ('BC03_M2L_V', float)]))
bc03color = bc03color[np.bitwise_and(bc03color['log_age'] > 9.0, bc03color['log_age'] < 10.)]

ax.scatter(at_mass, m2l_filter, c=np.log10(m2l), s=5, lw=0)
ax.scatter(bc03color['log_age'], bc03color['BC03_M2L_V'], color='red', s=5, lw=0)

ax.set_xlabel('at_mass')
# ax.set_ylabel('log M/L_{4020 AA}')
ax.set_ylabel('M/L_{%s}' % os.path.splitext(os.path.basename(config['filter_m2l']))[0])

### Plot BC03 over M/L ###
bc03color = np.loadtxt('bc2003_hr_m62_chab_ssp.4color', usecols=(0, 4, 5),
                       dtype=np.dtype([('log_age', float), ('BC03_M2L_B', float), ('BC03_M2L_V', float)]))
dt = np.dtype([('fname', 'S40'), ('age', float), ('metallicity', float)])
bc03_base = np.loadtxt('/Users/william/mestrado/BasesDir/Base.bc03.Padova1994.chab.m62', dtype=dt, usecols=(0, 1, 2),
                       skiprows=1)

aux = list()
for i_file in range(len(bc03color)):
    arq_base = np.loadtxt('/Users/william/mestrado/BasesDir/%s' % bc03_base['fname'][i_file],
                          dtype=np.dtype([('l', float), ('f', float)]))
    aux.append(1 / arq_base['f'][arq_base['l'] == 5500])

plt.figure(2)
plt.clf()
plt.scatter(aux, bc03color['BC03_M2L_V'])

i_met = int(np.argwhere(bt.metBase == 0.02))
plt.scatter(aux, m2l_ssp[i_met][1:], color='red', alpha=.2)
