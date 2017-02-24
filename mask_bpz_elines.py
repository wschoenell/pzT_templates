__author__ = 'william'

import numpy as np
import matplotlib.pyplot as plt

tpl = np.loadtxt('templates/eB11.list', np.str)
interpolations = np.linspace(0, len(tpl) - 1, (len(tpl) - 1) * 8 + 1)
i_interpolations = np.array(interpolations, dtype=np.int)

mask_eline = {
    'Scd_B10.sed': [[1875, 1925], [2290, 2370], [3650, 3900], [4891, 5109], [6490, 6810], [7072, 7185], [8980, 9160],
                    [10240, 10450], [10720, 10975], [12750, 12900], [18630, 18890], [26100, 26430]],
    'Sbc_B10.sed': [[4891, 5109], [6460, 6827], [8980, 9290]],
    'SB1_B10.sed': [[1875, 1925], [2290, 2370], [3650, 3900], [4891, 5109], [5820, 5920], [6245, 6346], [6490, 6810],
                    [7072, 7185], [7270, 7370], [8940, 9240], [10220, 10460], [10740, 11031], [12735, 12910],
                    [18635, 18865], [26130, 26360]],
    'SB2_B10.sed': [[1875, 1925], [2290, 2370], [3650, 3900], [4891, 5109], [5820, 5920], [6245, 6346], [6490, 6810],
                    [7072, 7185], [7285, 7375], [8960, 9275], [10255, 10395], [10765, 10975], [12770, 12880],
                    [18645, 18855], [26145, 26356]],
    'SB3_B10.sed': [[1875, 1925], [2290, 2370], [2780, 2815], [3633, 4104], [4435, 4500], [4802, 5109], [5820, 5920],
                    [6245, 6346], [6490, 6810], [7072, 7185], [7306, 7366], [8984, 9120], [8965, 9610], [10245, 10410],
                    [10605, 11025], [12700, 12935], [16320, 16610], [18595, 18900], [19325, 19590], [20435, 20705],
                    [21110, 21815], [26123, 26390]],
    'SB11_A_0_l.sed': [[1875, 1925], [2290, 2370], [2780, 2815], [3633, 3800], [4070, 4150], [4300, 4400], [4800, 4900],
                       [4925, 5050], [6490, 6810]],
    'SB11_B10.sed': [[2780, 2815], [3650, 3900], [4891, 5109], [5820, 5920], [6245, 6346], [6490, 6810], [7072, 7185]]
}
# 'SB1_B10.sed': [[3650, 3900], [4891, 5109], [5820, 5920], [6490, 6810], [7072, 7185]],
# 'SB3_B10.sed': [[3650, 3900], [4891, 5109], [6490, 6810], [7072, 7185]],

fig = plt.figure(1, figsize=(240/72., 2*240/72.))


plotpars = {'legend.fontsize': 8,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'font.size': 10,
            'axes.titlesize': 12,
            'lines.linewidth': 0.5,
            'font.family': 'Times New Roman',
            #             'figure.subplot.left': 0.08,
            #             'figure.subplot.bottom': 0.08,
            #             'figure.subplot.right': 0.97,
            #             'figure.subplot.top': 0.95,
            #             'figure.subplot.wspace': 0.42,
            #             'figure.subplot.hspace': 0.1,
            'image.cmap': 'GnBu',
            }
plt.rcParams.update(plotpars)

fig.set_tight_layout(True)
fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

plt.clf()

for i_template, f_template in enumerate(tpl):
    print f_template
    s1 = np.loadtxt('templates/%s' % f_template).T

    # plt.clf()
    plt.plot(s1[0], s1[1] + 1 * i_template, color='black')

    mask = np.zeros_like(s1[0], dtype=np.bool)

    try:
        for aux_eline in mask_eline[f_template]:
            aux_mask = np.bitwise_and(s1[0] >= aux_eline[0], s1[0] <= aux_eline[1])
            mask += aux_mask
            plt.plot(s1[0][aux_mask], s1[1][aux_mask] + 1 * i_template, color='red')
            # plt.plot(s1[0][aux_mask], np.interp(s1[0], s1[0][~aux_mask], s1[1][~aux_mask])[aux_mask], color='green')
    except KeyError:
        print 'Passed %s.' % f_template
        pass

    aux_out = np.interp(s1[0], s1[0][~mask], s1[1][~mask])
    # plt.plot(s1[0][~mask], s1[1][~mask], color='blue')
    # plt.xlim(1500, 30000)
    plt.xlim(2000, 9000)
    # plt.savefig('no_elines_plots/template_%i.png' % i_template)
    # np.savetxt('no_elines/%s' % f_template, np.array([s1[0], aux_out]).T)

    plt.draw()
    # raw_input('Next...')

plt.ylim(-.01, 17)
plt.xticks([2000+i*1750 for i in range(5)])
# plt.yticks([])
plt.xlabel('$\lambda [\AA]$')
plt.ylabel('$F_\lambda$ + const.')


plt.savefig('paper/figures/templates_noelines.pdf')

#
#
# templates = []
#
# for i in range(len(interpolations) - 1):
#     i_0 = i_interpolations[i]
#     i_1 = i_interpolations[i] + 1
#     i_frac = interpolations[i] - i_interpolations[i]
#     print i_0, i_1, interpolations[i]
#
#     s1 = np.loadtxt('templates/' + tpl[i_0]).T
#     s2 = np.loadtxt('templates/' + tpl[i_1]).T
#     s2_ = np.interp(s1[0], s2[0], s2[1])
#     s_int = s1.copy()
#     s_int[1] = s1[1] * (1 - i_frac) + s2_ * i_frac
#     this_template = np.empty(len(s_int[0]), dtype=np.dtype([('wl', np.float), ('flux', np.float)]))
#     this_template['wl'] = s_int[0]
#     this_template['flux'] = s_int[1]
#     templates.append(this_template)
# s2 = np.loadtxt('templates/' + tpl[i_1], dtype=np.dtype([('wl', np.float), ('flux', np.float)]))
# templates.append(s2)
# np.savetxt( 'template_int_%3.2f.txt' % (interpolations[i+1]), s2 )
# return templates
