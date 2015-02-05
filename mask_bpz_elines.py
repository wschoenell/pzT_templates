__author__ = 'william'

import numpy as np
import matplotlib.pyplot as plt

tpl = np.loadtxt('templates/eB11.list', np.str)
interpolations = np.linspace(0, len(tpl) - 1, (len(tpl) - 1) * 8 + 1)
i_interpolations = np.array(interpolations, dtype=np.int)

mask_eline = {'Scd_B10.sed': [[3650, 3900], [4891, 5109], [6490, 6810], [7072, 7185]],
              'Sbc_B10.sed': [[4891, 5109]],
              'SB1_B10.sed': [[3650, 3900], [4891, 5109], [5820, 5920], [6245, 6346], [6490, 6810], [7072, 7185]],
              'SB2_B10.sed': [[3650, 3900], [4891, 5109], [5820, 5920], [6245, 6346], [6490, 6810], [7072, 7185]],
              'SB3_B10.sed': [[3633, 4104], [4435, 4500], [4802, 5109], [5820, 5920], [6245, 6346], [6490, 6810], [7072, 7185], [7306, 7366]],
              'SB11_A_0_l.sed': [[3633, 3800], [4070, 4150], [4300, 4400], [4800, 4900], [4925, 5050], [6490, 6810]],
              'SB11_B10.sed': [[3650, 3900], [4891, 5109], [5820, 5920], [6245, 6346], [6490, 6810], [7072, 7185]]
             }
              # 'SB1_B10.sed': [[3650, 3900], [4891, 5109], [5820, 5920], [6490, 6810], [7072, 7185]],
              #'SB3_B10.sed': [[3650, 3900], [4891, 5109], [6490, 6810], [7072, 7185]],

for f_template in tpl:
    print f_template
    s1 = np.loadtxt('templates/%s' % f_template).T

    plt.clf()
    plt.plot(s1[0], s1[1], color='black')

    mask = np.zeros_like(s1[0], dtype=np.bool)

    try:
        for aux_eline in mask_eline[f_template]:
            aux_mask = np.bitwise_and(s1[0] >= aux_eline[0], s1[0] <= aux_eline[1])
            mask += aux_mask
            plt.plot(s1[0][aux_mask], s1[1][aux_mask], color='red')
            plt.plot(s1[0][aux_mask], np.interp(s1[0], s1[0][~aux_mask], s1[1][~aux_mask])[aux_mask], color='green')
    except KeyError:
        print 'Passed %s.' % f_template
        pass

    aux_out = np.interp(s1[0], s1[0][~mask], s1[1][~mask])
    plt.plot(s1[0], aux_out, color='blue')
    np.savetxt('no_elines/%s' % f_template, np.array([s1[0], aux_out]).T)

    plt.xlim(3000,9000)

    raw_input('Next...')



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
#return templates