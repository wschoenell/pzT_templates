import numpy as np

z = {'z_ini': 1e-3, 'z_fin': 1e-1, 'z_delta': 1e-3,}

aux_zfin = z['z_fin']
for i in range(20):
    # print "%6.4f, %6.4f, %6.4f" % (z['z_ini'], z['z_fin'], z['z_delta'])
    arange = np.arange(z['z_ini'], z['z_fin'] + z['z_delta'], z['z_delta'])[:100]
    # print len(arange), arange
    print "    'z_ini': %6.5f, 'z_fin': %6.5f, 'z_delta': %6.5f," % (z['z_ini'], z['z_fin'], z['z_delta'])
    aux = z['z_fin']
    z['z_fin'] += aux_zfin
    z['z_ini'] = aux + z['z_delta']

