import bpz_tools as B
import numpy as np
import matplotlib.pyplot as plt

z = np.arange(1e-3, .2, 1e-3)
m_0 = 14 #np.arange(14, 18, .1)
template_names = np.loadtxt('SED/eB11.list', dtype='S20')

p_i = B.prior(z, m_0, "SM", 11, 0)

plt.clf()
[plt.plot(z, p_i[:,i], label=template_names[i].split('.')[0]) for i in range(p_i.shape[1])]
plt.title('m_0 = %3.2f' % m_0)
plt.xlabel('z')
plt.ylabel('p(z,T|m_0)')
plt.legend()
plt.savefig('prior_bpz.pdf')