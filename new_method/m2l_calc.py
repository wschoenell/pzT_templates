import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pystarlight.util.constants import L_sun
from pystarlight.util.constants import d_sun

from config import home


def calc_m2l_filter(lam, lum, transm):
    """

    :param lam:  Lambda of filter/SSP
    :param lum:  CSP/SSP Luminosity (M_\odot/L_\odot)
    :param transm:  Filter transmission curve
    :return:
    """


m2l_bc03 = np.loadtxt("%s/mestrado/bc03/models/Padova1994/chabrier/bc2003_hr_m62_chab_ssp.4color" % home,
                      usecols=(5, 6), dtype=np.dtype([('M2L_V', np.float), ('fMcor', np.float)]))

plt.figure(1)
plt.clf()
plt.scatter(range(len(m2l_bc03)), m2l_bc03['M2L_V'], s=5)

data_dir = '%s/workspace/pystarlight/data/solar_spectrum/' % home
sun_spec = fits.getdata(data_dir + 'sun_reference_stis_002.fits')
filtercurve = np.loadtxt("%s/Downloads/M2L_filter/V.dat" % home,
                         dtype=np.dtype([('wl', np.float), ('transmission', np.float)]))

# These will be filled on the for loop:
sun_flux = None
filter_transm = None
aux_m2l = list()

for i_ssp in np.arange(1, len(m2l_bc03) + 1) + 1:
    bc03_filename = "%s/BasesDir/bc2003_hr_m62_chab_ssp_%03i.spec" % (home, i_ssp)
    bc03_spec = np.loadtxt(bc03_filename, dtype=np.dtype([('wl', np.float), ('l2m', np.float)]))

    # Resample filter and sun to the BC03 wavelengths
    if filter_transm is None:
        filter_transm = np.interp(bc03_spec["wl"], filtercurve["wl"], filtercurve["transmission"])
    if sun_flux is None:
        sun_flux = np.interp(bc03_spec["wl"], sun_spec["WAVELENGTH"], sun_spec["FLUX"])

    # Calculate the integral of l2m in the filter (this is bolometric!)
    # l2m_filter = np.trapz(bc03_spec["l2m"] * bc03_spec["wl"] * filter_transm, bc03_spec["wl"]) / \
    #              np.trapz((bc03_spec["wl"] ** -1) * filter_transm, bc03_spec["wl"])
    l2m_filter = np.trapz(bc03_spec["l2m"] * filter_transm,
                          bc03_spec["wl"])  # / np.trapz(filter_transm, bc03_spec["wl"])
    # Apply the bolometric -> filter correction
    # solar_luminosity_filter = np.trapz(sun_flux * bc03_spec["wl"] * filter_transm, bc03_spec["wl"]) /\
    #                           np.trapz(filter_transm * (bc03_spec["wl"] ** -1), bc03_spec["wl"])  # in erg/s/cm2
    solar_luminosity_filter = np.trapz(sun_flux * filter_transm, bc03_spec["wl"]) #/\
                             # np.trapz(filter_transm, bc03_spec["wl"])  # in erg/s/cm2
    # from erg/s/cm2 to L_sun:
    solar_luminosity_filter *= (4 * np.pi * d_sun ** 2) / L_sun

    # solar_luminosity_filter = 0.113

    l2m_filter /= solar_luminosity_filter

    aux_m2l.append(m2l_bc03['fMcor'][i_ssp - 2] / l2m_filter)

    plt.scatter(i_ssp + 1, aux_m2l[-1], s=5, c="black", alpha=.2)
    print i_ssp + 1, i_ssp - 1

plt.figure(2)
plt.clf()
plt.scatter(range(len(m2l_bc03["M2L_V"])), (m2l_bc03["M2L_V"] - np.array(aux_m2l)) /  m2l_bc03["M2L_V"])
