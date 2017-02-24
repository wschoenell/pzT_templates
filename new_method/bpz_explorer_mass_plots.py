import os

import h5py
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import tables
from bpz_explorer.plots import PlotBPZ

from config import cat_version

alhambra_fit_file = h5py.File('kk_alhambra_fit_fl1.hdf5')
file_dir = '/Users/william/data/alhambra_images'
catalog = alhambra_fit_file['bpz_catalog']

spec_file = os.path.expanduser('~/workspace/pzT_templates/templates/eB11.list')
template_names = [x.split('_')[0] for x in np.loadtxt(spec_file, dtype="S20")]

pdf_file = '/Users/william/data/alhambra_gold_feb2016/alhambragold_added_%s_1e-4_B13v6_eB11.h5' % cat_version
pdf = tables.File(pdf_file)

# tmp:
field = 2
pointing = 1
ccd = 2
# mask = np.bitwise_and(np.bitwise_and(catalog['Field'] == field, catalog['Pointing'] == pointing), catalog['CCD'] == ccd)
mask = np.ones(len(catalog['Field']), dtype=bool)
# mask = np.bitwise_and(catalog['Field'] == field, mask)
mask = np.bitwise_and(mask, catalog["stell"] < .4)
mask = np.bitwise_and(mask, catalog["MS"] > 0)
mask = np.bitwise_and(mask, catalog["MS"] > 0)
mask = np.bitwise_and(mask, catalog['F814W'] < 22.764)
# mask = np.bitwise_and(mask, catalog['F814W'] > 18)
mask = np.bitwise_and(mask, pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["odds"] > .8)
mask = np.bitwise_and(mask, pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Mabs"] < -17)
mask = np.bitwise_and(mask, pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["zml"] > 0.05)
mask = np.bitwise_and(mask, pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["chi2"] < .5)
catalog = catalog[mask]
# tmp end
gal_parameters_bins = alhambra_fit_file["gal_parameters_bins"]
gal_parameters_likelihood = alhambra_fit_file["gal_parameters_likelihood"][mask, :, :]

i_par = int(np.argwhere(alhambra_fit_file['gal_parameters_names'].value == 'mass'))
### Plot mass W vs. mass Tx.
aux_mass = [np.average(np.log10(gal_parameters_bins[..., i_par]), weights=gal_parameters_likelihood[..., i_par][i]) for
            i in range(len(catalog))]
plt.figure(1)
plt.clf()
plt.scatter(catalog["MS"], aux_mass, c=pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Tml"][mask])
plt.plot([7, 12], [7, 12])
plt.xlim(7, 12)
plt.ylim(7, 12)
plt.xlabel("BPZ - Taylor")
plt.ylabel("Willy")
plt.colorbar()
plt.draw()
plt.figure(4)
plt.clf()
plt.hist(catalog["MS"] - aux_mass, bins=50)
plt.title("%3.2f +/- %3.2f" % (np.mean(catalog["MS"] - aux_mass), np.std(catalog["MS"] - aux_mass)))
plt.xlabel("BPZ - Willy")
plt.draw()

### Plot check abs mags
pzt = np.zeros((len(alhambra_fit_file["gal_alhambra_seq_id"]), len(pdf.root.z), len(pdf.root.xt)), "float")
# pzt = np.zeros((n_galaxies, len(h5file.root.z), len(h5file.root.xt)), "float")
# for j, x in enumerate(h5file.root.Posterior[:pzt.shape[0]]):
i_gal = 0
for j in alhambra_fit_file["gal_alhambra_seq_id"]:
    gz = pdf.root.goodz[j]
    gt = pdf.root.goodt[j]
    if pdf.root.Posterior[j].sum() > 0:
        pzt[i_gal][np.outer(gz, gt)] += (pdf.root.Posterior[j] / pdf.root.Posterior[j].sum())
    i_gal += 1

plt.figure(3)
plt.clf()
aux_absmag = pdf.root.Absolute_Magnitude_zT_for_m0eq20[:100] - 20 + \
             pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["m0"][:, np.newaxis, np.newaxis]
pzt = pzt[:, :aux_absmag.shape[1], :]
pzt /= pzt.sum(axis=(1, 2))[:, np.newaxis, np.newaxis]
absmag = np.average(aux_absmag, weights=pzt, axis=(1, 2))[mask]
plt.scatter(pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Mabs"][mask], absmag,
            c=pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Tml"][mask])
plt.plot([-20, -10], [-20, -10])
plt.colorbar()
plt.xlabel('BPZ cat')
plt.ylabel('BPZ pdf')
plt.draw()

old_mass = None
old_absmag = None
# for i_gal in np.argsort(catalog["area"])[::-1]:
for i_gal in np.argsort((catalog["MS"] - aux_mass) ** 2)[::-1]:
# for i_gal in np.argsort((pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Mabs"][mask] - absmag) ** 2)[
#              ::-1]:
    img_file = '%s/f0%dp0%d_OPTICAL_%d.png' % (file_dir, catalog["Field"][i_gal], catalog["Pointing"][i_gal],
                                               catalog["CCD"][i_gal])
    img = mpimg.imread(img_file)[::-1, ...]

    bpz_plot = PlotBPZ(img, catalog, pdf, bpz_template_names=template_names, i_figure=2)
    bpz_plot.plot_dossier(i_gal)
    bpz_plot.figure.axes[3].hist(np.log10(gal_parameters_bins[..., i_par]),
                                 weights=gal_parameters_likelihood[i_gal, ..., i_par], bins=20, normed=True)
    # bpz_plot.figure.axes[3].hist(np.log10(alhambra_fit_file["gal_parameters_bins"][..., i_par]),
    bpz_plot.figure.axes[3].plot([catalog[i_gal]["MS"]], [.2], "*", color="yellow")

    plt.figure(1)
    if old_mass is not None:
        plt.plot(old_mass[0], old_mass[1], '*', color="red")
    plt.plot(catalog[i_gal]["MS"], aux_mass[i_gal], '*', color="yellow")
    old_mass = [catalog[i_gal]["MS"], aux_mass[i_gal]]
    plt.draw()

    plt.figure(3)
    if old_absmag is not None:
        plt.plot(old_absmag[0], old_absmag[1], '*', color="red")
    plt.plot(pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Mabs"][mask][i_gal], absmag[i_gal], '*',
             color="yellow")
    old_absmag = [pdf.root.bpz[alhambra_fit_file["gal_alhambra_seq_id"].value]["Mabs"][mask][i_gal], absmag[i_gal]]
    plt.draw()

    raw_input('delta_M = %3.2f. ENTER for next...' % (aux_mass[i_gal] - catalog[i_gal]["MS"]))
