__author__ = 'william'
import numpy as np
import h5py
import matplotlib.pyplot as plt

def get_axes(figure):
    axes = []
    axes.append(figure.add_axes([0.1, 0.5, 0.266666666667, .4]))
    axes.append(figure.add_axes([0.366666666667, 0.5, 0.266666666667, .4], sharey=axes[-1]))
    axes[-1].yaxis.set_ticklabels([])
    axes.append(figure.add_axes([0.633333333333, 0.5, 0.266666666667, .4], sharey=axes[-1]))
    axes[-1].yaxis.set_ticklabels([])
    axes.append(figure.add_axes([0.1, 0.1, 0.8, .35]))
    # axes[-1].yaxis.set_ticklabels([])

    return axes

def make_hist(xmin, xmax, dx, parameter, weights, axis, prior=None):
    bins = np.arange(xmin, xmax, dx)
    aux_hst, aux_bins = np.histogram(parameter, weights=weights, bins=bins)
    aux_hst /= np.max(aux_hst)
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    axis.plot(pts, aux_hst, linewidth=2, color='black')

    # if prior is not None:
    #     aux_hst, aux_bins = np.histogram(prior, bins=bins)
    #     aux_hst /= np.max(aux_hst)
    #     left = np.array(aux_bins[:-1])
    #     right = np.array(aux_bins[1:])
    #     pts = left + (right - left) / 2
    #     axis.plot(pts, aux_hst, linewidth=2, color='blue')

    axis.set_ylim(-.01, 1.01)


file_bpzmags = h5py.File('bpz_templates_mags.hdf5', 'r')
bpz_template_magnitudes = file_bpzmags.get('/bpz_template_magnitudes')
filter_lambdas = file_bpzmags.get('/filter_lambdas')

# file_fit = h5py.File('/Users/william/Downloads/bpz_fit_nointerp_newmask_chi2tx.hdf5', 'r')
file_fit = h5py.File('/Users/william/Downloads/bpz_fit_nointerp_newmask_chi2tx_CCM_test.hdf5', 'r')
f_processed = h5py.File('test.hdf5', 'r')
likelihood = f_processed.get('/likelihood')
parameters = f_processed.get('/parameters')
model_magnitudes = file_fit.get('/model_magnitudes')

f_aux_templates = np.loadtxt('../templates/eB11.list', dtype='S20')

n_z, n_templates, n_fits = likelihood.shape

i_z = 0
for i_template in range(n_templates):
    fig = plt.figure(1)
    fig.clf()
    axes = get_axes(fig)

    # Age historgram
    age_min = 8
    age_max = 10.1
    make_hist(age_min, age_max, 0.1, parameters['age'][i_z, i_template], likelihood[i_z, i_template], axes[0])
    axes[0].set_title('$\\langle \log t \\rangle_\star$')
    # axes[0].xaxis.get_ticklabels()[-1].set_visible(False)
    axes[0].set_xlim(8, 10.1)

    # a_v histogram
    av_min = 0
    av_max = 2
    make_hist(0, 2, 0.2, parameters['a_v'][i_z, i_template], likelihood[i_z, i_template], axes[1])
    axes[1].set_title('$A_V$')
    axes[1].xaxis.get_ticklabels()[0].set_visible(False)
    axes[1].xaxis.get_ticklabels()[-1].set_visible(False)

    # M/L ratio histogram
    make_hist(2, 5, 0.2, parameters['m2l'][i_z, i_template], likelihood[i_z, i_template], axes[2])
    axes[2].set_title('$\log M / L (\lambda = 4020 \\AA)$')
    axes[2].xaxis.get_ticklabels()[-1].set_visible(False)

    alpha = likelihood[i_z, i_template] / likelihood[i_z, i_template].max()
    for i_model in range(len(model_magnitudes[i_z, i_template])):
        axes[-1].plot(filter_lambdas, model_magnitudes[i_z, i_template, i_model] + np.mean(
            bpz_template_magnitudes[i_z, i_template] - model_magnitudes[i_z, i_template, i_model]), color='black',
                 alpha=alpha[i_model])
    axes[-1].plot(filter_lambdas, bpz_template_magnitudes[i_z, i_template], color='red', lw=2, alpha=0.8)
    axes[-1].set_xlabel('$\\lambda [\\AA]$')
    axes[-1].text(8000, 24.5, f_aux_templates[i_template].split('.')[0])
    plt.ylabel('magnitudes')
    plt.ylim(20.5, 25)

    plt.savefig('template_%i.png' % i_template)

    raw_input('Next...')