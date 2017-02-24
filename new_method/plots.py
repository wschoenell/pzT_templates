import os
from magal.plots.mosaic import get_mosaic
from pystarlight.util.base import StarlightBase

from config import config, save_dir, prefix
from common import load_filters
from model import Model

__author__ = 'william'
import numpy as np
import h5py
import matplotlib.pyplot as plt

i_z = 0

paper = False


# Change figure size to A4
# matplotlib.rcParams.update({'figure.figsize': (16, 25)})

# config = {'magerr': .05}

# TODO: Binar Z no Z_base e colocar estrelas em cada um dos elementos da base em Z

def get_axes(figure):
    axes = []
    axes.append(figure.add_axes([0.1, 0.5, 0.2, .4]))
    axes.append(figure.add_axes([0.3, 0.5, 0.2, .4], sharey=axes[-1]))
    axes[-1].yaxis.set_ticklabels([])
    axes.append(figure.add_axes([0.5, 0.5, 0.2, .4], sharey=axes[-1]))
    axes[-1].yaxis.set_ticklabels([])
    axes.append(figure.add_axes([0.7, 0.5, 0.2, .4], sharey=axes[-1]))
    axes[-1].yaxis.set_ticklabels([])
    axes.append(figure.add_axes([0.1, 0.1, 0.8, .35]))
    # axes[-1].yaxis.set_ticklabels([])

    return axes


def make_hist(xmin, xmax, dx, parameter, weights, axis, color='black', prior=None, bins=None, label=None, alpha=1):
    if bins is None:
        bins = np.arange(xmin, xmax, dx)
    aux_hst, aux_bins = np.histogram(parameter, weights=weights, bins=bins)
    aux_hst /= np.max(aux_hst)
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts = left + (right - left) / 2
    axis.plot(pts, aux_hst, linewidth=2, color=color, label=label, alpha=alpha)

    # if prior is not None:
    #     aux_hst, aux_bins = np.histogram(prior, bins=bins)
    #     aux_hst /= np.max(aux_hst)
    #     left = np.array(aux_bins[:-1])
    #     right = np.array(aux_bins[1:])
    #     pts = left + (right - left) / 2
    #     axis.plot(pts, aux_hst, linewidth=2, color='blue')

    axis.set_ylim(-.01, 1.01)


def load_files(config):
    # config['fit_bpz_outfile'] = 'bpz_fit_magerr_%s_%6.4f_lib_%s_z%.2f_%.2f_nwalk_%i_nsamp_%i_nstep_%i.hdf5' % (
    #     'single' if config['lg_frac_young_min'] == 0 else 'double',
    #     # config['fit_bpz_outfile'] = 'specmag_bpz_fit_magerr_%3.2f_lib_%s_z%.2f_%.2f_nwalk_%i_nsamp_%i_nstep_%i.hdf5' % (
    #     config['magerr'], os.path.splitext(os.path.basename(config['bpz_library']))[0], config['z_ini'],
    #     config['z_fin'],
    #     config['n_walkers'], config['n_samples'], config['n_steps'])
    # if config["chi2_spectrum"]:
    #     config['fit_bpz_outfile'] = "specmag_" + config['fit_bpz_outfile']

    prefix, suffix = os.path.splitext(config['fit_bpz_outfile'])
    config['fit_bpz_post'] = prefix + '_post.hdf5'
    file_fit = h5py.File('%s/%s' % (save_dir, config['fit_bpz_outfile']), 'r')
    template_magnitudes = file_fit.get('/template_magnitudes')
    _, filter_lambdas = load_filters(config)
    f_processed = h5py.File('%s/%s' % (save_dir, config['fit_bpz_post']), 'r')
    likelihood = f_processed.get('/likelihood')
    parameters = f_processed.get('/parameters')
    m2l_spec = f_processed.get('/m2l_spec')
    model_magnitudes = file_fit.get('/model_magnitudes')
    # Plot filters used for colors on template_colorspace
    filters_plot = [['F_458', 'F_768'], ['F_365', 'F_644']]
    fnumbers = np.sort(config['filters'].keys())
    plot_fid = [[int(np.argwhere(fnumbers == j)) for j in i] for i in filters_plot]
    # f_aux_templates = np.loadtxt('../templates/eB11.list', dtype='S20')
    true_values = None
    if 'bpz_library_true_values' in config:
        aux = np.loadtxt(config['bpz_library_true_values'], skiprows=1, usecols=config['bpz_library_true_values_cols'])
        true_values = np.zeros(len(aux), dtype=np.dtype([('age', float), ('metallicity', float), ('AV', float)]))
        true_values['age'] = np.log10(aux[:, 0])
        true_values['metallicity'] = np.log10(aux[:, 1])
        true_values['AV'] = aux[:, 2]
    n_z, n_templates, n_fits = likelihood.shape
    # templates = np.array(range(0, n_templates))
    templates = np.array(range(0, n_templates, config['n_interpolations'] + 1))
    n_templates = len(templates)

    return template_magnitudes, filter_lambdas, likelihood, parameters, model_magnitudes, plot_fid, true_values, n_templates, templates


f_color = {0: 'blue', 1: 'red', 2: 'green', 3: 'blue', 4: 'brown', 5: 'purple', 6: 'yellow', 7: 'grey', 8: 'red'}

if not paper:
    template_magnitudes, filter_lambdas, likelihood, parameters, model_magnitudes, plot_fid, true_values, n_templates, templates = load_files(
        config)

    # bins for the Z histogram
    metbins = np.log10(np.unique(np.sort(parameters['metallicity'])))
    metbins = np.insert(metbins, 0, metbins[0] - .3)
    metbins = np.append(metbins, metbins[-1] + .3)

    for i_template in templates:
        fig = plt.figure(1)
        fig.clf()
        axes = get_axes(fig)

        # Age historgram
        age_min = 6.4
        age_max = 10.2
        make_hist(age_min - .1, age_max + .1, 0.1, parameters['at_mass'][i_z, i_template], likelihood[i_z, i_template],
                  axes[0])
        if true_values is not None:
            axes[0].scatter(true_values['age'][i_template], 0.5, c='yellow', s=100, marker="*")
        axes[0].set_title('$\\langle \log t \\rangle_\star$')
        # axes[0].xaxis.get_ticklabels()[-1].set_visible(False)
        # axes[0].set_xlim(age_min, age_max)

        # a_v histogram
        make_hist(config['AV_min'] - .1, config['AV_max'] + .1, 0.2, parameters['a_v'][i_z, i_template],
                  likelihood[i_z, i_template], axes[1])
        if true_values is not None:
            axes[1].scatter(true_values['AV'][i_template], 0.5, c='yellow', s=100, marker="*")
        axes[1].set_title('$A_V$')
        axes[1].xaxis.get_ticklabels()[0].set_visible(False)
        axes[1].xaxis.get_ticklabels()[-1].set_visible(False)

        # M/L ratio histogram
        aux = parameters['m2l'][i_z, i_template]
        make_hist(np.nanmin(aux), np.nanmax(aux), .1, aux, likelihood[i_z, i_template], axes[2],
                  label=os.path.basename(config['filter_m2l']).split('.')[0])
        # for l in range(len(config['m2l_lambdas'])):
        #     print('l = ' + str(l))
        #     make_hist(np.nanmin(m2l_spec[i_z, i_template, :, l]) - 0.1,
        #               np.nanmax(m2l_spec[i_z, i_template, :, l]) + .1, 0.1,
        #               np.log10(m2l_spec[i_z, i_template, :, l]), likelihood[i_z, i_template], axes[2],
        #               label=str(config['m2l_lambdas'][l]), color=f_color[l], alpha=.5)
        # axes[2].legend()
        axes[2].set_title('$\log M / L$')
        axes[2].xaxis.get_ticklabels()[-1].set_visible(False)

        # Metallicity histogram
        aux1, aux2 = np.log10(0.004), np.log10(0.05)
        # np.min(np.log10(parameters['metallicity'][i_z, i_template])), np.max(
        # np.log10(parameters['metallicity'][i_z, i_template]))
        make_hist(aux1 - .3, aux2 + .3, .2, np.log10(parameters['metallicity'][i_z, i_template]),
                  likelihood[i_z, i_template], axes[3], bins=metbins)
        if true_values is not None:
            axes[3].scatter(true_values['metallicity'][i_template], 0.5, c='yellow', s=100, marker="*")
        axes[3].set_title('$\log Z_\star$')
        axes[3].xaxis.get_ticklabels()[-1].set_visible(False)

        alpha = likelihood[i_z, i_template] / (2 * likelihood[i_z, i_template].max())
        # for i_model in range(len(model_magnitudes[i_z, i_template])):
        #     axes[-1].plot(filter_lambdas, model_magnitudes[i_z, i_template, i_model] + np.mean(
        #         bpz_template_magnitudes[i_z, i_template] - model_magnitudes[i_z, i_template, i_model]), color='black',
        #                     alpha=alpha[i_model])


        # Fluxes:
        # for i_model in range(len(model_magnitudes[i_z, i_template])):
        #     axes[-1].plot(filter_lambdas, 10 ** (-.4 * (model_magnitudes[i_z, i_template, i_model] + np.mean(
        #         bpz_template_magnitudes[i_z, i_template] - model_magnitudes[i_z, i_template, i_model]))), color='black',
        #                   alpha=alpha[i_model])
        # axes[-1].plot(filter_lambdas, 10 ** (-.4 * bpz_template_magnitudes[i_z, i_template]), color='red', lw=2,
        #               alpha=0.8)
        #
        # axes[-1].plot(filter_lambdas, 10 ** (-.4 * bpz_template_magnitudes[i_z, i_template]), color='red', lw=2,
        #               alpha=0.8)
        # axes[-1].errorbar(filter_lambdas, 10 ** (-.4 * bpz_template_magnitudes[i_z, i_template]), color='red',yerr=10**(-.4*1/.05))
        # # Plot filters used for colors on template_colorspace
        # for id in np.ravel(plot_fid):
        #     axes[-1].scatter(filter_lambdas[id], 10 ** (-.4 * bpz_template_magnitudes[i_z, i_template])[id], c='yellow',
        #                      s=100, marker="*")  # , lw=2, alpha=0.8)
        # plt.ylabel('flux/arbitrary units')  # magnitudes')

        # Magnitudes:
        for i_model in range(len(model_magnitudes[i_z, i_template])):
            axes[-1].plot(filter_lambdas, model_magnitudes[i_z, i_template, i_model] + np.mean(
                template_magnitudes[i_z, i_template] - model_magnitudes[i_z, i_template, i_model]), color='black',
                          alpha=alpha[i_model])

        axes[-1].plot(filter_lambdas, template_magnitudes[i_z, i_template], color='red', lw=2, alpha=0.8)

        axes[-1].plot(filter_lambdas, template_magnitudes[i_z, i_template], color='red', lw=2, alpha=0.8)
        axes[-1].errorbar(filter_lambdas, template_magnitudes[i_z, i_template], color='red',
                          yerr=0.05)

        plt.draw()
        # Plot filters used for colors on template_colorspace
        for id in np.ravel(plot_fid):
            axes[-1].scatter(filter_lambdas[id], template_magnitudes[i_z, i_template][id], c='yellow', s=100,
                             marker="*")  # , lw=2, alpha=0.8)
        plt.draw()
        plt.ylabel('magnitudes')

        axes[-1].set_xlabel('$\\lambda [\\AA]$')
        axes[-1].invert_yaxis()
        # axes[-1].text(8000, 24.5, f_aux_templates[i_template/3].split('.')[0])
        # plt.ylim(20.5, 25)

        plt.title('template: %.2f' % (i_template / 4.))
        plt.savefig('plots/%s/tempalte_fit_%s.pdf' % (prefix, str(i_template)))

        print 'Done template # %i' % i_template
        # raw_input('Next...')

#### PLOT FOR THE PAPER ######

# err_color = {3: 'red', 1: 'green', 0.1: 'magenta', 0.5: 'cyan', 0.05: 'blue'}
err_color = {0.5: 'red', 0.1: 'green', 0.05: 'blue'}


def weighted_avg_and_std(values, weights):
    """
    From: http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy

    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)  # Fast and numerically precise
    return average, np.sqrt(variance)


stdev = dict()

if paper:

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

    i_z = 0
    i_label = 0
    fig = None
    for err in [0.5, 0.1, 0.05]:  # [0.01, 0.05, 0.1]: 0.1,
        stdev[err] = dict()
        print 'magerr = ', err
        config['magerr'] = err
        template_magnitudes, filter_lambdas, likelihood, parameters, model_magnitudes, plot_fid, true_values, n_templates, templates = load_files(
            config)
        i_ax = 0
        if not fig:
            fig = plt.figure(11, figsize=(2 * 240 / 72., 3 * 240 / 72.))
            # fig = plt.figure(1)
            fig.clf()
            fig = get_mosaic(n_templates, 5, i_fig=1)
            axes = fig.axes
        for i_template in templates:
            print 'i_template', i_template

            # For stdev versus template graph
            for k in ['at_mass', 'a_v', 'm2l', 'metallicity']:
                if not k in stdev[err]:
                    stdev[err][k] = dict()
                    stdev[err][k]['avg'] = list()
                    stdev[err][k]['stdev'] = list()
                if k in ['metallicity', 'm2l']:
                    a, s = weighted_avg_and_std(np.log10(parameters[k][i_z, i_template]), likelihood[i_z, i_template])
                else:
                    a, s = weighted_avg_and_std(parameters[k][i_z, i_template], likelihood[i_z, i_template])
                stdev[err][k]['avg'].append(a)
                stdev[err][k]['stdev'].append(s)

            # Age historgram
            age_min = 6.0
            age_max = 10.5
            make_hist(age_min, age_max, 0.1, parameters['at_mass'][i_z, i_template], likelihood[i_z, i_template],
                      axes[i_ax], color=err_color[err], alpha=.5)
            # axes[0].xaxis.get_ticklabels()[-1].set_visible(False)
            axes[0].xaxis.get_ticklabels()[0].set_visible(False)
            axes[i_ax].set_xlim(8, 10.1)
            axes[i_ax].yaxis.set_ticklabels([])

            # remove xlabels from all but the last line.
            if err != 0.05:
                axes[i_ax].xaxis.get_ticklabels()[0].set_visible(False)

            i_ax += 1
            # a_v histogram
            av_min = -.15
            av_max = 2
            make_hist(0, 2, 0.2, parameters['a_v'][i_z, i_template], likelihood[i_z, i_template], axes[i_ax],
                      color=err_color[err], alpha=.5)
            # axes[i_ax].xaxis.get_ticklabels()[0].set_visible(False)
            # axes[i_ax].xaxis.get_ticklabels()[-1].set_visible(False)
            axes[i_ax].yaxis.set_ticklabels([])
            axes[i_ax].xaxis.set_ticks([0.5, 1.0, 1.5])
            axes[i_ax].set_xlim(0, 2)

            i_ax += 1
            # M/L ratio histogram
            make_hist(-2.7, .6, 0.2, np.log10(parameters['m2l'][i_z, i_template]), likelihood[i_z, i_template],
                      axes[i_ax], color=err_color[err], alpha=.5)
            axes[i_ax].xaxis.get_ticklabels()[-1].set_visible(False)
            axes[i_ax].yaxis.set_ticklabels([])
            axes[i_ax].xaxis.set_ticks([float("%3.1f" % i) for i in np.linspace(-2.5, .4, 4)])

            i_ax += 1
            # Metallicity ratio histogram
            make_hist(-2.4, -1.31, .2, np.log10(parameters['metallicity'][i_z, i_template]), likelihood[i_z, i_template],
                      axes[i_ax], color=err_color[err], alpha=.5)
            axes[i_ax].xaxis.get_ticklabels()[-1].set_visible(False)
            axes[i_ax].yaxis.set_ticklabels([])
            axes[i_ax].xaxis.set_ticks([float("%3.1f" % i) for i in [-2.39794001, -2.09691001, -1.69897, -1.30103]])
            # np.linspace(-3.5, -1.1, 4)])

            # Fit
            i_ax += 1
            if err == 0.05:
                alpha = likelihood[i_z, i_template] / likelihood[i_z, i_template].max()
                for i_model in range(len(model_magnitudes[i_z, i_template])):
                    axes[i_ax].plot(filter_lambdas, model_magnitudes[i_z, i_template, i_model] + np.mean(
                        template_magnitudes[i_z, i_template] - model_magnitudes[i_z, i_template, i_model]),
                                    color='black',
                                    alpha=alpha[i_model])
                axes[i_ax].plot(filter_lambdas, template_magnitudes[i_z, i_template], color='red', lw=2, alpha=0.8)
                ylim = axes[i_ax].get_ylim()
                axes[i_ax].yaxis.set_ticklabels([])
                axes[i_ax].set_ylim(ylim[1], ylim[0])

                # axes[i_ax].text(8000, 24.5, f_aux_templates[i_label].split('.')[0])
                i_label += 1
                # plt.ylabel('magnitudes')
                axes[i_ax].xaxis.set_ticks([4000, 6500, 9000])

            i_ax += 1

        axes[i_ax - 5].set_xlabel('$\\langle \log t_\star \\rangle_L$')
        axes[i_ax - 4].set_xlabel('$A_V$')
        axes[i_ax - 3].set_xlabel('$\log M / L_%s$' % os.path.basename(config['filter_m2l']).split('.')[0])
        axes[i_ax - 2].set_xlabel('$\log Z_\star$')
        axes[i_ax - 1].set_xlabel('$\\lambda [\\AA]$')
    plt.savefig('plots/%s/template_fit_paper.pdf' % prefix)

i_plot = 1
plt.figure(1)
plt.clf()

err_color = {3: 'red', 1: 'green', 0.1: 'magenta', 0.5: 'cyan', 0.05: 'blue'}

for key in ['at_mass', 'a_v', 'm2l', 'metallicity']:
    plt.subplot(2, 2, i_plot)
    for err in [3, 1, 0.1, 0.5, 0.05]:
        plt.plot(np.arange(len(stdev[err][key]['avg'])) + 1, stdev[err][key]['stdev'], label='%f' % err)
    i_plot += 1
    plt.title(key)
    plt.legend(loc=2)
    plt.xlabel('i_template')
    plt.ylabel('stdev')

plt.savefig('plots/%s/stdev.pdf' % prefix)
