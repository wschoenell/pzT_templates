import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
from pystarlight.io.readfilter import readfilterfile
from pystarlight.util.base import StarlightBase
from pystarlight.util.redenninglaws import Cardelli_RedLaw
from scipy.spatial.qhull import ConvexHull
from common import load_filters, mag_in_z, calc_LsunFilter
from config import config, save_dir, prefix
from model import Model

plot_chains = False # False
plot_properties = False  # True
plot_colorspace = False
plot_inout = True


def plot_hull(c1, c2, label=None):
    points = np.array([c1, c2]).T
    hull = ConvexHull(points)
    for simplex in hull.simplices:
        if label:
            plt.plot(points[simplex, 0], points[simplex, 1], label=label,
                     color='black')  # , config['z_colors'][str(z_bpz)]
        else:
            plt.plot(points[simplex, 0], points[simplex, 1], color='black')  # , config['z_colors'][str(z_bpz)])


# Load fit:
file_fit = h5py.File('%s/%s' % (save_dir, config['fit_bpz_outfile']), 'r')
f_processed = h5py.File('%s/%s' % (save_dir, config['fit_bpz_post']), 'r')
# Create plots dir:
if not os.path.exists('plots/%s' % prefix):
    os.mkdir('plots/%s' % prefix)
# For Base plots...
basefile = None  # '/Users/william/data/BasesDir/Base.BC03.N'
if basefile is not None:
    bt = np.loadtxt(basefile, usecols=(1, 2), skiprows=1)

# For Random model plots:
random_file = config['bpz_library']


def mosaic_colorspace(i_fig):
    figure = plt.figure(i_fig)
    axis = figure.add_axes([0.1, 0.8, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.8, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.7, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.7, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.6, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.6, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.5, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.5, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.4, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.4, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.3, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.3, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.2, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.2, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis = figure.add_axes([0.1, 0.1, 0.7, 0.1])
    axis = figure.add_axes([0.8, 0.1, 0.1, 0.1], sharey=figure.axes[-1])
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    return figure


## Metallicity fix
bt_fit = StarlightBase(config['base_file'], config['base_path'])
bt_fit.ageBase[bt_fit.ageBase == 0] = 1e5
aux = np.log10(bt_fit.metBase)
aux2 = (aux[1:] - aux[0:-1]) / 2
met_min = aux - np.append(aux2[0], aux2)
met_min[0] += -0.00001
met_max = aux + np.append(aux2, aux2[-1])
met_max[-1] += 0.00001

# if random_file is not None:
#     rt = np.loadtxt(random_file, dtype=np.dtype(
#         [('fname', 'S20'), ('log t0_young', float), ('log (tau_young / t0_young)', float), ('log t0_old', float),
#          ('log (tau_old / t0_old)', float),
#          ('Mfrac_young', float), ('a_v', float), ('metallicity', float)]))
#
#     rt['log (tau_young / t0_young)'] = np.log10(rt['log (tau_young / t0_young)'] / rt['log t0_young'])
#     rt['log (tau_old / t0_old)'] = np.log10(rt['log (tau_old / t0_old)'] / rt['log t0_old'])
#     rt['Mfrac_young'] = rt['Mfrac_young'] # np.log10(rt['Mfrac_young'])
#     rt['log t0_young'] = np.log10(rt['log t0_young'])
#     rt['log t0_old'] = np.log10(rt['log t0_old'])
#     rt['metallicity'] = np.log10(rt['metallicity'])
#
#     i_norm = np.argmin((bt_fit.l_ssp - config['lambda_norm']) ** 2)
#     l2m_norm_ssp = np.empty((len(bt_fit.metBase), len(bt_fit.ageBase)))
#     for i_met in range(len(bt_fit.metBase)):
#         for i_age in range(len(bt_fit.ageBase)):
#             l2m_norm_ssp[i_met, i_age] = bt_fit.f_ssp[i_met, i_age, i_norm]

# For Random model plots:
if random_file is not None:
    rt = np.loadtxt(random_file, dtype=np.dtype(
        [('fname', 'S20'), ('log t0_young', float), ('log (tau_young / t0_young)', float), ('log t0_old', float),
         ('log (tau_old / t0_old)', float),
         ('Mfrac_young', float), ('a_v', float), ('metallicity', float)]))

    rt['log (tau_young / t0_young)'] = np.log10(rt['log (tau_young / t0_young)'] / rt['log t0_young'])
    rt['log (tau_old / t0_old)'] = np.log10(rt['log (tau_old / t0_old)'] / rt['log t0_old'])
    rt['Mfrac_young'] = rt['Mfrac_young']
    rt['log t0_young'] = np.log10(rt['log t0_young'])
    rt['log t0_old'] = np.log10(rt['log t0_old'])
    rt['metallicity'] = np.log10(rt['metallicity'])

    i_norm = np.argmin((bt_fit.l_ssp - config['lambda_norm']) ** 2)
    l2m_norm_ssp = np.empty((len(bt_fit.metBase), len(bt_fit.ageBase)))
    for i_met in range(len(bt_fit.metBase)):
        for i_age in range(len(bt_fit.ageBase)):
            l2m_norm_ssp[i_met, i_age] = bt_fit.f_ssp[i_met, i_age, i_norm]

# Plot chains
full_chi2 = -2 * file_fit['full_lnprob'].value / len(config['filters'])
i_z = 0
param_names = ['log t0_young', 'log (tau_young / t0_young)', 'log t0_old', 'log (tau_old / t0_old)', 'Mfrac_young',
               'a_v', 'metallicity']
param_priors = [[np.log10(config['t0_young_min']), np.log10(config['t0_young_max'])],
                [np.log10(config['tau_t0_young_min']), np.log10(config['tau_t0_young_max'])],
                [np.log10(config['t0_old_min']), 10.2],
                [np.log10(config['tau_t0_old_min']), np.log10(config['tau_t0_old_max'])],
                [-.1, 1.1],
                [config['AV_min'], config['AV_max']],
                [np.min(file_fit['full_chain'][:, :, :, :, -1]) - .2,
                 np.max(file_fit['full_chain'][:, :, :, :, -1]) + .2]]

n_templates = file_fit['full_chain'].shape[1]
n_chains = file_fit['full_chain'].shape[2]
n_steps = file_fit['full_chain'].shape[3]

if plot_chains:
    plt.figure(1, figsize=(11.69, 8.27))
    for i_template in range(n_templates):
        metallicity = []
        for met_chains in file_fit['full_chain'][i_z, i_template, :, :, 6]:
            metallicity.append(
                [np.log10(bt_fit.metBase[np.argwhere(np.bitwise_and(met >= met_min, met <= met_max))[0, 0]]) for met in
                 met_chains])
        metallicity = np.array(metallicity)
        # sys.exit()
        plt.clf()
        fig = mosaic_colorspace(i_fig=2)
        axes = fig.axes
        i_axis = 0
        for i_param in range(len(param_names)):
            ax = axes[i_axis]
            i_axis += 1
            # if 'tau' in param_names[i_param]:
            #     ax.plot(range(file_fit['full_chain'].shape[2]), np.log10(file_fit['full_chain'][i_z, i_template, :, i_param]))
            # else:


            # Plot chains
            for i_chain in range(n_chains):
                if 'metallicity' in param_names[i_param]:
                    ax.plot(range(n_steps), metallicity[i_chain], color='black', alpha=0.3)
                elif 'Mfrac' in param_names[i_param]:
                    ax.plot(range(n_steps), 10 ** file_fit['full_chain'][i_z, i_template, i_chain, :, i_param],
                            color='black', alpha=0.3)
                else:
                    ax.plot(range(n_steps), file_fit['full_chain'][i_z, i_template, i_chain, :, i_param], color='black',
                            alpha=0.3)

            # Plot "correct" value, if exists for SSPs
            if basefile is not None:
                if 'log t0' in param_names[i_param]:
                    print param_names[i_param]
                    ax.plot(range(n_steps), np.ones(n_steps) * np.log10(bt[i_template][0]), color='red', alpha=0.5)
                ax.set_ylim(param_priors[i_param][0] * .9, param_priors[i_param][1] * 1.1)
                if 'metallicity' in param_names[i_param]:
                    print param_names[i_param]
                    ax.plot(range(n_steps), np.ones(n_steps) * np.log10(bt[i_template][1]), color='red', alpha=0.5)
                    # ax.set_ylim(param_priors[i_param][0]*1.2, param_priors[i_param][1]*1.2)
            if random_file is not None:
                ax.plot(range(n_steps), np.ones(n_steps) * rt[param_names[i_param]][i_template], color='red', alpha=0.5)

            lim = ax.get_xlim()
            xpos = lim[0] + (lim[1] - lim[0]) * .6
            lim = ax.get_ylim()
            ypos = lim[0] + (lim[1] - lim[0]) * .8
            ax.text(xpos, ypos, param_names[i_param])

            ax = axes[i_axis]
            i_axis += 1
            if 'metallicity' in param_names[i_param]:
                ax.hist(np.ravel(metallicity),
                        range=[param_priors[i_param][0], param_priors[i_param][1]], orientation='horizontal')
            elif 'Mfrac' in param_names[i_param]:
                ax.hist(10 ** np.ravel(file_fit['full_chain'][i_z, i_template, :, :, i_param]),
                        range=[param_priors[i_param][0], param_priors[i_param][1]], orientation='horizontal')
            else:
                ax.hist(np.ravel(file_fit['full_chain'][i_z, i_template, :, :, i_param]),
                        range=[param_priors[i_param][0], param_priors[i_param][1]], orientation='horizontal')

        # Chi2 plot
        # chi2_range = np.percentile(-2*full_lnprob[i_z, i_template, i_chain], [20, 100])
        chi2_range = np.array([0, 2])
        ax = axes[-2]
        for i_chain in range(n_chains):
            ax.plot(range(n_steps), full_chi2[i_z, i_template, i_chain], color='black', alpha=0.3)
        ax.set_ylim(tuple(chi2_range))

        # Plot label:
        lim = ax.get_xlim()
        xpos = lim[0] + (lim[1] - lim[0]) * .6
        lim = ax.get_ylim()
        ypos = lim[0] + (lim[1] - lim[0]) * .8
        ax.text(xpos, ypos, 'chi2')

        # right histogram
        ax = axes[-1]
        ax.hist(np.ravel(full_chi2[i_z, i_template]), range=chi2_range, orientation='horizontal')

        # Title and save
        if basefile is not None:
            axes[0].set_title('%.2f %.2f acceptance: %.2f' % (np.log10(bt[i_template][0]), np.log10(bt[i_template][1]),
                                                              np.mean(
                                                                  file_fit['acceptance_fraction'][i_z, i_template, :])))
        print('full_chain: %i' % i_template)
        plt.savefig('plots/%s/fullchain_%s.pdf' % (prefix, str(i_template)))

if plot_properties:
    figure = plt.figure(2, figsize=(11.69, 8.27))
    plt.clf()
    for i_template in range(n_templates):

        plt.clf()

        ax_hist_upp_left = figure.add_axes([0.25, 0.75, 0.325, 0.15])
        ax_hist_upp_right = figure.add_axes([0.575, 0.75, 0.325, 0.15])
        ax_hist_low = figure.add_axes([0.1, 0.1, 0.15, 0.65])
        ax_left_low = figure.add_axes([0.25, 0.1, 0.325, 0.65], sharey=figure.get_axes()[-1])
        figure.get_axes()[-1].yaxis.set_visible(False)
        ax_right_low = figure.add_axes([0.575, 0.1, 0.325, 0.65], sharey=figure.get_axes()[-1])
        figure.get_axes()[-1].yaxis.set_visible(False)

        # i_template = 0
        # get_sfh params = bt_fit, t0_young, tau_young, t0_old, tau_old, frac_young
        # model_params = t0_young, tauT_young, t0_old, tauT_old, frac_young, a_v, metallicity

        # plt.subplot(1, 2, 1)  # a_v vs age
        # plt.plot(f_processed['parameters']['a_v'][i_z, i_template], f_processed['parameters']['at_mass'][i_z, i_template],
        #          '.', color='black', alpha=.3)
        ax_left_low.scatter(f_processed['parameters']['a_v'][i_z, i_template],
                            f_processed['parameters']['at_mass'][i_z, i_template],
                            c=np.log10(f_processed['parameters']['metallicity'][i_z, i_template]), alpha=.3)
        # plt.colorbar()


        if basefile is not None:
            plt.plot(0, np.log10(bt[i_template][0]), '.', color='red')
            plt.errorbar(0, np.log10(bt[i_template][0]), xerr=0.1, yerr=0.2, color='red')
        if random_file is not None:
            # t0_young, tau_young, t0_old, tau_old, frac_youn
            sfh = Model.get_sfh(bt_fit, 10 ** rt['log t0_young'][i_template],
                                10 ** rt['log t0_young'][i_template] * 10 ** rt['log (tau_young / t0_young)'][
                                    i_template],
                                10 ** rt['log t0_old'][i_template],
                                10 ** rt['log t0_old'][i_template] * 10 ** rt['log (tau_old / t0_old)'][i_template],
                                rt['Mfrac_young'][i_template])

            i_met = np.argwhere(
                np.bitwise_and(rt['metallicity'][i_template] >= met_min, rt['metallicity'][i_template] <= met_max))[
                0, 0]
            age = np.sum(sfh * l2m_norm_ssp[i_met] * np.log10(bt_fit.ageBase)) / np.sum(sfh * l2m_norm_ssp[i_met])
            ax_left_low.plot(rt['a_v'][i_template], age, '.', color='red')
            ax_left_low.errorbar(rt['a_v'][i_template], age, xerr=0.1, yerr=0.2, color='red')

        ax_left_low.set_xlabel('AV')

        # plt.subplot(1, 2, 2)  # age vs m2l
        ax_right_low.scatter(np.log10(np.ravel(f_processed['m2l_spec'][:, :, :, 0][i_z, i_template])),
                             f_processed['parameters']['at_mass'][i_z, i_template],
                             c=np.log10(f_processed['parameters']['metallicity'][i_z, i_template]), alpha=.3)

        ax_right_low.plot(np.log10(np.average(1. / l2m_norm_ssp[i_met, :], weights=sfh)), age, '.', color='red')
        ax_right_low.errorbar(np.log10(np.average(1. / l2m_norm_ssp[i_met, :], weights=sfh)), age, xerr=0.1, yerr=0.1,
                              color='red')

        ax_right_low.set_xlabel('log m/l 4020')

        bins = np.arange(np.min(f_processed['parameters']['at_mass'][i_z, i_template]),
                         np.min(f_processed['parameters']['at_mass'][i_z, i_template]), 0.1)
        ax_hist_low.hist(f_processed['parameters']['at_mass'][i_z, i_template], normed=True, bins=20,
                         orientation='horizontal', alpha=.3, color='black')
        ax_hist_low.hist(f_processed['parameters']['at_mass'][i_z, i_template],
                         weights=f_processed['likelihood'][i_z, i_template], normed=True, bins=20,
                         orientation='horizontal', alpha=.3, color='red')
        ax_hist_low.invert_xaxis()
        ax_hist_low.set_ylabel('<age>')

        # plt.subplot(1, 2, 2)        # metallicity vs age
        # plt.plot(f_processed['parameters']['metallicity'][i_z, i_template], f_processed['parameters']['at_mass'][i_z, i_template],
        #          '.', color='black', alpha=.3)
        #
        # # if basefile is not None:
        # #     plt.plot(0, np.log10(bt[i_template][0]), '.', color='red')
        # #     plt.errorbar(0, np.log10(bt[i_template][0]), xerr=0.1, yerr=0.2, color='red')
        # # if random_file is not None:
        # #     # t0_young, tau_young, t0_old, tau_old, frac_youn
        # #     sfh = Model.get_sfh(bt_fit, 10 ** rt['log t0_young'][i_template],
        # #                         10 ** rt['log t0_young'][i_template] * 10 ** rt['log (tau_young / t0_young)'][
        # #                             i_template],
        # #                         10 ** rt['log t0_old'][i_template],
        # #                         10 ** rt['log t0_old'][i_template] * 10 ** rt['log (tau_old / t0_old)'][i_template],
        # #                         rt['Mfrac_young'][i_template])
        # #
        # #     i_met = np.argwhere(
        # #         np.bitwise_and(rt['metallicity'][i_template] >= met_min, rt['metallicity'][i_template] <= met_max))[
        # #         0, 0]
        # #     age = np.sum(sfh * l2m_norm_ssp[i_met] * np.log10(bt_fit.ageBase)) / np.sum(sfh * l2m_norm_ssp[i_met])
        # #     plt.plot(rt['a_v'][i_template], age, '.', color='red')
        # #     plt.errorbar(rt['a_v'][i_template], age, xerr=0.1, yerr=0.2, color='red')
        #
        # plt.xlabel('Z')
        # plt.ylabel('<age>')


        plt.savefig('plots/%s/properties_%s.pdf' % (prefix, str(i_template)))
        print('properties: %f' % i_template)

# convexhull
if plot_colorspace:
    plt.figure(3, figsize=(11.69, 8.27))
    plt.clf()

    filters_plot = [['F_458', 'F_768'], ['F_365', 'F_644']]
    scatter_color = 'metallicity'

    i_z = 0

    fnumbers = np.sort(config['filters'].keys())
    fid = [[int(np.argwhere(fnumbers == j)) for j in i] for i in filters_plot]

    # file_fit = h5py.File('bpz_fit_full_newmask_kk_test.hdf5', 'r')
    # For the A_V Arrow:
    e_c1 = Cardelli_RedLaw([4580])[0] - Cardelli_RedLaw([7680])[0]
    e_c2 = Cardelli_RedLaw([3650])[0] - Cardelli_RedLaw([6440])[0]

    # Open model
    likelihood = f_processed.get('/likelihood')
    if scatter_color == 'metallicity':
        parameter_color = np.log10(f_processed.get('/parameters')[scatter_color])
    else:
        parameter_color = f_processed.get('/parameters')[scatter_color]

    vmin, vmax = np.min(parameter_color), np.max(parameter_color)  # For a constant colorbar

    # Open BPZ magnitudes
    templates = np.loadtxt(config['bpz_library'], np.str)
    aux = templates[:, 0] if len(templates.shape) > 1 else templates
    templates_data = [
        np.loadtxt('%s/%s' % (config['bpz_library_dir'], t), dtype=np.dtype([('lambda', np.float), ('flux', np.float)]))
        for
        t in aux]
    interpolations = np.linspace(0., len(templates) - 1,
                                 len(templates) + (len(templates) - 1) * config['n_interpolations'])

    # Plot model magnitudes convex hull
    # FIXME: model_magnitudes = file_fit['model_magnitudes']


    l = []
    for i in range(file_fit['model_magnitudes'].shape[1]):
        l += [i]
        # l += [None, i]

    for plot_template in l:
        # plot_template = 36
        # plt.cm
        plt.clf()
        if plot_template is None:
            plt.clf()
        else:
            # For the Age Arrow:
            ages_base = np.log10(bt_fit.ageBase)
            filters, _ = load_filters(config)
            i_met = int(np.argwhere(bt_fit.metBase == (bt[:, 1][plot_template] if basefile is not None else 0.02)))
            mags = np.array(
                [mag_in_z(bt_fit.l_ssp, bt_fit.f_ssp[i_met, i_model], config['z_ini'], filters) for i_model in
                 range(bt_fit.nAges)])
            color1_base = np.ravel(mags[:, fid[0][0]] - mags[:, fid[0][1]])
            color2_base = np.ravel(mags[:, fid[1][0]] - mags[:, fid[1][1]])
        # model_magnitudes = np.load('model_magnitudes.npz')['arr_0']
        model_magnitudes = np.copy(file_fit['model_magnitudes'])
        template_magnitudes = np.copy(file_fit['template_magnitudes'])
        n_z, n_templates, n_fits, n_filters = model_magnitudes.shape
        if plot_template is not None:
            aux_shape = list(model_magnitudes.shape)
            aux_shape[1] = 1
            model_magnitudes = model_magnitudes[:, plot_template, :, :]
            model_magnitudes = model_magnitudes.reshape(aux_shape)
            aux_shape = list(template_magnitudes.shape)
            aux_shape[1] = 1
            template_magnitudes = template_magnitudes[:, plot_template, :]
            template_magnitudes = template_magnitudes.reshape(aux_shape)
            interpolations_colors = [interpolations[plot_template]]
            point_color = np.array([interpolations[plot_template]] * n_fits)
        else:
            interpolations_colors = interpolations
            # point_color =[]
            # for i_template in range(n_templates):
            #     point_color += [interpolations[i_template]]*n_fits
            # point_color = np.array(point_color)
            point_color = 'gray'
        n_z, n_templates, n_fits, n_filters = model_magnitudes.shape
        color1 = np.ravel(model_magnitudes[i_z, :, :, fid[0][0]] - model_magnitudes[i_z, :, :, fid[0][1]])
        color2 = np.ravel(model_magnitudes[i_z, :, :, fid[1][0]] - model_magnitudes[i_z, :, :, fid[1][1]])

        if plot_template is None:
            plot_hull(color1, color2)
            plt.scatter(color1, color2, c=point_color, marker='.', alpha=0.2)
        else:
            plt.scatter(color1, color2, c=parameter_color[i_z, plot_template, :], alpha=.7, vmin=vmin,
                        vmax=vmax)  # , marker='.', alpha=1)
            plt.colorbar()

        # plt.plot(color1, color2, '.')
        # Plot (BPZ) template magnitudes
        color1 = np.ravel(template_magnitudes[i_z, :, fid[0][0]] - template_magnitudes[i_z, :, fid[0][1]])
        color2 = np.ravel(template_magnitudes[i_z, :, fid[1][0]] - template_magnitudes[i_z, :, fid[1][1]])
        if plot_template is None:
            plt.scatter(color1, color2, c=interpolations_colors, s=40, marker='*')
        else:
            plt.scatter(color1, color2, c='red', s=80, marker='*')
            plt.errorbar(color1, color2, color='red', xerr=0.07, yerr=0.07)
            # A_V Arrow:
            plt.arrow(color1[0], color2[0], dx=e_c1 * 0.5, dy=e_c2 * 0.5)
            # Age Arrow:
            if basefile is not None:
                ages_templates = np.log10(bt[:, 0])
                x, y = np.interp(np.log10(bt[:, 0])[plot_template], ages_base, color1_base), np.interp(
                    ages_templates[plot_template], ages_base, color2_base)
                dx, dy = np.interp(ages_templates[plot_template] + 0.1, ages_base, color1_base) - x, np.interp(
                    ages_templates[plot_template] + 0.1, ages_base, color2_base) - y
                plt.arrow(x, y, dx, dy, color='red')
                dx, dy = np.interp(ages_templates[plot_template] - 0.1, ages_base, color1_base) - x, np.interp(
                    ages_templates[plot_template] - 0.1, ages_base, color2_base) - y
                plt.arrow(x, y, dx, dy, color='red')
            plt.plot(color1_base, color2_base, color='blue')

            plt.xlim(color1[0] - .5, color1[0] + .5)
            plt.ylim(color2[0] - .5, color2[0] + .5)
        # for i in range(len(interpolations)):
        #     plt.text(color1[i], color2[i], interpolations[i], color='blue')
        # plt.colorbar()

        plt.xlabel("%s - %s" % (filters_plot[0][0], filters_plot[0][1]))

        plt.ylabel("%s - %s" % (filters_plot[1][0], filters_plot[1][1]))

        if plot_template is not None:
            plt.title('%.2f' % interpolations[plot_template])
        plt.savefig('plots/%s/%s_colorspace_%s.png' % (prefix, scatter_color, str(plot_template)))
        print 'Done %s.' % str(plot_template)

        # raw_input('Next...')

# file_fit.close()


# Plot age-color for solar metallicity
if basefile is not None:
    plt.figure(3)
    plt.clf()
    filters, _ = load_filters(config)
    i_met = int(np.argwhere(bt_fit.metBase == 0.02))
    mags = np.array(
        [mag_in_z(bt_fit.l_ssp, bt_fit.f_ssp[i_met, i_model], config['z_ini'], filters) for i_model in
         range(bt_fit.nAges)])
    color1 = np.ravel(mags[:, fid[0][0]] - mags[:, fid[0][1]])
    color2 = np.ravel(mags[:, fid[1][0]] - mags[:, fid[1][1]])
    plt.plot(color1, color2)
    ages_base = np.log10(bt_fit.ageBase)
    ages = np.log10(bt[:, 0])
    for age in ages:
        x, y = np.interp(age, ages_base, color1), np.interp(age, ages_base, color2)
        dx, dy = np.interp(age + 0.1, ages_base, color1) - x, np.interp(age + 0.1, ages_base, color2) - y
        plt.arrow(x, y, dx, dy, color='red')
    plt.xlabel("%s - %s" % (filters_plot[0][0], filters_plot[0][1]))
    plt.ylabel("%s - %s" % (filters_plot[1][0], filters_plot[1][1]))

if plot_inout:

    # aux_m2l_lambdas_i = [int(np.argwhere(bt_fit.l_ssp == i)) for i in config['m2l_lambdas']]
    # m2l_ssp_spec = np.zeros((len(bt_fit.metBase), len(bt_fit.ageBase), len(config["m2l_lambdas"])))
    # for i_met in range(len(bt_fit.metBase)):
    #     for i_age in range(len(bt_fit.ageBase)):
    #         m2l_ssp_spec[i_met, i_age] = 1. / bt_fit.f_ssp[i_met, i_age, aux_m2l_lambdas_i]
    # bt_fit.Mstars[i_met, i_age] / bt_fit.f_ssp[i_met, i_age, aux_m2l_lambdas_i]

    filter_m2l = readfilterfile(config['filter_m2l'], ).data
    aux_l = np.arange(filter_m2l['lambda'].min(), filter_m2l['lambda'].max())
    filter_m2l_new = np.empty(len(aux_l), dtype=filter_m2l.dtype)
    filter_m2l_new['lambda'] = aux_l
    filter_m2l_new['transm'] = np.interp(aux_l, filter_m2l['lambda'], filter_m2l['transm'])
    filter_m2l = filter_m2l_new
    LsunFilter = calc_LsunFilter(filter_m2l)

    l2m_ssp = np.empty((len(bt_fit.metBase), len(bt_fit.ageBase)))
    l2m_ssp_spec = np.zeros((len(bt_fit.metBase), len(bt_fit.ageBase), len(config["m2l_lambdas"])))
    for i_met in range(len(bt_fit.metBase)):
        for i_age in range(len(bt_fit.ageBase)):
            aux = np.interp(filter_m2l['lambda'], bt_fit.l_ssp, bt_fit.f_ssp[i_met, i_age])
            l2m_ssp[i_met, i_age] = 1 / (bt_fit.Mstars[i_met, i_age] * LsunFilter / np.trapz(aux * filter_m2l['transm'],
                                                                                             filter_m2l['lambda']))

    plt.figure(4)
    # at_flux:
    at_flux = np.average(f_processed['parameters']['at_flux'][i_z, :], weights=f_processed['likelihood'][i_z, :],
                         axis=1)
    # at_mass:
    at_mass = np.average(f_processed['parameters']['at_mass'][i_z, :], weights=f_processed['likelihood'][i_z, :],
                         axis=1)
    mean_lgmetallicity = np.average(np.log10(f_processed['parameters']['metallicity'][i_z, :]),
                                    weights=f_processed['likelihood'][i_z, :], axis=1)

    # mean_m2l = np.array([np.average(f_processed['m2l_spec'][i_z, :, :, i], weights=f_processed['likelihood'][i_z],
    #                                 axis=1) for i in range(len(config["m2l_lambdas"]))])

    mean_m2l = np.average(f_processed['parameters']['m2l'][i_z], weights=f_processed['likelihood'][i_z], axis=1)

    std_m2l = np.sqrt(
        np.average((np.log10(f_processed['parameters']['m2l'][i_z]) - np.log10(mean_m2l[:, np.newaxis])) ** 2,
                   weights=f_processed['likelihood'][i_z], axis=1))

    mean_av = np.average(f_processed['parameters']['a_v'][i_z, :], weights=f_processed['likelihood'][i_z, :], axis=1)
    real_at_flux = list()
    real_at_mass = list()
    real_m2l = list()
    for i_template in range(n_templates):
        sfh = Model.get_sfh(bt_fit, 10 ** rt['log t0_young'][i_template],
                            10 ** rt['log t0_young'][i_template] * 10 ** rt['log (tau_young / t0_young)'][i_template],
                            10 ** rt['log t0_old'][i_template],
                            10 ** rt['log t0_old'][i_template] * 10 ** rt['log (tau_old / t0_old)'][i_template],
                            rt['Mfrac_young'][i_template])

        i_met = np.argwhere(
            np.bitwise_and(rt['metallicity'][i_template] >= met_min, rt['metallicity'][i_template] <= met_max))[
            0, 0]
        # at_flux
        real_at_flux.append(
            np.sum(sfh * l2m_norm_ssp[i_met] * np.log10(bt_fit.ageBase)) / np.sum(sfh * l2m_norm_ssp[i_met]))
        # at_mass:
        real_at_mass.append(np.sum(sfh * np.log10(bt_fit.ageBase)) / np.sum(sfh))
        # real_m2l.append([np.average(m2l_ssp_spec[i_met, :, i], weights=sfh) for i in range(len(config["m2l_lambdas"]))])
        real_m2l.append(1 / np.sum(l2m_ssp[i_met] * sfh))

    real_at_flux = np.array(real_at_flux)
    real_at_mass = np.array(real_at_mass)
    real_m2l = np.array(real_m2l)

    mask = np.ones(len(real_at_flux), dtype=np.bool)
    # mask = np.bitwise_and()
    # mask = np.where(real_at_flux > 8)

    plt.clf()

    # at_flux
    # plt.subplot(4, 2, 1)
    plt.title('magerr = %f - magnoise = %f' % (config['magerr'], config['mag_noise']))

    plt.subplot(4, 2, 1)
    plt.scatter(real_at_flux[mask], at_flux[mask], c=rt['metallicity'][mask], alpha=.3)
    plt.plot([-4, 11], [-4, 11])
    plt.xlabel('age_real')
    plt.ylabel('age_mean')
    plt.xlim(5.6, 10.2)
    plt.ylim(5.6, 10.2)
    plt.grid()

    plt.subplot(4, 2, 2)
    diff = at_flux[mask] - real_at_flux[mask]
    plt.scatter(real_at_flux[mask], diff, alpha=.3, c=rt['metallicity'][mask])  # color='black')
    plt.plot([6, 10.5], np.ones(2) * np.average(diff), color='red')
    [plt.plot([6, 10.5], np.ones(2) * x, color='black') for x in np.percentile(diff, [16, 84])]
    plt.plot([6, 10.5], np.ones(2) * np.average(diff) - np.std(diff), color='magenta')
    plt.plot([6, 10.5], np.ones(2) * np.average(diff) + np.std(diff), color='magenta')
    plt.plot([5.5, 10.5], [.2, .2], color='blue')
    plt.plot([5.5, 10.5], [-.2, -.2], color='blue')
    plt.xlim(5.6, 10.2)
    plt.ylim(-.5, .5)
    plt.xlabel('age_real')
    plt.ylabel('mean - real')
    plt.grid()
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.text(xmin + (xmax - xmin) / 10, ymax - (ymax - ymin) / 10,
             'd = %3.2f s = %3.2f' % (np.mean(diff), np.std(diff)))

    # at_mass
    plt.subplot(4, 2, 3)
    plt.scatter(real_at_mass[mask], at_mass[mask], c=rt['metallicity'][mask], alpha=.3)
    plt.plot([-4, 11], [-4, 11])
    plt.xlabel('at_mass real')
    plt.ylabel('at_mass mean')
    plt.xlim(5.6, 10.2)
    plt.ylim(5.6, 10.2)
    plt.grid()

    plt.subplot(4, 2, 4)
    diff = at_mass[mask] - real_at_mass[mask]
    plt.scatter(real_at_mass[mask], diff, alpha=.3, c=rt['metallicity'][mask])  # color='black')
    plt.plot([6, 10.5], np.ones(2) * np.average(diff), color='red')
    [plt.plot([6, 10.5], np.ones(2) * x, color='black') for x in np.percentile(diff, [16, 84])]
    plt.plot([6, 10.5], np.ones(2) * np.average(diff) - np.std(diff), color='magenta')
    plt.plot([6, 10.5], np.ones(2) * np.average(diff) + np.std(diff), color='magenta')
    plt.plot([5.5, 10.5], [.2, .2], color='blue')
    plt.plot([5.5, 10.5], [-.2, -.2], color='blue')
    plt.xlim(5.6, 10.2)
    plt.ylim(-.5, .5)
    plt.xlabel('at_mass real')
    plt.ylabel('mean - real')
    plt.grid()
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.text(xmin + (xmax - xmin) / 10, ymax - (ymax - ymin) / 10,
             'd = %3.2f s = %3.2f' % (np.mean(diff), np.std(diff)))

    plt.subplot(4, 2, 7)
    plt.scatter(rt['a_v'][mask], mean_av[mask], c=rt['metallicity'][mask], alpha=.3)
    plt.plot([-.1, 2], [-.1, 2])
    plt.xlabel('AV_real')
    plt.ylabel('AV_mean')
    # plt.xlim(5, 11)
    # plt.ylim(5, 11)
    plt.grid()

    plt.subplot(4, 2, 8)
    diff = mean_av[mask] - rt['a_v'][mask]
    plt.scatter(rt['a_v'][mask], diff, alpha=.3, color='black')
    plt.plot([-.1, 2], np.ones(2) * np.average(diff), color='red')
    [plt.plot([-.1, 2], np.ones(2) * x, color='black') for x in np.percentile(diff, [16, 84])]
    plt.plot([-.1, 2], np.ones(2) * np.average(diff) - np.std(diff), color='magenta')
    plt.plot([-.1, 2], np.ones(2) * np.average(diff) + np.std(diff), color='magenta')
    # plt.xlim(5, 11)
    # plt.ylim(-1, 1)
    plt.xlabel('AV_real')
    plt.ylabel('mean - real')
    plt.grid()
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.text(xmin + (xmax - xmin) / 10, ymax - (ymax - ymin) / 10,
             'd = %3.2f s = %3.2f' % (np.mean(diff), np.std(diff)))

    plt.subplot(4, 2, 6)

    diff = np.log10(mean_m2l)[mask] - np.log10(real_m2l)[mask]
    plt.scatter(np.log10(real_m2l)[mask], diff, alpha=.3, c=rt['a_v'][mask])  # color='black')
    # plt.xlim(5, 11)
    # plt.ylim(-1, 1)
    plt.plot([-2, 2], np.ones(2) * np.average(diff), color='red')
    [plt.plot([-2, 2], np.ones(2) * x, color='black') for x in np.percentile(diff, [16, 84])]
    plt.plot([-2, 2], np.ones(2) * np.average(diff) - np.std(diff), color='magenta')
    plt.plot([-2, 2], np.ones(2) * np.average(diff) + np.std(diff), color='magenta')
    plt.xlim(-2, .6)
    plt.xlabel('log M2L_real')
    plt.ylabel('log (mean) - log (real)')
    plt.grid()
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.text(xmin + (xmax - xmin) / 10, ymax - (ymax - ymin) / 10,
             'd = %3.2f s = %3.2f' % (np.mean(diff), np.std(diff)))

    plt.subplot(4, 2, 5)
    # plt.scatter(np.log10(mean_m2l)[mask], std_m2l[mask], alpha=.3, color='black')
    plt.scatter(np.log10(mean_m2l)[mask], np.log10(real_m2l[mask]), alpha=.3, color='black')
    plt.xlabel('log M2L_real')
    plt.plot([-.4, .6], [-.4, .6], color='blue')
    # plt.ylabel('log stdev ')
    plt.ylabel('log mean')
    plt.grid()

    plt.savefig('plots/%s/inout.pdf' % prefix)

    plt.figure(5)
    plt.title('magerr = %f - magnoise = %f' % (config['magerr'], config['mag_noise']))
    plt.clf()
    diff = rt['metallicity'] - mean_lgmetallicity
    plt.scatter(rt['metallicity'], diff, alpha=.3, c=real_at_flux)  # color='black')
    # plt.xlim(5, 11)
    # plt.ylim(-1, 1)
    plt.plot([-2.5, -1], np.ones(2) * np.average(diff), color='red')
    [plt.plot([-2.5, -1], np.ones(2) * x, color='black') for x in np.percentile(diff, [16, 84])]
    plt.plot([-2.5, -1], np.ones(2) * np.average(diff) - np.std(diff), color='magenta')
    plt.plot([-2.5, -1], np.ones(2) * np.average(diff) + np.std(diff), color='magenta')
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.text(xmin + (xmax - xmin) / 10, ymax - (ymax - ymin) / 10,
             'd = %3.2f s = %3.2f' % (np.mean(diff), np.std(diff)))
    # plt.xlim(-1, 1)

    plt.savefig('plots/%s/inout_metallicity.pdf' % prefix)
