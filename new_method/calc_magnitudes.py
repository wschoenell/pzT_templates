model_magnitudes = np.empty((z_len, templates_len, n_samples, magnitudes_len))
for i_z in z_fit:
    print 'Evaluating magnitudes for z = %.2f' % z_fit[i_z]
    for i_t in range(templates_len):
        model = Model(templates_magnitudes[i_z, i_t], filters, z[i_z], config['base_file'], config['base_path'],
                      config['taylor_file'], config['magerr'])
        for i_sample, (
        lg_t0_young, lg_tauT_young, lg_t0_old, lg_tauT_old, lg_frac_young, a_v, lg_metallicity) in enumerate(
                model_parameters[i_z, i_t]):
            # print i_sample, i_z, i_t
            aux_mags = model.get_mags(10 ** lg_t0_young, 10 ** (lg_tauT_young + lg_t0_young),
                                      10 ** lg_t0_old, 10 ** (lg_tauT_old + lg_t0_old),
                                      10 ** lg_frac_young, a_v, 10 ** lg_metallicity)
            model_magnitudes[i_z, i_t, i_sample] = aux_mags