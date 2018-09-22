
class make_config():

    def __init__(self, dict_files, section):

        self.name = name
        self.bad_frames = bad_frames  # user defined
        self.dark_frames = dark_frames

        # number of filters will decide number of parameters here
        self.filters = filters

        # sky subtraction
        self.sky_min = sky_min  # user defined
        self.sky_max = sky_max  # user defined

        # dark parameters
        self.aper_lim = aper_lim
        self.cent_x = cent_x
        self.cent_y = cent_y
        self.dark_radii = dark_radii  # user defined

        # xregister parameters UV

        self.xref_ref = xreg_ref  # user defined
        self.xreg_lim = xreg_lim  # user defined
        self.xreg_xwindow = xreg_xwindow  # user defined
        self.xreg_ywindow = xreg_ywindow  # user defined

        # xregister parameters HA

        self.xwin_ha = xwin_ha  # user defined
        self.ywin_ha = ywin_ha  # user defined


        # add PSF parameters
if __name__ == '__main__':
        # if not os.path.exists(config_file):

    config_file = config_dir + 'test.cfg'
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())

    config.add_section('basic')
    config.set('basic', 'work_dir', data_dir)
    config.set('basic', 'reporsitory_dir', '/home/sourabh/ULIRG_v2/')
    config.set('basic', 'psf_dir', '${reporsitory_dir}PSF/')
    config.set('basic', 'dark_dir', '${reporsitory_dir}DARK_files_ISR/DARK_ULIRG/')
    config.set('basic', 'flt_files', '${reporsitory_dir}scripts/FLT_list.txt')
    config.set('basic', 'dq', '1000')
    config.set('basic', 'dq_bad', '0.0')

    # sky param
    config.set('basic', 'sky_min', '400')
    config.set('basic', 'sky_max', '480')

    # <<<<<<<<<<<<<<<<<<<<<DARK parameters >>>>>>>>>>>>>>>>>>>>>
    # sextractor config file
    config.set('basic', 'sex_config', '${reporsitory_dir}scripts/sex_default.conf')

    ###
    config.set('basic', 'hot_pix_x', '698,669')
    config.set('basic', 'hot_pix_y', '420,410')

    # dark sub parameters
    config.set('basic', 'aper_lim', '510')
    config.set('basic', 'cent_x', '512')
    config.set('basic', 'cent_y', '512')
    config.set('basic', 'a_lim ', '5')
    config.set('basic', 'k_lim', '0.02')
    config.set('basic', 'del_a', '0.1')
    config.set('basic', 'del_k', '0.001')
    config.set('basic', 'exp_dark', '1000')
    config.set('basic', 'nx', '1024')
    config.set('basic', 'ny', '1024')
    config.set('basic', 'drizzle_scale', '0.04')

    filt_list = sorted(np.unique(dict_files['filtername']))
    UV_filters = [125, 140, 150, 160]
    HA_filters = [775, 782]
    print ('\n...number of filters', len(filt_list))
    print ('\n...', filt_list)
    filt_sec = [str(x).lower() for x in filt_list]
    # filters =

    for i in range(len(np.unique(galaxies))):
        section = 'NULIRG%s' % (i + 1)
        config.add_section(section)
        params, params_gal = basic_params(config_dir + 'ULIRG_params.cfg', 'basic', section)
        # for j in range(len(params_gal.keys())):
        # config.set(section, params_gal[section].keys()[j], params_gal[section].value()[j])
        config.set(section, 'name', sorted(np.unique(galaxies))[i])
        config.set(section, 'bad_frames', params_gal['bad_frames'])  # user defined
        config.set(section, 'dark_frames', params_gal['dark_frames'])

        # number of filters will decide number of parameters here
        for k in range(4):

            config.set(section, filt_sec[k].replace('lp', ''), params_gal[filt_sec[k].replace('lp', '')])

            # sky subtraction
        config.set(section, 'sky_min', params_gal['sky_min'])  # user defined
        config.set(section, 'sky_max', params_gal['sky_max'])  # user defined

        # dark parameters
        config.set(section, 'aper_lim', params_gal['aper_lim'])
        config.set(section, 'cent_x', params_gal['cent_x'])
        config.set(section, 'cent_y', params_gal['cent_y'])
        config.set(section, 'dark_radii', params_gal['dark_radii'])  # user defined

        # xregister parameters UV

        config.set(section, 'xref_ref', params_gal['xreg_ref'])  # user defined
        config.set(section, 'xreg_lim', params_gal['xreg_lim'])  # user defined
        config.set(section, 'xreg_xwindow', params_gal['xreg_xwindow'])  # user defined
        config.set(section, 'xreg_ywindow', params_gal['xreg_ywindow'])  # user defined

        # xregister parameters HA

        config.set(section, 'xwin_ha', params_gal['xwin_ha'])  # user defined
        config.set(section, 'ywin_ha', params_gal['ywin_ha'])  # user defined

    print (config_file)
    with open(config_file, 'w') as configfile:
        config.write(configfile)
