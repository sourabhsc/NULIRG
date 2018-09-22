import imports_all
from astropy.table import Table, Column, MaskedColumn

filt = ['775', '782']


def main(config_file):
    tab = Table.read(flt_files, format='csv')

    for i in range(5):

        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)

        primary_dir = params['data_dir'] + gal_name + '/'

        for j in range(2):
            if tab['galaxy'] == gal_name and tab['filtername'] == 'F775W':
                angle = float(tab['orientation'])
            print (angle)
            psf = params['script_dir'] + 'PSF_%s_gal%s.fits' % (filt[j], i + 1)
            psf_rot = psf.replace('gal%s', 'gal%s_rotate')
            iraf.roate(psf, psf_rot, -angle, interpolant='nearest')


if __name__ == '__main__':
    main(config_file)
