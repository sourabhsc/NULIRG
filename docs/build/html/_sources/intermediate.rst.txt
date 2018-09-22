*****
INTER
*****



Intermediate files
==================


==== ================ ==================== ====================================== ======================= ================================= ==============
col0          process           subprocess                                   name                location                       parent code     sub_module
==== ================ ==================== ====================================== ======================= ================================= ==============
   0  sky subtraction                   --                   prefix+_sky_flt.fits /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py   sky_dark_sub
   1               --              sky sub                      FLT_sky_gal_#.txt  /DIR/INTERMEDIATE_TXT/ /SOURCE/galaxy_sky_dark_sub_v4.py   sky_dark_sub
   2 dark subtraction     sextract catalog                    prefix+_catalog.cat  /DIR/INTERMEDIATE_TXT/ /SOURCE/galaxy_sky_dark_sub_v4.py sextractor_seg
   3               --              seg_map                   prefix+_seg_map.fits /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py sextractor_seg
   4               --           masked_flt                      prefix+masked_flt /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py    galaxy_mask
   5               --             flt_mask                        prefix+flt_mask /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py    galaxy_mask
   6               --               smooth                      prefix+flt_smooth /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py     smooth_gal
   7               --        selected dark                    prefix+dark_**.fits /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py       dark_sub
   8               -- dark minimizer value               prefix+dark_**_diff.fits /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py       dark_sub
   9               --       minimized dark                    prefix+dark_**.fits /DIR/INTERMEDIATE_FITS/ /SOURCE/galaxy_sky_dark_sub_v4.py       dark_sub
  10               --        dark sub data                   FLT_drk_v3_gal_#.txt  /DIR/INTERMEDIATE_TXT/ /SOURCE/galaxy_sky_dark_sub_v4.py       dark_sub
  11               --       dark minimizer                   prefix+minimizer.png  /DIR/INTERMEDIATE_PNG/ /SOURCE/galaxy_sky_dark_sub_v4.py       dark_sub
  12               --    dark verification                prefix+verification.png  /DIR/INTERMEDIATE_PNG/ /SOURCE/galaxy_sky_dark_sub_v4.py       dark_sub
  13         aligning           rotate SCI     prefix+_FILTER_SUB_rotate_flt.fits /DIR/INTERMEDIATE_FITS/     /SOURCE/galaxy_xreg_rotate.py         rotate
  14               --            rotate DQ  prefix+_FILTER_SUB_rotate_flt_DQ.fits /DIR/INTERMEDIATE_FITS/     /SOURCE/galaxy_xreg_rotate.py         rotate
  15               --           rotate ERR prefix+_FILTER_SUB_rotate_flt_err.fits /DIR/INTERMEDIATE_FITS/     /SOURCE/galaxy_xreg_rotate.py         rotate
  16               --            shift SCI           prefix+_FILTER_SUB_shift.txt  /DIR/INTERMEDIATE_TXT/     /SOURCE/galaxy_xreg_rotate.py        xregflt
  17               --            shift ERR       prefix+_FILTER_SUB_shift_err.txt  /DIR/INTERMEDIATE_TXT/     /SOURCE/galaxy_xreg_rotate.py        xregflt
  18               --          xregflt SCI        prefix+_FILTER_SUB_xregflt.fits /DIR/INTERMEDIATE_FITS/     /SOURCE/galaxy_xreg_rotate.py        xregflt
  19               --          xregflt ERR    prefix+_FILTER_SUB_xregflt_err.fits /DIR/INTERMEDIATE_FITS/     /SOURCE/galaxy_xreg_rotate.py        xregflt
  20    header_rename           shifted DQ               prefix+_shifted_err.fits /DIR/INTERMEDIATE_FITS/   /SOURCE/galaxy_header_rename.py        shifted
  21               --          shifted ERR                prefix+_shifted_DQ.fits /DIR/INTERMEDIATE_FITS/   /SOURCE/galaxy_header_rename.py        shifted
  22               --               allext                    prefix+_allext.fits                   /DIR/   /SOURCE/galaxy_header_rename.py        shifted
  23          drizzle                   --                gal#_UV_FILTER_drz.fits                   /DIR/    /SOURCE/galaxy_astrodrizzle.py   astrodrizzle
  24               --        drz intermed.                         multiple files     /DIR/DRIZZLE_INTER/    /SOURCE/galaxy_astrodrizzle.py   astrodrizzle
  25         psfmatch             psfmatch           gal#_UV_FILTER_psfmatch.fits                   /DIR/        /SOURCE/galaxy_psfmatch.py            psf
  26               HA                   HA                           gal#_HA.fits                   /DIR/       /SOURCE/galaxy_ha_cut_v3.py             ha
  27               --                  cut gal#_HA_FILTERha_scale_04_cut_drc.fits               /DIR/DRC/       /SOURCE/galaxy_ha_cut_v3.py         ha_cut
  28               --                align    gal#_HA_FILTERha_type_UV_align.fits /DIR/INTERMEDIATE_FITS/       /SOURCE/galaxy_ha_cut_v3.py        ha_xreg
  29               --                  sky gal#_HA_FILTERha_scale_04_sky_drc.fits               /DIR/DRC/       /SOURCE/galaxy_ha_cut_v3.py         ha_sky
==== ================ ==================== ====================================== ======================= ================================= ==============
