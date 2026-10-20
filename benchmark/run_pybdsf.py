#!/usr/bin/env python
# -*- coding: utf-8 -*-

import timeit, bdsf, sys

# 127 MB
fitsname = 'NGC6048_9-MFS-image.fits'

# 1.5 GB
#fitsname = 'image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'

start = timeit.default_timer()
img = bdsf.process_image(fitsname,
                         thresh_isl = 4.0,
                         thresh_pix = 5.0,
                         rms_box = (150, 15),
                         rms_map = True,
                         mean_map = 'zero',
                         ini_method = 'intensity',
                         adaptive_rms_box = True,
                         adaptive_thresh = 150,
                         rms_box_bright = (60, 15),
                         group_by_isl = False,
                         group_tol = 10.0,
                         output_opts = True,
                         output_all = True,
                         atrous_do = True,
                         atrous_jmax = 4,
                         flagging_opts = True,
                         flag_maxsize_fwhm = 0.5,
                         advanced_opts = True,
                         blank_limit = None,
                         frequency = 143.65e6,
                         clobber = True)

stop = timeit.default_timer()
total_timeHi = stop - start
mins, secs  = divmod(total_timeHi, 60)
hours, mins = divmod(mins, 60)
sys.stdout.write("PyBDSF running time: %d:%d:%d\n" % (hours, mins, secs))