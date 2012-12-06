"""
        Does miscellaneous jobs at the end, which assumes all other tasks are run.
"""

import numpy as N
import os
from image import *
import mylogger, os
try:
    import matplotlib.pyplot as pl
    import matplotlib.cm as cm
    has_pl = True
except ImportError:
    has_pl = False
import functions as func

class Op_cleanup(Op):
    """  """
    def __call__(self, img):

        mylog = mylogger.logging.getLogger("PyBDSM.Cleanup")

        ### plotresults for all gaussians together
        if img.opts.plot_allgaus and has_pl:
            pl.figure()
            pl.title('All gaussians including wavelet images')
            allgaus = img.gaussians
            if hasattr(img, 'atrous_gaussians'):
              for gg in img.atrous_gaussians:
                allgaus += gg

            for g in allgaus:
              ellx, elly = func.drawellipse(g)
              pl.plot(ellx, elly, 'r')

            from math import log10
            bdir = img.basedir + '/misc/'
            if not os.path.isdir(bdir): os.mkdir(bdir)
            im_mean = img.clipped_mean
            im_rms = img.clipped_rms
            low = 1.1*abs(img.min_value)
            low1 = 1.1*abs(N.min(im_mean-im_rms*5.0))
            if low1 > low: low = low1
            vmin = log10(im_mean-im_rms*5.0 + low)
            vmax = log10(im_mean+im_rms*15.0 + low)
            im = N.log10(img.ch0 + low)

            pl.imshow(N.transpose(im), origin='lower', interpolation='nearest',vmin=vmin, vmax=vmax, \
                      cmap=cm.gray); pl.colorbar()
            pl.savefig(bdir+'allgaussians.png')
            pl.close()

        img.completed_Ops.append('cleanup')



