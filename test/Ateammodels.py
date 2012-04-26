
import pylab as pl
import bdsm, pyfits
import numpy as N
import os, subprocess

from bdsm.FITS import Op_loadFITS
from bdsm.collapse import Op_collapse
from bdsm.preprocess import Op_preprocess
from bdsm.rmsimage import Op_rmsimage
from bdsm.threshold import Op_threshold
from bdsm.islands import Op_islands
import bdsm.functions as func
from bdsm.analysis import plotresults

""" Try blindly running bdsm to see if boxsize is ok, so fitting doesnt hang. Then try various segmenting algorithms which dont
depend on rms ? """

dir = "A-team/"
ls = subprocess.Popen(["ls",dir], stdout=subprocess.PIPE).communicate()[0]
ls = ls.split('\n')

files = []; rmsbox = []
chain = [Op_loadFITS(), Op_collapse(), Op_preprocess(), Op_rmsimage(), Op_threshold(), Op_islands()]

#ls = ['subim.fits']
bms = [(0.0015, 0.0015, 0.0)]
dir=''
for ifile, file in enumerate(ls):
  op = subprocess.Popen(["file",dir+file], stdout=subprocess.PIPE).communicate()[0]
  if "FITS image data" in op:
    print 'Processing ', file
    img = bdsm.execute(chain, {'fits_name': file, 'thresh':"hard", 'solnname' : 'new', 'beam' : bms[ifile]}), 'indir' : dir})
    files.append(file)
    rmsbox.append(img.opts.rms_box)

    thr = img.clipped_rms
    op1, markers1 = func.watershed(img.image, thr=thr*3.)


    pl.figure()
    pl.suptitle(img.filename)
    pl.subplot(2,2,1); pl.imshow(N.transpose(img.image), origin='lower', interpolation='nearest', vmin=-7*thr, vmax=15*thr); pl.title('Image')
    pl.subplot(2,2,2); pl.imshow(N.transpose(op1), origin='lower', interpolation='nearest'), pl.title('watershed1')
    pl.subplot(2,2,3); pl.imshow(N.transpose(markers1), origin='lower', interpolation='nearest'), pl.title('markers1')
    pl.subplot(2,2,4); plotresults(img, newfig=False, cbar=False)
    pl.savefig(dir+file+'_watershed.png')

  else:
    print dir+file+' is not a FITS file !!'





