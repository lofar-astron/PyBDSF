
import pylab as pl
import bdsf, pyfits
import numpy as N
import os, subprocess

from bdsf.FITS import Op_loadFITS
from bdsf.collapse import Op_collapse
from bdsf.preprocess import Op_preprocess
from bdsf.rmsimage import Op_rmsimage
from bdsf.threshold import Op_threshold
from bdsf.islands import Op_islands
import bdsf.functions as func
from bdsf.analysis import plotresults

""" Try blindly running bdsf to see if boxsize is ok, so fitting doesnt hang. Then try various segmenting algorithms which dont
depend on rms ? """

directory = "A-team/"
ls = subprocess.Popen(["ls",directory], stdout=subprocess.PIPE).communicate()[0]
ls = ls.split('\n')

files = []; rmsbox = []
chain = [Op_loadFITS(), Op_collapse(), Op_preprocess(), Op_rmsimage(), Op_threshold(), Op_islands()]

#ls = ['subim.fits']
bms = [(0.0015, 0.0015, 0.0)]
directory=''
for ifile, filename in enumerate(ls):
  op = subprocess.Popen(["file",directory+filename], stdout=subprocess.PIPE).communicate()[0]
  if "FITS image data" in op:
    print 'Processing ', filename
    img = bdsf.execute(chain, {'fits_name': filename, 'thresh':"hard", 'solnname' : 'new', 'beam' : bms[ifile]), 'indir' : directory})
    files.append(filename)
    rmsbox.append(img.opts.rms_box)

    thr = img.clipped_rms
    op1, markers1 = func.watershed(img.image, thr=thr*3.)


    pl.figure()
    pl.suptitle(img.filename)
    pl.subplot(2,2,1); pl.imshow(N.transpose(img.image), origin='lower', interpolation='nearest', vmin=-7*thr, vmax=15*thr); pl.title('Image')
    pl.subplot(2,2,2); pl.imshow(N.transpose(op1), origin='lower', interpolation='nearest'), pl.title('watershed1')
    pl.subplot(2,2,3); pl.imshow(N.transpose(markers1), origin='lower', interpolation='nearest'), pl.title('markers1')
    pl.subplot(2,2,4); plotresults(img, newfig=False, cbar=False)
    pl.savefig(directory+filename+'_watershed.png')

  else:
    print directory+filename+' is not a FITS file !!'





