
import matplotlib.cm as cm
import scipy.ndimage as nd
from bdsm.const import fwsig
from bdsm.gausfit import Op_gausfit as gg
import bdsm.functions as func
from _cbdsm import MGFunction
from _cbdsm import lmder_fit, dn2g_fit, dnsg_fit
import numpy as N
from copy import deepcopy as cp

for isl in img.islands:
  #isl = img.islands[153]
  if isl.ngaus > 1:
    thr = isl.mean + img.opts.thresh_pix*isl.rms
    im = isl.image; mask = isl.mask_active; av = img.clipped_mean; im1 = cp(im)
    ind = N.array(N.where(~mask)).transpose()
    ind = [tuple(coord) for coord in ind if im[tuple(coord)] > thr]
    n, m = isl.shape; iniposn = []; inipeak = []
    for c in ind:
      goodlist = [im[i,j] for i in range(c[0]-1,c[0]+2) for j in range(c[1]-1,c[1]+2) \
                   if i>=0 and i<n and j>=0 and j<m and (i,j) != c]
      peak = N.sum(im[c] > goodlist) == len(goodlist)
      if peak:
        iniposn.append(c); inipeak.append(im[c])
    nmulsrc = len(iniposn) 
    if nmulsrc > 1:
      markers = N.zeros(im.shape, int)
      markers[0,0] = 1
      for ipk in range(nmulsrc):
        pk = inipeak[ipk]; x, y = iniposn[ipk]
        markers[int(round(x)), int(round(y))] = ipk+2
      im2 = N.zeros(im.shape, int)
      im1 = im1 - im1.min()
      im1 = im1/im1.max()*255
      im1 = 255-im1
      nd.watershed_ift(N.array(im1, N.uint8), markers, output = im2)
      fcn = MGFunction(im, isl.mask_active, 1)
      fit = lmder_fit
      gg1 = gg()
      for ipk in range(nmulsrc):
        ind = ipk+2
        mom = func.momanalmask_gaus(im, im2, ind, 1.0, True)
        indd = N.where(im2==ind)
        mom[3] = 3.0; mom[4]=3.0
        g = [float(N.max(im[indd])), int(round(mom[1])), int(round(mom[2])), mom[3]/fwsig, mom[4]/fwsig, mom[5]]
        gg1.add_gaussian(fcn, g, dof = isl.size_active)
        print g
      fit(fcn, final=0, verbose=True)
      print fcn.parameters
      import pylab as pl
      pl.figure()
      pl.subplot(2,2,1);pl.imshow(N.transpose(im), interpolation='nearest', origin='lower'); pl.title(str(isl.island_id))
      pl.subplot(2,2,2);pl.imshow(N.transpose(im1), interpolation='nearest', origin='lower'); pl.title(str(isl.island_id))
      pl.subplot(2,2,3);pl.imshow(N.transpose(im2), interpolation='nearest', origin='lower'); pl.title(str(isl.island_id))
      for g in fcn.parameters:
            A, x1, x2, s1, s2, th = g
            s1, s2 = map(abs, [s1, s2])
            if s1 < s2:   # s1 etc are sigma
              ss1=s2; ss2=s1; th1 = divmod(th+90.0, 180)[1]
            else:
              ss1=s1; ss2=s2; th1 = divmod(th, 180)[1]
            c = [A, x1, x2, ss1, ss2, th1]
            x, y = N.indices(isl.shape)
            x2, y2 = func.drawellipse(c)
            #x2 = x2 + isl.origin[0]; y2 = y2 + isl.origin[1]
            pl.subplot(2,2,4); pl.plot(x2, y2, '-r')
            pl.imshow(N.transpose(im), origin='lower', interpolation='nearest')



import matplotlib.cm as cm
import scipy.ndimage as nd
import numpy as N
from bdsm.const import fwsig
from bdsm.gausfit import Op_gausfit as gg
import bdsm.functions as func
from _cbdsm import MGFunction
from _cbdsm import lmder_fit, dn2g_fit, dnsg_fit
image = N.zeros((100,100))
markers = N.zeros(image.shape, int)
op1 = N.zeros(image.shape, int)
op2 = N.zeros(image.shape, int)
x, y = N.indices(image.shape)
peaks = [2.0, 8.0, 8.0, 2.0]
posns = [(30, 20), (50, 20), (30, 70), (50, 70)]
bmaj = [2.0, 12.0, 2.0, 12.0]
brat = [2.0, 2.0, 2.0, 2.0]
markers[10,10] = 1
for ig in range(len(peaks)):
  g = peaks[ig]*N.exp(-0.5*((x-posns[ig][0])*(x-posns[ig][0])+(y-posns[ig][1])*(y-posns[ig][1]))  \
                  /(bmaj[ig]*bmaj[ig]/brat[ig]))
  image = image + g
  markers[int(round(posns[ig][0])), int(round(posns[ig][1]))] = ig+2

image1 = image - image.min()
image1 = image1/image1.max()*255
image1 = 255-image1
nd.watershed_ift(N.array(image1, N.uint8), markers, output = op1)
pl.figure();pl.imshow(N.transpose(image), interpolation='nearest', origin='lower'); pl.title('orig'); pl.colorbar()
pl.figure();pl.imshow(N.transpose(image1), interpolation='nearest', origin='lower'); pl.title('input1'); pl.colorbar()
pl.figure();pl.imshow(N.transpose(op1), interpolation='nearest', origin='lower'); pl.title('output1'); pl.colorbar()
pl.figure();pl.imshow(N.transpose(markers), interpolation='nearest', origin='lower'); pl.title('markers'); pl.colorbar()



