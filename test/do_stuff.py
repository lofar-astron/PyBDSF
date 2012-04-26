


"""make watershed images for each island in isls """

def do_ws(isls, crms):
    import bdsm.functions as func
    import os, subprocess
    import pylab as pl
    import numpy as N

    thr = crms
    for isl in isls:
    
      image = isl.image*~isl.mask_active
      op1, markers1 = func.watershed(image, thr=thr*3.)
    
      pl.figure()
      pl.suptitle('Island '+str(isl.island_id))
      pl.subplot(2,2,1); pl.imshow(N.transpose(image), origin='lower', interpolation='nearest', vmin=-7*thr, vmax=15*thr); pl.title('Image')
      pl.subplot(2,2,2); pl.imshow(N.transpose(op1*~isl.mask_active), origin='lower', interpolation='nearest'); pl.title('watershed1')
      pl.subplot(2,2,3); pl.imshow(N.transpose(markers1*~isl.mask_active), origin='lower', interpolation='nearest'); pl.title('markers1')

def open_isl(isls, crms):
    import pylab as pl
    import scipy.ndimage as nd
    import numpy as N

    thr = crms
    ft1 = N.array(((1,0,1), (0,1,0), (1,0,1)), int)
    ft2 = N.array(((0,1,0), (1,1,1), (0,1,0)), int)
    ft3 = N.ones((3,3), int)
    ft5 = N.ones((5,5), int)
    for isl in isls:
      ma = ~isl.mask_active
      open1 = nd.binary_opening(ma, ft1)
      open2 = nd.binary_opening(ma, ft2)
      open3 = nd.binary_opening(ma, ft3)
      open5 = nd.binary_opening(ma, ft5)

      pl.figure()
      pl.suptitle('Island '+str(isl.island_id))
      pl.subplot(2,2,1); pl.imshow(N.transpose(isl.image), origin='lower', interpolation='nearest'); pl.title('Image')
      pl.subplot(2,2,2); pl.imshow(N.transpose(ma), origin='lower', interpolation='nearest'); pl.title('mask')
      pl.subplot(2,2,3); pl.imshow(N.transpose(open3), origin='lower', interpolation='nearest'); pl.title('open 3x3')
      pl.subplot(2,2,4); pl.imshow(N.transpose(open5), origin='lower', interpolation='nearest'); pl.title('open 5x5')
      #pl.subplot(2,2,3); pl.imshow(N.transpose(open1), origin='lower', interpolation='nearest'); pl.title('open diag')
      #pl.subplot(2,2,4); pl.imshow(N.transpose(open2), origin='lower', interpolation='nearest'); pl.title('open str')
      pl.savefig('cyga_p_w12_bigisl_'+str(isl.island_id)+'_open.png')

    


