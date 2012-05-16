
""" Has all the miscellaneous complicated modules which are not general enough for functions.py """

def isl_tosplit(isl, img):
    """ Splits an island and sends back parameters """
    import functions as func
    import numpy as N

    size_extra5 = img.opts.splitisl_size_extra5
    frac_bigisl3 = img.opts.splitisl_frac_bigisl3

    connected, count = func.connect(isl.mask_active)
    index = 0
    n_subisl3, labels3, isl_pixs3 = func.open_isl(isl.mask_active, 3)
    n_subisl5, labels5, isl_pixs5 = func.open_isl(isl.mask_active, 5)
    isl_pixs3, isl_pixs5 = N.array(isl_pixs3), N.array(isl_pixs5)
    
                                # take open 3 or 5 
    open3, open5 = False, False
    if n_subisl3 > 0 and isl_pixs3 != None:                                 # open 3 breaks up island
      max_sub3 = N.max(isl_pixs3)
      if max_sub3 < frac_bigisl3 : open3 = True       # if biggest sub island isnt too big
    if n_subisl5 > 0 and isl_pixs5 != None:                                 # open 5 breaks up island
      max_sub5 = N.max(isl_pixs5)                     # if biggest subisl isnt too big OR smallest extra islands add upto 10 %
      if (max_sub5 < 0.75*max_sub3) or (N.sum(N.sort(isl_pixs5)[:len(isl_pixs5)-n_subisl3]) > size_extra5):
        open5 = True
                                # index=0 => dont split
    if open5: index = 5; n_subisl = n_subisl5; labels = labels5
    else:
      if open3: index = 3; n_subisl = n_subisl3; labels = labels3
      else: index = 0
    convex_def =  func.convexhull_deficiency(isl) 
    #print 'CONVEX = ',convex_def

    if img.opts.plot_islands:
        import matplotlib.pyplot as pl
        #import pylab as pl
        pl.figure()
        pl.suptitle('Island '+str(isl.island_id) + ' ' + repr(img.waveletimage))
        pl.subplot(2,2,1); pl.imshow(N.transpose(isl.image*~isl.mask_active), origin='lower', interpolation='nearest'); pl.title('Image')
        pl.subplot(2,2,2); pl.imshow(N.transpose(labels3), origin='lower', interpolation='nearest'); pl.title('labels3')
        pl.subplot(2,2,3); pl.imshow(N.transpose(labels5), origin='lower', interpolation='nearest'); pl.title('labels5')

    if index == 0: return [index, n_subisl5, labels5]
    else: return [index, n_subisl, labels]

