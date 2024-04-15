"""Module shapelets

This will do all the shapelet analysis of islands in an image
"""
from __future__ import absolute_import

from .image import *
from .islands import *
from .shapelets import *
from . import mylogger
from . import statusbar
from . import multi_proc as mp
import itertools
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass
from . import functions as func
from .gausfit import find_bbox


class Op_shapelets(Op):
    """ Get the image and mask from each island and send it to
    shapelet programs which can then also be called seperately """

    def __call__(self, img):

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Shapefit")
        bar = statusbar.StatusBar('Decomposing islands into shapelets ...... : ', 0, img.nisl)
        opts = img.opts
        if img.opts.shapelet_do:
            if img.nisl == 0:
                mylog.warning("No islands found. Skipping shapelet decomposition.")
                img.completed_Ops.append('shapelets')
                return

            if not opts.quiet:
                bar.start()

            # Set up multiproccessing. First create a simple copy of the Image
            # object that contains the minimal data needed.
            opts_dict = opts.to_dict()
            img_simple = Image(opts_dict)
            img_simple.pixel_beamarea = img.pixel_beamarea
            img_simple.pixel_beam = img.pixel_beam
            img_simple.thresh_pix = img.thresh_pix
            img_simple.minpix_isl = img.minpix_isl
            img_simple.clipped_mean = img.clipped_mean
            img_simple.shape = img.ch0_arr.shape

            # Now call the parallel mapping function. Returns a list of
            # [beta, centre, nmax, basis, cf] for each island
            shap_list = mp.parallel_map(func.eval_func_tuple,
                        zip(itertools.repeat(self.process_island),
                        img.islands, itertools.repeat(img_simple),
                        itertools.repeat(opts)), numcores=opts.ncores,
                        bar=bar)

            for id, isl in enumerate(img.islands):
                beta, centre, nmax, basis, cf = shap_list[id]
                isl.shapelet_beta=beta
                isl.shapelet_centre=centre
                isl.shapelet_posn_sky=img.pix2sky(centre)
                isl.shapelet_posn_skyE=[0.0, 0.0, 0.0]
                isl.shapelet_nmax=nmax
                isl.shapelet_basis=basis
                isl.shapelet_cf=cf

            img.completed_Ops.append('shapelets')


    def process_island(self, isl, img, opts=None):
        """Processes a single island.

        Returns shapelet parameters.
        """
        if opts is None:
            opts = img.opts
        if opts.shapelet_gresid:
            shape = img.shape
            thresh= opts.fittedimage_clip
            model_gaus = N.zeros(shape, dtype=N.float32)
            for g in isl.gaul:
                C1, C2 = g.centre_pix
                b = find_bbox(thresh*isl.rms, g)
                bbox = N.s_[max(0, int(C1-b)):min(shape[0], int(C1+b+1)),
                            max(0, int(C2-b)):min(shape[1], int(C2+b+1))]
                x_ax, y_ax = N.mgrid[bbox]
                ffimg = func.gaussian_fcn(g, x_ax, y_ax)
                model_gaus[bbox] = model_gaus[bbox] + ffimg
            arr = isl.image - isl.islmean - model_gaus[tuple(isl.bbox)]
        else:
            arr = isl.image - isl.islmean
        mask = isl.mask_active
        basis = opts.shapelet_basis
        beam_pix = img.pixel_beam()
        mode = opts.shapelet_fitmode
        if mode != 'fit':
            mode = ''
        fixed = (0,0,0)
        (beta, centre, nmax) = self.get_shapelet_params(arr, mask, basis, beam_pix, fixed, N.array(isl.origin), mode)

        cf = decompose_shapelets(arr, mask, basis, beta, centre, nmax, mode)

        return [beta, tuple(N.array(centre) + N.array(isl.origin)), nmax, basis, cf]


    def get_shapelet_params(self, image, mask, basis, beam_pix, fixed, ori, mode, beta=None, cen=None, nmax=None):
        """ This takes as input an image, its mask (false=valid), basis="cartesian"/"polar",
        fixed=(i,j,k) where i,j,k =0/1 to calculate or take as fixed for (beta, centre, nmax),
        beam_pix has the beam in (pix_fwhm, pix_fwhm, deg),
        beta (the scale), cen (centre of basis expansion), nmax (max order). The output
        is an updated set of values of (beta, centre, nmax). If fixed is 1 and the value is not
        specified as an argument, then fixed is taken as 0."""
        from math import sqrt, log, floor
        from . import functions as func
        import numpy as N

        if fixed[0]==1 and beta is None: fixed[0]=0
        if fixed[1]==1 and cen is None: fixed[1]=0
        if fixed[2]==1 and nmax is None: fixed[2]=0

        if fixed[0]*fixed[1]==0:
            (m1, m2, m3)=func.moment(image, mask)

        if fixed[0]==0:
            try:
                beta = sqrt(m3[0]*m3[1])*2.0
            except ValueError:
                beta = 0.5
            if beta == 0.0:
                beta = 0.5
        if fixed[1]==0:
            cen=m2
        if fixed[2]==0:
            (n, m)=image.shape
            nmax=int(round(sqrt(1.0*n*n+m*m)/beam_pix[1]))-1
            nmax=min(max(nmax*2+2,10),10)                      # totally ad hoc
            npix = N.product(image.shape)-N.sum(mask)
            if nmax*nmax >= n*m : nmax = int(floor(sqrt(npix-1)))  # -1 is for when n*m is a perfect square
            if mode == 'fit':   # make sure npara <= npix
                nmax_max = int(round(0.5*(-3+sqrt(1+8*npix))))
                nmax=min(nmax, nmax_max)

            betarange=[0.5,sqrt(beta*max(n,m))]  # min, max

            if fixed[1]==0:
                cen=shape_findcen(image, mask, basis, beta, nmax, beam_pix) # + check_cen_shapelet
                #print 'First Centre = ',cen,N.array(cen)+ori

                from time import time
                t1 = time()
                if fixed[0]==0:
                    beta, err=shape_varybeta(image, mask, basis, beta, cen, nmax, betarange, plot=False)
                t2 = time()
                #print 'TIME ',t2-t1, '\n'
                #print 'Final Beta = ',beta, err

        if fixed[1]==0 and fixed[0]==0:
            cen=shape_findcen(image, mask, basis, beta, nmax, beam_pix) # + check_cen_shapelet
            #print 'Final Cen = ',N.array(cen)+ori

        return beta, cen, nmax
