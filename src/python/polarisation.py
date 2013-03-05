"""Module polarisation.

This module finds the Q, U, and V fluxes, the total, linear, and circular
polarisation fractions and the linear polarisation angle of each source identified
by gaul2srl. The position angle is defined from North, with positive angles
towards East.

"""

from image import *
from islands import *
from gausfit import Gaussian
from gaul2srl import *
from preprocess import Op_preprocess
from rmsimage import Op_rmsimage
from threshold import Op_threshold
from islands import Op_islands
from gausfit import Op_gausfit

from gaul2srl import Op_gaul2srl
from make_residimage import Op_make_residimage
from const import fwsig
import mylogger
import numpy as N
import functions as func
import statusbar

### Insert polarization attributes into Gaussian and Source classes
Gaussian.total_flux_Q        = Float(doc="Total flux density (Jy), Stokes Q", colname='Total_Q',
                                   units='Jy')
Gaussian.total_fluxE_Q       = Float(doc="Error in total flux density (Jy), Stokes Q", colname='E_Total_Q',
                                   units='Jy')
Gaussian.total_flux_U        = Float(doc="Total flux density (Jy), Stokes U", colname='Total_U',
                                   units='Jy')
Gaussian.total_fluxE_U       = Float(doc="Error in total flux density (Jy), Stokes U", colname='E_Total_U',
                                   units='Jy')
Gaussian.total_flux_V        = Float(doc="Total flux density (Jy), Stokes V", colname='Total_V',
                                   units='Jy')
Gaussian.total_fluxE_V       = Float(doc="Error in total flux density (Jy), Stokes V", colname='E_Total_V',
                                   units='Jy')
Gaussian.lpol_fraction       = Float(doc="Linear polarisation fraction",
                                   colname='Linear_Pol_frac', units=None)
Gaussian.lpol_fraction_loerr   = Float(doc="Linear polarisation fraction low error",
                                   colname='Elow_Linear_Pol_frac', units=None)
Gaussian.lpol_fraction_hierr   = Float(doc="Linear polarisation fraction high error",
                                   colname='Ehigh_Linear_Pol_frac', units=None)
Gaussian.cpol_fraction       = Float(doc="Circular polarisation fraction",
                                   colname='Circ_Pol_Frac', units=None)
Gaussian.cpol_fraction_loerr   = Float(doc="Circular polarisation fraction low error",
                                   colname='Elow_Circ_Pol_Frac', units=None)
Gaussian.cpol_fraction_hierr   = Float(doc="Circular polarisation fraction high error",
                                   colname='Ehigh_Circ_Pol_Frac', units=None)
Gaussian.tpol_fraction       = Float(doc="Total polarisation fraction",
                                   colname='Total_Pol_Frac', units=None)
Gaussian.tpol_fraction_loerr   = Float(doc="Total polarisation fraction low error",
                                   colname='Elow_Total_Pol_Frac', units=None)
Gaussian.tpol_fraction_hierr   = Float(doc="Total polarisation fraction high error",
                                   colname='Ehigh_Total_Pol_Frac', units=None)
Gaussian.lpol_angle          = Float(doc="Polarisation angle (deg from North towards East)",
                                   colname='Linear_Pol_Ang', units='deg')
Gaussian.lpol_angle_err      = Float(doc="Polarisation angle error (deg)",
                                   colname='E_Linear_Pol_Ang', units='deg')

Source.total_flux_Q        = Float(doc="Total flux density (Jy), Stokes Q", colname='Total_Q',
                                   units='Jy')
Source.total_fluxE_Q       = Float(doc="Error in total flux density (Jy), Stokes Q", colname='E_Total_Q',
                                   units='Jy')
Source.total_flux_U        = Float(doc="Total flux density (Jy), Stokes U", colname='Total_U',
                                   units='Jy')
Source.total_fluxE_U       = Float(doc="Error in total flux density (Jy), Stokes U", colname='E_Total_U',
                                   units='Jy')
Source.total_flux_V        = Float(doc="Total flux density (Jy), Stokes V", colname='Total_V',
                                   units='Jy')
Source.total_fluxE_V       = Float(doc="Error in total flux density (Jy), Stokes V", colname='E_Total_V',
                                   units='Jy')
Source.lpol_fraction       = Float(doc="Linear polarisation fraction",
                                   colname='Linear_Pol_frac', units=None)
Source.lpol_fraction_loerr   = Float(doc="Linear polarisation fraction low error",
                                   colname='Elow_Linear_Pol_frac', units=None)
Source.lpol_fraction_hierr   = Float(doc="Linear polarisation fraction high error",
                                   colname='Ehigh_Linear_Pol_frac', units=None)
Source.cpol_fraction       = Float(doc="Circular polarisation fraction",
                                   colname='Circ_Pol_Frac', units=None)
Source.cpol_fraction_loerr   = Float(doc="Circular polarisation fraction low error",
                                   colname='Elow_Circ_Pol_Frac', units=None)
Source.cpol_fraction_hierr   = Float(doc="Circular polarisation fraction high error",
                                   colname='Ehigh_Circ_Pol_Frac', units=None)
Source.tpol_fraction       = Float(doc="Total polarisation fraction",
                                   colname='Total_Pol_Frac', units=None)
Source.tpol_fraction_loerr   = Float(doc="Total polarisation fraction low error",
                                   colname='Elow_Total_Pol_Frac', units=None)
Source.tpol_fraction_hierr   = Float(doc="Total polarisation fraction high error",
                                   colname='Ehigh_Total_Pol_Frac', units=None)
Source.lpol_angle          = Float(doc="Polarisation angle (deg from North towards East)",
                                   colname='Linear_Pol_Ang', units='deg')
Source.lpol_angle_err      = Float(doc="Polarisation angle error (deg)",
                                   colname='E_Linear_Pol_Ang', units='deg')

class Op_polarisation(Op):
    """ Finds the flux in each Stokes and calculates the polarisation fraction
    and angle.

    Fluxes are calculated by summing all nonmasked pixels assigned to
    the Gaussian. If a pixel contains contributions from two or more
    Gaussians, its flux is divided between the Gaussians by the ratio of
    fluxes that they contribute to the pixel. Errors on the fluxes are
    derived by summing the same pixels in the rms maps in quadrature.
    The results are stored in the Gaussian and Source structures.

    Fits are also done to the polarized intensity (PI) image to
    determine if there are any islands of emission that lie outside
    those found in the I image. If there are, they are fit and the
    process above is done for them too.

    For linearly polarised emission, the signal and noise add
    vectorially, giving a Rice distribution (Vinokur 1965) instead of a
    Gaussian one. To correct for this, a bias is estimated and removed
    from the polarisation fraction using the same method used for the
    NVSS catalog (see ftp://ftp.cv.nrao.edu/pub/nvss/catalog.ps). Errors
    on the linear and total polarisation fractions and polarisation
    angle are estimated using the debiased polarised flux and standard
    error propagation. See Sparks & Axon (1999) for a more detailed
    treatment.

    Prerequisites: module gaul2srl should be run first."""

    def __call__(self, img):
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Polarisatn")
        if img.opts.polarisation_do:
          mylog.info('Extracting polarisation properties for all sources')

          # Run gausfit and gual2srl on PI image to look for polarized sources
          # undetected in I
          fit_PI = img.opts.pi_fit
          n_new = 0
          ch0_pi = N.sqrt(img.ch0_Q**2 + img.ch0_U**2)
          img.ch0_pi = ch0_pi

          if fit_PI:
              from . import _run_op_list
              mylogger.userinfo(mylog, "\nChecking PI image for new sources")

              mask = img.mask
              minsize = img.opts.minpix_isl

              # Set up image object for PI image.
              pi_chain, pi_opts = self.setpara_bdsm(img)
              pimg = Image(pi_opts)
              pimg.beam = img.beam
              pimg.pixel_beam = img.pixel_beam
              pimg.pixel_beamarea = img.pixel_beamarea
              pimg.log = 'PI.'
              pimg.pix2beam = img.pix2beam
              pimg.beam2pix = img.beam2pix
              pimg.pix2sky = img.pix2sky
              pimg.sky2pix = img.sky2pix
              pimg.pix2coord = img.pix2coord
              pimg.wcs_obj = img.wcs_obj
              pimg.mask = mask
              pimg.ch0 = ch0_pi
              pimg._pi = True

              success = _run_op_list(pimg, pi_chain)
              if not success:
                  return

              img.pi_islands = pimg.islands
              img.pi_gaussians = pimg.gaussians
              img.pi_sources = pimg.sources

              # Now check for new sources in the PI image that are not
              # found in the Stokes I image. If any new sources are found,
              # adjust their IDs to follow after those found in I.
              new_isl = []
              new_src = []
              n_new_src = 0
              isl_id = img.islands[-1].island_id
              src_id = img.sources[-1].source_id
              gaus_id = img.gaussians[-1].gaus_num
              for pi_isl in pimg.islands:
                  new_sources = []
                  for pi_src in pi_isl.sources:
                      if img.pyrank[int(img.sky2pix(pi_src.posn_sky_max)[0]),
                                    int(img.sky2pix(pi_src.posn_sky_max)[1])] == -1:
                          src_id += 1
                          pi_src._pi = True
                          pi_src.island_id = isl_id
                          pi_src.source_id = src_id
                          pi_src.spec_indx = N.NaN
                          pi_src.e_spec_indx = N.NaN
                          pi_src.spec_norm = N.NaN
                          pi_src.specin_flux = [N.NaN]
                          pi_src.specin_fluxE = [N.NaN]
                          pi_src.specin_freq = [N.NaN]
                          pi_src.specin_freq0 = N.NaN
                          new_sources.append(pi_src)
                          new_src.append(pi_src)
                          n_new_src += 1
                          for g in pi_src.gaussians:
                              gaus_id += 1
                              g.gaus_num = gaus_id
                  if len(new_sources) > 0:
                      isl_id += 1
                      pi_isl.sources = new_sources
                      pi_isl.island_id = isl_id
                      pi_isl._pi = True
                      new_isl.append(pi_isl)

              n_new = len(new_isl)
              mylogger.userinfo(mylog, "New sources found in PI image", '%i (%i total)' %
                                (n_new, img.nsrc+n_new))

          if n_new > 0:
              img.islands += new_isl
              img.sources += new_src
              img.nsrc += n_new_src

          bar = statusbar.StatusBar('Calculating polarisation properties ....  : ', 0, img.nsrc)
          if img.opts.quiet == False:
              bar.start()

          for isl in img.islands:
            isl_bbox = isl.bbox
            ch0_I = img.ch0[isl_bbox]
            ch0_Q = img.ch0_Q[isl_bbox]
            ch0_U = img.ch0_U[isl_bbox]
            ch0_V = img.ch0_V[isl_bbox]
            ch0_images = [ch0_I, ch0_Q, ch0_U, ch0_V]

            for i, src in enumerate(isl.sources):
                # For each source, assume the morphology does not change
                # across the Stokes cube. This assumption allows us to fit
                # the Gaussians of each source to each Stokes image by
                # simply fitting only the overall normalizations of the
                # individual Gaussians.
                #
                # First, fit all source Gaussians to each Stokes image:
                x, y = N.mgrid[isl_bbox]
                gg = src.gaussians
                fitfix = N.ones(len(gg)) # fit only normalization
                srcmask = isl.mask_active
                total_flux = N.zeros((4, len(fitfix))) # array of fluxes: N_Stokes x N_Gaussians
                errors = N.zeros((4, len(fitfix))) # array of fluxes: N_Stokes x N_Gaussians

                for sind, image in enumerate(ch0_images):
                    if (sind==0 and hasattr(src, '_pi')) or sind > 0: # Fit I only for PI sources
                        p, ep = func.fit_mulgaus2d(image, gg, x, y, srcmask, fitfix)
                        for ig in range(len(fitfix)):
                            center_pix = (p[ig*6 + 1], p[ig*6 + 2])
                            bm_pix = N.array([img.pixel_beam(location=center_pix)[0], img.pixel_beam(location=center_pix)[1], img.pixel_beam(location=center_pix)[2]])
                            total_flux[sind, ig] = p[ig*6]*p[ig*6+3]*p[ig*6+4]/(bm_pix[0]*bm_pix[1])
                        p = N.insert(p, N.arange(len(fitfix))*6+6, total_flux[sind])
                        if sind > 0:
                            rms_img = img.rms_QUV[sind-1]
                        else:
                            rms_img = img.rms
                        if len(rms_img.shape) > 1:
                            rms_isl = rms_img[isl.bbox].mean()
                        else:
                            rms_isl = rms_img
                        errors[sind] = func.get_errors(img, p, rms_isl)[6]

                # Now, assign fluxes to each Gaussian.
                src_flux_I = 0.0
                src_flux_Q = 0.0
                src_flux_U = 0.0
                src_flux_V = 0.0
                src_flux_I_err_sq = 0.0
                src_flux_Q_err_sq = 0.0
                src_flux_U_err_sq = 0.0
                src_flux_V_err_sq = 0.0

                for ig, gaussian in enumerate(src.gaussians):
                    flux_I = total_flux[0, ig]
                    flux_I_err = abs(errors[0, ig])
                    flux_Q = total_flux[1, ig]
                    flux_Q_err = abs(errors[1, ig])
                    flux_U = total_flux[2, ig]
                    flux_U_err = abs(errors[2, ig])
                    flux_V = total_flux[3, ig]
                    flux_V_err = abs(errors[3, ig])

                    if hasattr(src, '_pi'):
                        gaussian.total_flux = flux_I
                        gaussian.total_fluxE = flux_I_err
                    gaussian.total_flux_Q = flux_Q
                    gaussian.total_flux_U = flux_U
                    gaussian.total_flux_V = flux_V
                    gaussian.total_fluxE_Q = flux_Q_err
                    gaussian.total_fluxE_U = flux_U_err
                    gaussian.total_fluxE_V = flux_V_err

                    if hasattr(src, '_pi'):
                        src_flux_I += flux_I
                        src_flux_I_err_sq += flux_I_err**2
                    src_flux_Q += flux_Q
                    src_flux_U += flux_U
                    src_flux_V += flux_V
                    src_flux_Q_err_sq += flux_Q_err**2
                    src_flux_U_err_sq += flux_U_err**2
                    src_flux_V_err_sq += flux_V_err**2

                    # Calculate and store polarisation fractions and angle for each Gaussian in the island
                    # For this we need the I flux, which we can just take from g.total_flux and src.total_flux
                    flux_I = gaussian.total_flux
                    flux_I_err = gaussian.total_fluxE
                    stokes = [flux_I, flux_Q, flux_U, flux_V]
                    stokes_err = [flux_I_err, flux_Q_err, flux_U_err, flux_V_err]

                    lpol_frac, lpol_frac_loerr, lpol_frac_hierr = self.calc_lpol_fraction(stokes, stokes_err) # linear pol fraction
                    lpol_ang, lpol_ang_err = self.calc_lpol_angle(stokes, stokes_err) # linear pol angle
                    cpol_frac, cpol_frac_loerr, cpol_frac_hierr = self.calc_cpol_fraction(stokes, stokes_err) # circular pol fraction
                    tpol_frac, tpol_frac_loerr, tpol_frac_hierr = self.calc_tpol_fraction(stokes, stokes_err) # total pol fraction

                    gaussian.lpol_fraction = lpol_frac
                    gaussian.lpol_fraction_loerr = lpol_frac_loerr
                    gaussian.lpol_fraction_hierr = lpol_frac_hierr
                    gaussian.cpol_fraction = cpol_frac
                    gaussian.cpol_fraction_loerr = cpol_frac_loerr
                    gaussian.cpol_fraction_hierr = cpol_frac_hierr
                    gaussian.tpol_fraction = tpol_frac
                    gaussian.tpol_fraction_loerr = tpol_frac_loerr
                    gaussian.tpol_fraction_hierr = tpol_frac_hierr
                    gaussian.lpol_angle = lpol_ang
                    gaussian.lpol_angle_err = lpol_ang_err

                # Store fluxes for each source in the island
                if hasattr(src, '_pi'):
                    src.total_flux = src_flux_I
                    src.total_fluxE = N.sqrt(src_flux_I_err_sq)
                src.total_flux_Q = src_flux_Q
                src.total_flux_U = src_flux_U
                src.total_flux_V = src_flux_V
                src.total_fluxE_Q = N.sqrt(src_flux_Q_err_sq)
                src.total_fluxE_U = N.sqrt(src_flux_U_err_sq)
                src.total_fluxE_V = N.sqrt(src_flux_V_err_sq)

                # Calculate and store polarisation fractions and angle for each source in the island
                # For this we need the I flux, which we can just take from g.total_flux and src.total_flux
                src_flux_I = src.total_flux
                src_flux_I_err = src.total_fluxE
                stokes = [src_flux_I, src_flux_Q, src_flux_U, src_flux_V]
                stokes_err = [src_flux_I_err, N.sqrt(src_flux_Q_err_sq), N.sqrt(src_flux_U_err_sq), N.sqrt(src_flux_V_err_sq)]

                lpol_frac, lpol_frac_loerr, lpol_frac_hierr = self.calc_lpol_fraction(stokes, stokes_err) # linear pol fraction
                lpol_ang, lpol_ang_err = self.calc_lpol_angle(stokes, stokes_err) # linear pol angle
                cpol_frac, cpol_frac_loerr, cpol_frac_hierr = self.calc_cpol_fraction(stokes, stokes_err) # circular pol fraction
                tpol_frac, tpol_frac_loerr, tpol_frac_hierr = self.calc_tpol_fraction(stokes, stokes_err) # total pol fraction

                src.lpol_fraction = lpol_frac
                src.lpol_fraction_loerr = lpol_frac_loerr
                src.lpol_fraction_hierr = lpol_frac_hierr
                src.cpol_fraction = cpol_frac
                src.cpol_fraction_loerr = cpol_frac_loerr
                src.cpol_fraction_hierr = cpol_frac_hierr
                src.tpol_fraction = tpol_frac
                src.tpol_fraction_loerr = tpol_frac_loerr
                src.tpol_fraction_hierr = tpol_frac_hierr
                src.lpol_angle = lpol_ang
                src.lpol_angle_err = lpol_ang_err
                if bar.started:
                    bar.increment()
          bar.stop()
          img.completed_Ops.append('polarisation')

  ####################################################################################
    def calc_lpol_fraction(self, stokes, err):
        """ Calculate linear polarisation fraction and error from:
            stokes = [I, Q, U, V] and err = [Ierr, Qerr, Uerr, Verr]

        """
        I, Q, U, V = stokes
        Ierr, Qerr, Uerr, Verr = err
        QUerr = N.mean([Qerr, Uerr])
        stokes_lpol = [I, Q, U, 0.0]
        err_lpol = [Ierr, Qerr, Uerr, 0.0]

        lfrac, loerr, uperr, Iup, Qup, Uup, Vup = self.estimate_err_frac_with_limits(stokes_lpol, err_lpol)

        # If all are detections, debias and use error propagation instead
        if not Iup and not Qup and not Uup:
            lpol = N.sqrt(Q**2 + U**2)
            lpol_debiased = self.debias(lpol, QUerr) # debias (to first order)
            if lpol_debiased > 0.0:
                lfrac = lpol_debiased / I
                dlfrac = lfrac * N.sqrt((Ierr/I)**2 + (Q*Qerr/lpol_debiased**2)**2 + (U*Uerr/lpol_debiased**2)**2)
            else:
                # if debiased fraction is consistent with zero, estimate a ballpark error with biased value
                lfrac = 0.0
                lpolsq = Q**2 + U**2
                dlfrac = N.sqrt(lpolsq) / I * N.sqrt((Ierr/I)**2 + (Q*Qerr/lpolsq)**2 + (U*Uerr/lpolsq)**2)
            loerr = dlfrac
            uperr = dlfrac

        lfrac, loerr, uperr = self.check_frac(lfrac, loerr, uperr)
        return lfrac, loerr, uperr


  ####################################################################################
    def calc_cpol_fraction(self, stokes, err):
        """ Calculate circular polarisation fraction and error from:
            stokes = [I, Q, U, V] and err = [Ierr, Qerr, Uerr, Verr]

        """
        I, Q, U, V = stokes
        Ierr, Qerr, Uerr, Verr = err
        stokes_cpol = [I, 0.0, 0.0, V]
        err_cpol = [Ierr, 0.0, 0.0, Verr]

        cfrac, loerr, uperr, Iup, Qup, Uup, Vup = self.estimate_err_frac_with_limits(stokes_cpol, err_cpol)

        # If all are detections, debias and use error propagation instead
        if not Iup and not Vup:
            cfrac = abs(V) / I
            dcfrac = cfrac * N.sqrt((Ierr/I)**2 + (Verr/V)**2)
            loerr = dcfrac
            uperr = dcfrac

        cfrac, loerr, uperr = self.check_frac(cfrac, loerr, uperr)
        return cfrac, loerr, uperr


  ####################################################################################
    def calc_tpol_fraction(self, stokes, err):
        """ Calculate total polarisation fraction and error from:
            stokes = [I, Q, U, V] and err = [Ierr, Qerr, Uerr, Verr]

        """
        I, Q, U, V = stokes
        Ierr, Qerr, Uerr, Verr = err
        QUerr = N.mean([Qerr, Uerr])

        tfrac, loerr, uperr, Iup, Qup, Uup, Vup = self.estimate_err_frac_with_limits(stokes, err)

        # If all are detections, debias and use error propagation instead
        if not Iup and not Qup and not Uup and not Vup:
            lpol = N.sqrt(Q**2 + U**2)
            lpol_debiased = self.debias(lpol, QUerr)
            tpol_debiased = N.sqrt(Q**2 + U**2 + V**2) - (lpol - lpol_debiased) # debias (to first order)
            if tpol_debiased > 0.0:
                tfrac = tpol_debiased / I
                dtfrac = tfrac * N.sqrt((Ierr/I)**2 + (Q*Qerr/tpol_debiased**2)**2 + (U*Uerr/tpol_debiased**2)**2 + (V*Verr/tpol_debiased**2)**2)
            else:
                # if debiased fraction is consistent with zero, estimate a ballpark error with biased value
                tfrac = 0.0
                tpolsq = Q**2 + U**2 + V**2
                dtfrac = N.sqrt(tpolsq) / I * N.sqrt((Ierr/I)**2 + (Q*Qerr/tpolsq)**2 + (U*Uerr/tpolsq)**2 + (V*Verr/tpolsq)**2)
            loerr = dtfrac
            uperr = dtfrac

        tfrac, loerr, uperr = self.check_frac(tfrac, loerr, uperr)
        return tfrac, loerr, uperr


  ####################################################################################
    def calc_lpol_angle(self, stokes, err, sig=3.0):
        """ Calculate linear polarisation angle and error (in degrees) from:
            stokes = [I, Q, U, V] and err = [Ierr, Qerr, Uerr, Verr]

        """
        I, Q, U, V = stokes
        Ierr, Qerr, Uerr, Verr = err
        if abs(Q) < sig*abs(Qerr) and abs(U) < sig*abs(Uerr):
            return 0.0, 0.0

        ang = 0.5 * N.arctan2(U, Q) * 180.0 / N.pi
        dang = 0.5 / (1.0 + (U/Q)**2) * N.sqrt((Uerr/Q)**2 + (U*Qerr/Q**2)**2) * 180.0 / N.pi

        return ang, dang


  ####################################################################################
    def debias(self, pflux, QUerr):
        """ Debiases the linearly polarised flux using the same method
            used for the NVSS catalog (see ftp://ftp.cv.nrao.edu/pub/nvss/catalog.ps).

        """
        data_table=N.array([[1.253,1.2530], [1.256,1.1560], [1.266,1.0660], [1.281,0.9814],
                            [1.303,0.9030], [1.330,0.8304], [1.364,0.7636], [1.402,0.7023],
                            [1.446,0.6462], [1.495,0.5951], [1.549,0.5486], [1.606,0.5064],
                            [1.668,0.4683], [1.734,0.4339], [1.803,0.4028], [1.875,0.3749],
                            [1.950,0.3498], [2.027,0.3273], [2.107,0.3070], [2.189,0.2888],
                            [2.272,0.2724], [2.358,0.2576], [2.444,0.2442], [2.532,0.2321],
                            [2.621,0.2212], [2.711,0.2112], [2.802,0.2021], [2.894,0.1938],
                            [2.986,0.1861], [3.079,0.1791], [3.173,0.1726], [3.267,0.1666],
                            [3.361,0.1610], [3.456,0.1557], [3.551,0.1509], [3.646,0.1463],
                            [3.742,0.1420], [3.838,0.1380], [3.934,0.1342], [4.031,0.1306]])

        pnorm = pflux / QUerr
        if pnorm <= data_table[0,0]:
            bias = data_table[0,1]
        else:
            if pnorm >= data_table[-1,0]:
                bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm**3)
                pnorm = pnorm - bias
                bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm**3)
            else:
                bias = N.interp(pnorm, data_table[:,0], data_table[:,1])

        pflux_debiased = pflux - bias * QUerr

        return pflux_debiased

    def check_frac(self, frac, loerr, uperr):
        if frac < 0.0:
            frac = 0.0
        if frac > 1.0:
            frac = 1.0
        if loerr < 0.0:
            loerr = frac
        if frac + uperr > 1.0:
            uperr = 1.0 - frac
        return frac, loerr, uperr

  ####################################################################################
    def setpara_bdsm(self, img):
        from types import ClassType, TypeType

        chain=[Op_preprocess, Op_rmsimage(), Op_threshold(), Op_islands(),
               Op_gausfit(), Op_gaul2srl(), Op_make_residimage()]

        opts = img.opts.to_dict()
        if img.opts.pi_thresh_isl != None:
            opts['thresh_isl'] = img.opts.pi_thresh_isl
        if img.opts.pi_thresh_pix != None:
            opts['thresh_pix'] = img.opts.pi_thresh_pix
        opts['thresh'] = 'hard'
        opts['polarisation_do'] = False
        opts['filename'] = ''
        opts['detection_image'] = ''

        ops = []
        for op in chain:
          if isinstance(op, (ClassType, TypeType)):
            ops.append(op())
          else:
            ops.append(op)

        return ops, opts

    def estimate_err_frac_with_limits(self, stokes, err, sig=3.0):
        """Estimate reasonable errors on polarization fraction when upper
        limits are present.

        """
        I, Q, U, V = stokes
        Ierr, Qerr, Uerr, Verr = err

        Iup = False
        Qup = False
        Uup = False
        Vup = False

        if abs(I) < sig * abs(Ierr):
            Iup = True
        if abs(Q) < sig * abs(Qerr):
            Q = 0.0
            Qup = True
        if abs(U) < sig * abs(Uerr):
            U = 0.0
            Uup = True
        if abs(V) < sig * abs(Verr):
            V = 0.0
            Vup = True

        pol = N.sqrt(Q**2 + U**2 + V**2)
        frac = pol / I
        if frac < 0.0:
            frac = 0.0
        if frac > 1.0:
            frac = 1.0

        if Iup:
            if Qup and Uup and Vup:
                frac = 0.0
                loerr = 0.0
                uperr = 1.0
            else:
                loerr = frac - N.sqrt((abs(Q) - Qerr)**2 + (abs(U) - Uerr)**2 + (abs(V) - Verr)**2) / abs(Ierr)
                uperr = 1.0 - frac
        else:
            loerr = frac - N.sqrt((abs(Q) - Qerr)**2 + (abs(U) - Uerr)**2 + (abs(V) - Verr)**2) / (I + Ierr)
            uperr = N.sqrt((abs(Q) + Qerr)**2 + (abs(U) + Uerr)**2 + (abs(V) + Verr)**2) / (I - Ierr) - frac

        if loerr < 0.0:
            loerr = frac
        if frac + uperr > 1.0:
            uperr = 1.0 - frac

        return frac, loerr, uperr, Iup, Qup, Uup, Vup


    def double_bbox(self, bbox, shape):
        """Expand bbox of the island by factor of 2

        bbox is isl.bbox
        shape is img.shape
        """
        def expand(bbox, shape):
            bbox_width = (bbox.stop - bbox.start)/2.0
            return slice(max(0, bbox.start - bbox_width), min(shape, bbox.stop + bbox_width))
        return map(expand, bbox, shape)
