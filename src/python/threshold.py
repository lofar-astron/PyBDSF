"""Module threshold.

Defines operation Op_threshold. If the option 'thresh' is defined
as 'fdr' then the value of thresh_pix is estimated using the
False Detection Rate algorithm (using the user defined value
of fdr_alpha). If thresh is None, then the false detection
probability is first calculated, and if the number of false source
pixels is more than fdr_ratio times the estimated number of true source
pixels, then FDR is chosen, else the hard threshold option is chosen.

Masked images aren't handled properly yet.
"""

import numpy as N
from image import Op, Image, NArray
from math import sqrt,pi,log
from scipy.special import erfc
import const
import mylogger


class Op_threshold(Op):
    """Calculates FDR threshold if necessary.

    Prerequisites: Module preprocess and rmsimage should be run first.
    """
    def __call__(self, img):
        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Threshold ")
        data = img.ch0
        mask = img.mask
        opts = img.opts
        size = N.product(img.ch0.shape)
        sq2  = sqrt(2)

	if img.opts.thresh is None:
            source_p = self.get_srcp(img)
	    cutoff = 5.0
	    false_p = 0.5*erfc(cutoff/sq2)*size
	    if false_p < opts.fdr_ratio*source_p:
                img.thresh = 'hard'
                mylogger.userinfo(mylog, "Expected 5-sigma-clipped false detection rate < fdr_ratio")
                mylogger.userinfo(mylog, "Using sigma-clipping ('hard') thresholding")
	    else:
                img.thresh = 'fdr'
                mylogger.userinfo(mylog, "Expected 5-sigma-clipped false detection rate > fdr_ratio")
                mylogger.userinfo(mylog, "Using FDR (False Detection Rate) thresholding")
            mylog.debug('%s %g' % ("Estimated number of source pixels (using sourcecounts.py) is ",source_p))
            mylog.debug('%s %g' % ("Number of false positive pixels expected for 5-sigma is ",false_p))
            mylog.debug("Threshold for pixels set to : "+str.swapcase(img.thresh))
        else:
            img.thresh = img.opts.thresh

	if img.thresh=='fdr':
            cdelt = img.wcs_obj.acdelt[:2]
	    bm = (img.beam[0], img.beam[1])
            area_pix = int(round(N.product(bm)/(abs(N.product(cdelt))* \
                                                  pi/(4.0*log(2.0)))))
	    s0 = 0
	    for i in range(area_pix):
                s0 +=  1.0/(i+1)
	    slope = opts.fdr_alpha/s0
                                                # sort erf of normalised image as vector
	    v = N.sort(0.5*erfc(N.ravel((data-img.mean)/img.rms)/sq2))[::-1]
	    for i,x in enumerate(v):
                if x < slope*i/size:
                    pcrit = x
                    break
            dumr1 = 1.0-2.0*pcrit;
            dumr = 8.0/3.0/pi*(pi-3.0)/(4.0-pi)
			    # approx for inv(erfc)
	    sigcrit = sqrt(-2.0/pi/dumr-log(1.0-dumr1*dumr1)/2.0+  \
                      sqrt((2.0/pi/dumr+log(1.0-dumr1*dumr1)/2.0)* \
                      (2.0/pi/dumr+log(1.0-dumr1*dumr1)/2.0)-      \
                      log(1.0-dumr1*dumr1)/dumr))*sq2
            if pcrit == 0.0:
                img.thresh = 'hard'
            else:
                img.thresh_pix = sigcrit
                mylogger.userinfo(mylog, "FDR threshold (replaces thresh_pix)", str(round(sigcrit, 4)))
        else:
            img.thresh_pix = opts.thresh_pix

        img.completed_Ops.append('threshold')
        return img

    def get_srcp(self, img):
        import sourcecounts as sc
        fwsig = const.fwsig
	cutoff = 5.0
	spin = -0.80
	freq = img.frequency
        bm = (img.beam[0], img.beam[1])
        cdelt = img.wcs_obj.acdelt[:2]
	x = 2.0*pi*N.product(bm)/abs(N.product(cdelt))/(fwsig*fwsig)*img.omega

	smin_L = img.clipped_rms*cutoff*((1.4e9/freq)**spin)
	scflux = sc.s
	scnum = sc.n
        index = 0
	for i,s in enumerate(scflux):
            if s < smin_L:
                index = i
                break
	n1 = scnum[index]; n2 = scnum[-1]
	s1 = scflux[index]; s2 = scflux[-1]
	alpha = 1.0-log(n1/n2)/log(s1/s2)
	A = (alpha-1.0)*n1/(s1**(1.0-alpha))
	source_p = x*A*((cutoff*img.clipped_rms)**(1.0-alpha)) \
                 /((1.0-alpha)*(1.0-alpha))

        return source_p
