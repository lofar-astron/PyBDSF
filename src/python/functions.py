# some functions

def poly(c,x):
    """ y = Sum { c(i)*x^i }, i=0,len(c)"""
    import numpy as N
    y=N.zeros(len(x))
    for i in range(len(c)):
        y += c[i]*(x**i)
    return y

def sp_in(c, x):
    """ Spectral index in freq-flux space """
    import numpy as N

    order = len(c)-1
    if order == 1:
      y = c[0]*N.power(x, c[1])
    else:
      if order == 2:
        y = c[0]*N.power(x, c[1])*N.power(x, c[2]*N.log(x))
      else:
        print 'Not yet implemented'

    return y

def wenss_fit(c,x):
    """ sqrt(c0*c0 + c1^2/x^2)"""
    import numpy as N
    y = N.sqrt(c[0]*c[0]+c[1]*c[1]/(x*x))
    return y

def nanmean(x):
    """ Mean of array with NaN """
    import numpy as N

    sum = N.nansum(x)
    n = N.sum(~N.isnan(x))

    if n > 0:
      mean = sum/n
    else:
      mean = float("NaN")

    return mean

def shapeletfit(cf, Bset, cfshape):
    """ The function """
    import numpy as N

    ordermax = Bset.shape[0]
    y = (Bset[0,0,::]).flatten()
    y = N.zeros(y.shape)
    index = [(i,j) for i in range(ordermax) for j in range(ordermax-i)]  # i=0->nmax, j=0-nmax-i
    for coord in index:
	linbasis = (Bset[coord[0], coord[1], ::]).flatten()
        y += cf.reshape(cfshape)[coord]*linbasis

    return y

def func_poly2d(ord,p,x,y):
    """ 2d polynomial.
    ord=0 : z=p[0]
    ord=1 : z=p[0]+p[1]*x+p[2]*y
    ord=2 : z=p[0]+p[1]*x+p[2]*y+p[3]*x*x+p[4]*y*y+p[5]*x*y
    ord=3 : z=p[0]+p[1]*x+p[2]*y+p[3]*x*x+p[4]*y*y+p[5]*x*y+
              p[6]*x*x*x+p[7]*x*x*y+p[8]*x*y*y+p[9]*y*y*y"""

    if ord == 0:
        z=p[0]
    if ord == 1:
        z=p[0]+p[1]*x+p[2]*y
    if ord == 2:
        z=p[0]+p[1]*x+p[2]*y+p[3]*x*x+p[4]*y*y+p[5]*x*y
    if ord == 3:
        z=p[0]+p[1]*x+p[2]*y+p[3]*x*x+p[4]*y*y+p[5]*x*y+\
          p[6]*x*x*x+p[7]*x*x*y+p[8]*x*y*y+p[9]*y*y*y
    if ord > 3:
        print " We do not trust polynomial fits > 3 "
	z = None

    return z

def func_poly2d_ini(ord, av):
    """ Initial guess -- assume flat plane. """

    if ord == 0:
        p0 = N.asarray([av])
    if ord == 1:
        p0 = N.asarray([av] + [0.0]*2)
    if ord == 2:
        p0 = N.asarray([av] + [0.0]*5)
    if ord == 3:
        p0 = N.asarray([av] + [0.0]*9)
    if ord > 3:
        p0 = None

    return p0

def ilist(x):
    """ integer part of a list of floats. """

    fn = lambda x : [int(round(i)) for i in x]
    return fn(x)

def cart2polar(cart, cen):
    """ convert cartesian coordinates to polar coordinates around cen. theta is
    zero for +ve xaxis and goes counter clockwise. cart is a numpy array [x,y] where
    x and y are numpy arrays of all the (>0) values of coordinates."""
    import math

    polar = N.zeros(cart.shape)
    pi = math.pi
    rad = 180.0/pi

    cc = N.transpose(cart)
    cc = (cc-cen)*(cc-cen)
    polar[0] = N.sqrt(N.sum(cc,1))
    th = N.arctan2(cart[1]-cen[1],cart[0]-cen[0])*rad
    polar[1] = N.where(th > 0, th, 360+th)

    return polar


def polar2cart(polar, cen):
    """ convert polar coordinates around cen to cartesian coordinates. theta is
    zero for +ve xaxis and goes counter clockwise. polar is a numpy array of [r], [heta]
    and cart is a numpy array [x,y] where x and y are numpy arrays of all the (>0)
    values of coordinates."""
    import math

    cart = N.zeros(polar.shape)
    pi = math.pi
    rad = 180.0/pi

    cart[0]=polar[0]*N.cos(polar[1]/rad)+cen[0]
    cart[1]=polar[0]*N.sin(polar[1]/rad)+cen[1]

    return cart

def gaus_pixval(g, pix):
    """ Calculates the value at a pixel pix due to a gaussian object g. """
    from const import fwsig, pi
    from math import sin, cos, exp

    cen = g.centre_pix
    peak = g.peak_flux
    bmaj_p, bmin_p, bpa_p = g.size_pix

    a4 = bmaj_p/fwsig; a5 = bmin_p/fwsig
    a6 = (bpa_p+90.0)*pi/180.0
    spa = sin(a6); cpa = cos(a6)
    dr1 = ((pix[0]-cen[0])*cpa + (pix[1]-cen[1])*spa)/a4
    dr2 = ((pix[1]-cen[1])*cpa - (pix[0]-cen[0])*spa)/a5
    pixval = peak*exp(-0.5*(dr1*dr1+dr2*dr2))

    return pixval

def atanproper(dumr, dx, dy):
    from math import pi

    ysign = (dy >= 0.0)
    xsign = (dx >= 0.0)
    if ysign and (not xsign): dumr = pi - dumr
    if (not ysign) and (not xsign): dumr = pi + dumr
    if (not ysign) and xsign: dumr = 2.0*pi - dumr

    return dumr

def gdist_pa(pix1, pix2, gsize):
    """ Computes FWHM in degrees in the direction towards second source, of an elliptical gaussian. """
    from math import atan, pi, sqrt, cos, sin, tan

    dx = pix2[0] - pix1[0]
    dy = pix2[1] - pix1[1]
    if dx == 0.0:
        val = pi/2.0
    else:
        dumr = atan(abs(dy/dx))
        val = atanproper(dumr, dx, dy)

    psi = val - (gsize[2]+90.0)/180.0*pi

    # convert angle to eccentric anomaly
    if approx_equal(gsize[1], 0.0):
        psi = pi/2.0
    else:
        psi=atan(gsize[0]/gsize[1]*tan(psi))
    dumr2 = gsize[0]*cos(psi)
    dumr3 = gsize[1]*sin(psi)
    fwhm = sqrt(dumr2*dumr2+dumr3*dumr3)

    return fwhm

def gaus_2d(c, x, y):
    """ x and y are 2d arrays with the x and y positions. """
    import math
    import numpy as N

    rad = 180.0/math.pi
    cs = math.cos(c[5]/rad)
    sn = math.sin(c[5]/rad)
    f1 = ((x-c[1])*cs+(y-c[2])*sn)/c[3]
    f2 = ((y-c[2])*cs-(x-c[1])*sn)/c[4]
    val = c[0]*N.exp(-0.5*(f1*f1+f2*f2))

    return val

def gaus_2d_itscomplicated(c, x, y, p_tofix, ind):
    """ x and y are 2d arrays with the x and y positions. c is a list (of lists) of gaussian parameters to fit, p_tofix
    are gaussian parameters to fix. ind is a list with 0, 1; 1 = fit; 0 = fix. """

    import math
    import numpy as N

    val = N.zeros(x.shape)
    indx = N.array(ind)
    if len(indx) % 6 != 0:
      print " Something wrong with the parameters passed - need multiples of 6 !"
    else:
      ngaus = len(indx)/6
      params = N.zeros(6*ngaus)
      params[N.where(indx==1)[0]] = c
      params[N.where(indx==0)[0]] = p_tofix
      for i in range(ngaus):
        gau = params[i*6:i*6+6]
        val = val + gaus_2d(gau, x, y)

    return val

def g2param(g, adj=False):
    """Convert gaussian object g to param list [amp, cenx, ceny, sigx, sigy, theta] """
    from const import fwsig
    from math import pi

    A = g.peak_flux
    if adj and hasattr(g, 'size_pix_adj'):
        sigx, sigy, th = g.size_pix_adj
    else:
        sigx, sigy, th = g.size_pix
    cenx, ceny = g.centre_pix
    sigx = sigx/fwsig; sigy = sigy/fwsig; th = th+90.0
    params = [A, cenx, ceny, sigx, sigy, th]

    return params

def g2param_err(g, adj=False):
    """Convert errors on gaussian object g to param list [Eamp, Ecenx, Eceny, Esigx, Esigy, Etheta] """
    from const import fwsig
    from math import pi

    A = g.peak_fluxE
    if adj and hasattr(g, 'size_pix_adj'):
        sigx, sigy, th = g.size_pix_adj
    else:
        sigx, sigy, th = g.size_pixE
    cenx, ceny = g.centre_pixE
    sigx = sigx/fwsig; sigy = sigy/fwsig
    params = [A, cenx, ceny, sigx, sigy, th]

    return params

def corrected_size(size):
    """ convert major and minor axis from sigma to fwhm and angle from horizontal to P.A. """

    from const import fwsig

    csize = [0,0,0]
    csize[0] = size[0]*fwsig
    csize[1] = size[1]*fwsig
    bpa = size[2]
    pa = bpa-90.0
    pa = pa % 360
    if pa < 0.0: pa = pa + 360.0
    if pa > 180.0: pa = pa - 180.0
    csize[2] = pa

    return csize

def drawellipse(g):
    import math
    import numpy as N
    from gausfit import Gaussian

    rad = 180.0/math.pi
    if isinstance(g, Gaussian):
      param = g2param(g)
    else:
      if isinstance(g, list) and len(g)>=6:
        param = g
      else:
        raise RuntimeError("Input to drawellipse neither Gaussian nor list")

    x2 = []; y2 = []
    size = [param[3], param[4], param[5]]
    size_fwhm = corrected_size(size)
    for th in range(0, 370, 10):
      x1=size_fwhm[0]*math.cos(th/rad)
      y1=size_fwhm[1]*math.sin(th/rad)
      x2.append(x1*math.cos(param[5]/rad)-y1*math.sin(param[5]/rad)+param[1])
      y2.append(x1*math.sin(param[5]/rad)+y1*math.cos(param[5]/rad)+param[2])
    x2 = N.array(x2); y2 = N.array(y2)

    return x2, y2

def drawsrc(src):
    import math
    import numpy as N
    import matplotlib.path as mpath
    Path = mpath.Path
    paths = []
    xmin = []
    xmax = []
    ymin = []
    ymax = []
    ellx = []
    elly = []
    for indx, g in enumerate(src.gaussians):
        gellx, gelly = drawellipse(g)
        ellx += gellx.tolist()
        elly += gelly.tolist()
    yarr = N.array(elly)
    minyarr = N.min(yarr)
    maxyarr = N.max(yarr)
    xarr = N.array(ellx)
    for i in range(10):
        inblock = N.where(yarr > minyarr + float(i)*(maxyarr-minyarr)/10.0)
        yarr = yarr[inblock]
        xarr = xarr[inblock]
        inblock = N.where(yarr < minyarr + float(i+1)*(maxyarr-minyarr)/10.0)
        xmin.append(N.min(xarr[inblock])-1.0)
        xmax.append(N.max(xarr[inblock])+1.0)
        ymin.append(N.mean(yarr[inblock]))
        ymax.append(N.mean(yarr[inblock]))

    xmax.reverse()
    ymax.reverse()
    pathdata = [(Path.MOVETO, (xmin[0], ymin[0]))]
    for i in range(10):
        pathdata.append((Path.LINETO, (xmin[i], ymin[i])))
        pathdata.append((Path.CURVE3, (xmin[i], ymin[i])))
    pathdata.append((Path.LINETO, ((xmin[9]+xmax[0])/2.0, (ymin[9]+ymax[0])/2.0+1.0)))
    for i in range(10):
        pathdata.append((Path.LINETO, (xmax[i], ymax[i])))
        pathdata.append((Path.CURVE3, (xmax[i], ymax[i])))
    pathdata.append((Path.LINETO, ((xmin[0]+xmax[9])/2.0, (ymin[0]+ymax[9])/2.0-1.0)))
    pathdata.append((Path.CLOSEPOLY, (xmin[0], ymin[0])))
    codes, verts = zip(*pathdata)
    path = Path(verts, codes)
    return path

def mask_fwhm(g, fac1, fac2, delc, shap):
    """ take gaussian object g and make a mask (as True) for pixels which are outside (less flux)
        fac1*FWHM and inside (more flux) fac2*FWHM. Also returns the values as well."""
    import math
    import numpy as N
    from const import fwsig

    x, y = N.indices(shap)
    params = g2param(g)
    params[1] -= delc[0]; params[2] -= delc[1]
    gau = gaus_2d(params, x, y)
    dumr1 = 0.5*fac1*fwsig
    dumr2 = 0.5*fac2*fwsig
    flux1= params[0]*math.exp(-0.5*dumr1*dumr1)
    flux2 = params[0]*math.exp(-0.5*dumr2*dumr2)
    mask = (gau <= flux1) * (gau > flux2)
    gau = gau * mask

    return mask, gau

def flatten(x):
    """flatten(sequence) -> list
    Taken from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def moment(x,mask=None):
    """
    Calculates first 3 moments of numpy array x. Only those values of x
    for which mask is False are used, if mask is given. Works for any
    dimension of x.
    """
    import numpy as N

    if mask == None:
        mask=N.zeros(x.shape, dtype=bool)
    m1=N.zeros(1)
    m2=N.zeros(x.ndim)
    m3=N.zeros(x.ndim)
    for i, val in N.ndenumerate(x):
        if not mask[i]:
            m1 += val
	    m2 += val*N.array(i)
	    m3 += val*N.array(i)*N.array(i)
    m2 /= m1
    m3 = N.sqrt(m3/m1-m2*m2)
    return m1, m2, m3

def fit_mask_1d(x, y, sig, mask, funct, do_err, order=0, p0 = None):
    """
    Calls scipy.optimise.leastsq for a 1d function with a mask.
    Takes values only where mask=False.
    """
    from scipy.optimize import leastsq
    from math import sqrt, pow
    import numpy as N
    import sys

    ind=N.where(~N.array(mask))[0]
    if len(ind) > 1:
      n=sum(mask)

      if isinstance(x, list): x = N.array(x)
      if isinstance(y, list): y = N.array(y)
      if isinstance(sig, list): sig = N.array(sig)
      xfit=x[ind]; yfit=y[ind]; sigfit=sig[ind]

      if p0 == None:
        if funct == poly:
           p0=N.array([0]*(order+1))
           p0[1]=(yfit[0]-yfit[-1])/(xfit[0]-xfit[-1])
           p0[0]=yfit[0]-p0[1]*xfit[0]
        if funct == wenss_fit:
           p0=N.array([yfit[N.argmax(xfit)]] + [1.])
        if funct == sp_in:
           ind1 = N.where(yfit > 0.)[0]
           if len(ind1) >= 2:
             low = ind1[0]; hi = ind1[-1]
             sp = N.log(yfit[low]/yfit[hi])/N.log(xfit[low]/xfit[hi])
             p0=N.array([yfit[low]/pow(xfit[low], sp), sp] + [0.]*(order-1))
           elif len(ind1) == 1:
             p0=N.array([ind1[0], -0.8] + [0.]*(order-1))
           else:
             return [0, 0], [0, 0]
      res=lambda p, xfit, yfit, sigfit: (yfit-funct(p, xfit))/sigfit
      try:
        (p, cov, info, mesg, flag)=leastsq(res, p0, args=(xfit, yfit, sigfit), full_output=True, warning=False)
      except TypeError:
        # This error means no warning argument is available, so redirect stdout to a null device
        # to suppress printing of (unnecessary) warning messages
        original_stdout = sys.stdout  # keep a reference to STDOUT
        sys.stdout = NullDevice()  # redirect the real STDOUT
        (p, cov, info, mesg, flag)=leastsq(res, p0, args=(xfit, yfit, sigfit), full_output=True)
        sys.stdout = original_stdout  # turn STDOUT back on

      if do_err:
        if cov != None:
          if N.sum(sig != 1.) > 0:
            err = N.array([sqrt(abs(cov[i,i])) for i in range(len(p))])
          else:
            chisq=sum(info["fvec"]*info["fvec"])
            dof=len(info["fvec"])-len(p)
            err = N.array([sqrt(abs(cov[i,i])*chisq/dof) for i in range(len(p))])
        else:
          p, err = [0, 0], [0, 0]
      else: err = [0]
    else:
      p, err = [0, 0], [0, 0]

    return p, err

def dist_2pt(p1, p2):
    """ Calculated distance between two points given as tuples p1 and p2. """
    from math import sqrt
    dx=p1[0]-p2[0]
    dy=p1[1]-p2[1]
    dist=sqrt(dx*dx + dy*dy)

    return dist

def std(y):
    """ Returns unbiased standard deviation. """
    from math import sqrt
    import numpy as N

    l=len(y)
    s=N.std(y)
    if l == 1:
        return s
    else:
        return s*sqrt(float(l)/(l-1))

def imageshift(image, shift):
    """ Shifts a 2d-image by the tuple (shift). Positive shift is to the right and upwards.
    This is done by fourier shifting. """
    import scipy
    from scipy import ndimage

    shape=image.shape

    f1=scipy.fft(image, shape[0], axis=0)
    f2=scipy.fft(f1, shape[1], axis=1)

    s=ndimage.fourier_shift(f2,shift, axis=0)

    y1=scipy.ifft(s, shape[1], axis=1)
    y2=scipy.ifft(y1, shape[0], axis=0)

    return y2.real

def trans_gaul(q):
    " transposes a tuple "
    y=[]
    if len(q) > 0:
        for i in range(len(q[0])):
            elem=[]
            for j in range(len(q)):
                elem.append(q[j][i])
            y.append(elem)
    return y

def momanalmask_gaus(subim, mask, isrc, bmar_p, allpara=True):
    """ Compute 2d gaussian parameters from moment analysis, for an island with
        multiple gaussians. Compute only for gaussian with index (mask value) isrc.
        Returns normalised peak, centroid, fwhm and P.A. assuming North is top.
    """
    from math import sqrt, atan, pi
    from const import fwsig
    import numpy as N
    N.seterr(all='ignore')

    m1 = N.zeros(2); m2 = N.zeros(2); m11 = 0.0; tot = 0.0
    mompara = N.zeros(6)
    n, m = subim.shape[0], subim.shape[1]
    index = [(i, j) for i in range(n) for j in range(m) if mask[i,j]==isrc]
    for coord in index:
        tot += subim[coord]
        m1 += N.array(coord)*subim[coord]
    mompara[0] = tot/bmar_p
    mompara[1:3] = m1/tot

    if allpara:
      for coord in index:
          co = N.array(coord)
          m2 += (co - mompara[1:3])*(co - mompara[1:3])*subim[coord]
          m11 += N.product(co - mompara[1:3])*subim[coord]

      mompara[3] = sqrt((m2[0]+m2[1]+sqrt((m2[0]-m2[1])*(m2[0]-m2[1])+4.0*m11*m11))/(2.0*tot))*fwsig
      mompara[4] = sqrt((m2[0]+m2[1]-sqrt((m2[0]-m2[1])*(m2[0]-m2[1])+4.0*m11*m11))/(2.0*tot))*fwsig
      dumr = atan(abs(2.0*m11/(m2[0]-m2[1])))
      dumr = atanproper(dumr, m2[0]-m2[1], 2.0*m11)
      mompara[5] = 0.5*dumr*180.0/pi - 90.0
      if mompara[5] < 0.0: mompara[5] += 180.0
    return mompara

def fit_gaus2d(data, p_ini, x, y, mask = None, err = None):
    """ Fit 2d gaussian to data with x and y also being 2d numpy arrays with x and y positions.
        Takes an optional error array and a mask array (True => pixel is masked). """
    from scipy.optimize import leastsq
    import numpy as N
    import sys

    if mask != None and mask.shape != data.shape:
        print 'Data and mask array dont have the same shape, ignoring mask'
        mask = None
    if err != None and err.shape != data.shape:
        print 'Data and error array dont have the same shape, ignoring error'
        err = None

    if mask == None: mask = N.zeros(data.shape, bool)
    g_ind = N.where(~N.ravel(mask))[0]

    if err == None:
        errorfunction = lambda p: N.ravel(gaus_2d(p, x, y) - data)[g_ind]
    else:
        errorfunction = lambda p: N.ravel((gaus_2d(p, x, y) - data)/err)[g_ind]
    try:
        p, success = leastsq(errorfunction, p_ini, warning=False)
    except TypeError:
        # This error means no warning argument is available, so redirect stdout to a null device
        # to suppress printing of warning messages
        original_stdout = sys.stdout  # keep a reference to STDOUT
        sys.stdout = NullDevice()  # redirect the real STDOUT
        p, success = leastsq(errorfunction, p_ini)
        sys.stdout = original_stdout  # turn STDOUT back on


    return p, success

def deconv(gaus_bm, gaus_c):
    """ Deconvolves gaus_bm from gaus_c to give gaus_dc.
        Stolen shamelessly from aips DECONV.FOR.
        All PA is in degrees."""
    from math import pi, cos, sin, atan, sqrt

    rad = 180.0/pi
    gaus_d = [0.0, 0.0, 0.0]

    phi_c = gaus_c[2]+900.0 % 180
    phi_bm = gaus_bm[2]+900.0 % 180
    maj2_bm = gaus_bm[0]*gaus_bm[0]; min2_bm = gaus_bm[1]*gaus_bm[1]
    maj2_c = gaus_c[0]*gaus_c[0]; min2_c = gaus_c[1]*gaus_c[1]
    theta=2.0*(phi_c-phi_bm)/rad
    cost = cos(theta)
    sint = sin(theta)

    rhoc = (maj2_c-min2_c)*cost-(maj2_bm-min2_bm)
    if rhoc == 0.0:
      sigic = 0.0
      rhoa = 0.0
    else:
      sigic = atan((maj2_c-min2_c)*sint/rhoc)   # in radians
      rhoa = ((maj2_bm-min2_bm)-(maj2_c-min2_c)*cost)/(2.0*cos(sigic))

    gaus_d[2] = sigic*rad/2.0+phi_bm
    dumr = ((maj2_c+min2_c)-(maj2_bm+min2_bm))/2.0
    gaus_d[0] = dumr-rhoa
    gaus_d[1] = dumr+rhoa
    error = 0
    if gaus_d[0] < 0.0: error += 1
    if gaus_d[1] < 0.0: error += 1

    gaus_d[0] = max(0.0,gaus_d[0])
    gaus_d[1] = max(0.0,gaus_d[1])
    gaus_d[0] = sqrt(abs(gaus_d[0]))
    gaus_d[1] = sqrt(abs(gaus_d[1]))
    if gaus_d[0] < gaus_d[1]:
      sint = gaus_d[0]
      gaus_d[0] = gaus_d[1]
      gaus_d[1] = sint
      gaus_d[2] = gaus_d[2]+90.0

    gaus_d[2] = gaus_d[2]+900.0 % 180
    if gaus_d[0] == 0.0:
      gaus_d[2] = 0.0
    else:
      if gaus_d[1] == 0.0:
        if (abs(gaus_d[2]-phi_c) > 45.0) and (abs(gaus_d[2]-phi_c) < 135.0):
          gaus_d[2] = gaus_d[2]+450.0 % 180

# errors
     #if rhoc == 0.0:
    #if gaus_d[0] != 0.0:
    #  ed_1 = gaus_c[0]/gaus_d[0]*e_1
    #else:
    #  ed_1 = sqrt(2.0*e_1*gaus_c[0])
    #if gaus_d[1] != 0.0:
    #  ed_2 = gaus_c[1]/gaus_d[1]*e_2
    #else:
    #  ed_2 = sqrt(2.0*e_2*gaus_c[1])
    #ed_3 =e_3
    #else:
    #  pass

    return gaus_d

def deconv2(gaus_bm, gaus_c):
    """ Deconvolves gaus_bm from gaus_c to give gaus_dc.
        Stolen shamelessly from Miriad gaupar.for.
        All PA is in degrees.

        Returns deconvolved gaussian parameters and flag:
   	 0   All OK.
     1   Result is pretty close to a point source.
	 2   Illegal result.

        """
    from math import pi, cos, sin, atan2, sqrt

    rad = 180.0/pi

    phi_c = gaus_c[2]+900.0 % 180.0
    phi_bm = gaus_bm[2]+900.0 % 180.0
    theta1 = phi_c / rad
    theta2 = phi_bm / rad
    bmaj1 = gaus_c[0]
    bmaj2 = gaus_bm[0]
    bmin1 = gaus_c[1]
    bmin2 = gaus_bm[1]

    alpha = ( (bmaj1*cos(theta1))**2 + (bmin1*sin(theta1))**2 -
              (bmaj2*cos(theta2))**2 - (bmin2*sin(theta2))**2 )
    beta = ( (bmaj1*sin(theta1))**2 + (bmin1*cos(theta1))**2 -
             (bmaj2*sin(theta2))**2 - (bmin2*cos(theta2))**2 )
    gamma = 2.0 * ( (bmin1**2-bmaj1**2)*sin(theta1)*cos(theta1) -
                  (bmin2**2-bmaj2**2)*sin(theta2)*cos(theta2) )

    s = alpha + beta
    t = sqrt((alpha-beta)**2 + gamma**2)
    limit = min(bmaj1, bmin1, bmaj2, bmin2)
    limit = 0.1*limit*limit

    if alpha < 0.0 or beta < 0.0 or s < t:
        if alpha < 0.0 or beta < 0.0:
            bmaj = 0.0
            bpa = 0.0
        else:
	        bmaj = sqrt(0.5*(s+t))
	        bpa = rad * 0.5 * atan2(-gamma, alpha-beta)
        bmin = 0.0
        if 0.5*(s-t) < limit and alpha > -limit and beta > -limit:
	        ifail = 1
        else:
            ifail = 2
    else:
        bmaj = sqrt(0.5*(s+t))
        bmin = sqrt(0.5*(s-t))
        if abs(gamma) + abs(alpha-beta) == 0.0:
	        bpa = 0.0
        else:
	        bpa = rad * 0.5 * atan2(-gamma, alpha-beta)
        ifail = 0
    return (bmaj, bmin, bpa), ifail


def get_errors(img, p, stdav, bm_pix=None):

    from const import fwsig
    from math import sqrt, log, pow, pi
    import mylogger
    import numpy as N

    mylog = mylogger.logging.getLogger("PyBDSM.Compute")

    if len(p) % 7 > 0:
      mylog.error("Gaussian parameters passed have to have 7n numbers")
    ngaus = len(p)/7
    errors = []
    for i in range(ngaus):
      pp = p[i*7:i*7+7]
      ### Now do error analysis as in Condon (and fBDSM)
      size = pp[3:6]
      size = corrected_size(size) # angle is now degrees CCW from +y-axis
      if size[0] == 0.0 or size[1] == 0.0:
        errors = errors + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      else:
        sq2 = sqrt(2.0)
        if bm_pix == None:
            bm_pix = N.array([img.pixel_beam[0]*fwsig, img.pixel_beam[1]*fwsig, img.pixel_beam[2]])
        dumr = sqrt(abs(size[0]*size[1]/(4.0*bm_pix[0]*bm_pix[1])))
        dumrr1 = 1.0+bm_pix[0]*bm_pix[1]/(size[0]*size[0])
        dumrr2 = 1.0+bm_pix[0]*bm_pix[1]/(size[1]*size[1])
        dumrr3 = dumr*pp[0]/stdav
        d1 = sqrt(8.0*log(2.0))
        d2 = (size[0]*size[0]-size[1]*size[1])/(size[0]*size[0])
        try:
            e_peak = pp[0]*sq2/(dumrr3*pow(dumrr1,0.75)*pow(dumrr2,0.75))
            e_maj=size[0]*sq2/(dumrr3*pow(dumrr1,1.25)*pow(dumrr2,0.25))
            e_min=size[1]*sq2/(dumrr3*pow(dumrr1,0.25)*pow(dumrr2,1.25))  # in fw
            pa_rad = size[2]*pi/180.0
            e_x0 = sqrt( (e_maj*N.sin(pa_rad))**2 + (e_min*N.cos(pa_rad))**2 ) / d1
            e_y0 = sqrt( (e_maj*N.cos(pa_rad))**2 + (e_min*N.sin(pa_rad))**2 ) / d1
            e_pa=2.0/(d2*dumrr3*pow(dumrr1,0.25)*pow(dumrr2,1.25))
            e_pa=e_pa*180.0/pi
            e_tot=pp[0]*sqrt(e_peak*e_peak/(pp[0]*pp[0])+(0.25/dumr/dumr)*(e_maj*e_maj/(size[0]*size[0])+e_min*e_min/(size[1]*size[1])))
        except:
            e_peak = 0.0
            e_x0 = 0.0
            e_y0 = 0.0
            e_maj = 0.0
            e_min = 0.0
            e_pa = 0.0
            e_tot = 0.0
        if abs(e_pa) > 180.0: e_pa=180.0  # dont know why i did this
        errors = errors + [e_peak, e_x0, e_y0, e_maj, e_min, e_pa, e_tot]

    return errors

def fit_chisq(x, p, ep, mask, funct, order):
    import numpy as N

    ind = N.where(N.array(mask)==False)[0]
    if order == 0:
      fit = [funct(p)]*len(p)
    else:
      fitpara, efit = fit_mask_1d(x, p, ep, mask, funct, True, order)
      fit = funct(fitpara, x)

    dev = (p-fit)*(p-fit)/(ep*ep)
    num = order+1
    csq = N.sum(dev[ind])/(len(fit)-num-1)

    return csq

def calc_chisq(x, y, ey, p, mask, funct, order):
    import numpy as N

    if order == 0:
      fit = [funct(y)]*len(y)
    else:
      fit = funct(p, x)

    dev = (y-fit)*(y-fit)/(ey*ey)
    ind = N.where(~N.array(mask))
    num = order+1
    csq = N.sum(dev[ind])/(len(mask)-num-1)

    return csq

def get_windowsize_av(S_i, rms_i, chanmask, K, minchan):
    import numpy as N

    av_window = N.arange(2, int(len(S_i)/minchan)+1)
    win_size = 0
    for window in av_window:
      fluxes, vars, mask = variance_of_wted_windowedmean(S_i, rms_i, chanmask, window)
      minsnr = N.min(fluxes[~mask]/vars[~mask])
      if minsnr > K*1.1:                ### K*1.1 since fitted peak can be less than wted peak
        win_size = window  # is the size of averaging window
        break

    return win_size

def variance_of_wted_windowedmean(S_i, rms_i, chanmask, window_size):
    from math import sqrt
    import numpy as N

    nchan = len(S_i)
    nwin = nchan/window_size
    wt = 1/rms_i/rms_i
    wt = wt/N.median(wt)
    fluxes = N.zeros(nwin); vars = N.zeros(nwin); mask = N.zeros(nwin, bool)
    for i in range(nwin):
      strt = i*window_size; stp = (i+1)*window_size
      if i == nwin-1: stp = nchan
      ind = N.arange(strt,stp)
      m = chanmask[ind]
      index = [arg for ii,arg in enumerate(ind) if not m[ii]]
      if len(index) > 0:
        s = S_i[index]; r = rms_i[index]; w = wt[index]
        fluxes[i] = N.sum(s*w)/N.sum(w)
        vars[i] = 1.0/sqrt(N.sum(1.0/r/r))
        mask[i] = N.product(m)
      else:
        fluxes[i] = 0
        vars[i] = 0
        mask[i] = True

    return fluxes, vars, mask

def fit_mulgaus2d(image, gaus, x, y, mask = None, fitfix = None, err = None, adj=False):
    """ fitcode : 0=fit all; 1=fit amp; 2=fit amp, posn; 3=fit amp, size """
    from scipy.optimize import leastsq
    import numpy as N
    import sys

    if mask != None and mask.shape != image.shape:
        print 'Data and mask array dont have the same shape, ignoring mask'
        mask = None
    if err != None and err.shape != image.shape:
        print 'Data and error array dont have the same shape, ignoring error'
        err = None
    if mask == None: mask = N.zeros(image.shape, bool)

    g_ind = N.where(~N.ravel(mask))[0]

    ngaus = len(gaus)
    if ngaus > 0:
      p_ini = []
      for g in gaus:
        p_ini = p_ini + g2param(g, adj)
      p_ini = N.array(p_ini)

      if fitfix == None: fitfix = [0]*ngaus
      ind = N.ones(6*ngaus)                                     # 1 => fit ; 0 => fix
      for i in range(ngaus):
        if fitfix[i] == 1: ind[i*6+1:i*6+6] = 0
        if fitfix[i] == 2: ind[i*6+3:i*6+6] = 0
        if fitfix[i] == 3: ind[i*6+1:i*6+3] = 0
      ind = N.array(ind)
      p_tofit = p_ini[N.where(ind==1)[0]]
      p_tofix = p_ini[N.where(ind==0)[0]]
      if err == None: err = N.ones(image.shape)

      errorfunction = lambda p, x, y, p_tofix, ind, image, err, g_ind: \
                     N.ravel((gaus_2d_itscomplicated(p, x, y, p_tofix, ind)-image)/err)[g_ind]
      try:
          p, success = leastsq(errorfunction, p_tofit, args=(x, y, p_tofix, ind, image, err, g_ind), warning=False)
      except TypeError:
          # This error means no warning argument is available, so redirect stdout to a null device
          # to suppress printing of warning messages
          original_stdout = sys.stdout  # keep a reference to STDOUT
          sys.stdout = NullDevice()  # redirect the real STDOUT
          p, success = leastsq(errorfunction, p_tofit, args=(x, y, p_tofix, ind, image, err, g_ind))
          sys.stdout = original_stdout  # turn STDOUT back on
    else:
      p, sucess = None, 1

    para = N.zeros(6*ngaus)
    para[N.where(ind==1)[0]] = p
    para[N.where(ind==0)[0]] = p_tofix

    for igaus in range(ngaus):
      para[igaus*6+3] = abs(para[igaus*6+3])
      para[igaus*6+4] = abs(para[igaus*6+4])

    return para, success

def gaussian_fcn(g, x1, x2):
    """Evaluate Gaussian on the given grid.

    Parameters:
    x1, x2: grid (as produced by numpy.mgrid f.e.)
    g: Gaussian object or list of Gaussian paramters
    """
    from math import radians, sin, cos
    from const import fwsig
    import numpy as N

    if isinstance(g, list):
        A, C1, C2, S1, S2, Th = g
    else:
        A = g.peak_flux
        C1, C2 = g.centre_pix
        S1, S2, Th = g.size_pix
    S1 = S1/fwsig; S2 = S2/fwsig; Th = Th + 90.0 # Define theta = 0 on x-axis

    th = radians(Th)
    cs = cos(th)
    sn = sin(th)

    f1 = ((x1-C1)*cs + (x2-C2)*sn)/S1
    f2 = (-(x1-C1)*sn + (x2-C2)*cs)/S2

    return A*N.exp(-(f1*f1 + f2*f2)/2)

def mclean(im1, c, beam):
    """ Simple image plane clean of one gaussian at posn c and size=beam """
    import numpy as N

    amp = im1[c]
    b1, b2, b3 = beam
    b3 += 90.0
    para = [amp, c[0], c[1], b1, b2, b3]
    x, y = N.indices(im1.shape)

    im = gaus_2d(para, x, y)
    im1 = im1-im

    return im1

def arrstatmask(im, mask):
    """ Basic statistics for a masked array. dont wanna use numpy.ma """
    import numpy as N

    ind = N.where(~mask)
    im1 = im[ind]
    av = N.mean(im1)
    std = N.std(im1)
    maxv = N.max(im1)
    x, y = N.where(im == maxv)
    xmax = x[0]; ymax = y[0]
    minv = N.min(im1)
    x, y = N.where(im == minv)
    xmin = x[0]; ymin = y[0]

    return (av, std, maxv, (xmax, ymax), minv, (xmin, ymin))

def get_maxima(im, mask, thr, shape, beam):
    """ Gets the peaks in an image """
    from copy import deepcopy as cp
    import numpy as N

    im1 = cp(im)
    ind = N.array(N.where(~mask)).transpose()
    ind = [tuple(coord) for coord in ind if im[tuple(coord)] > thr]
    n, m = shape; iniposn = []; inipeak = []
    for c in ind:
      goodlist = [im[i,j] for i in range(c[0]-1,c[0]+2) for j in range(c[1]-1,c[1]+2) \
                   if i>=0 and i<n and j>=0 and j<m and (i,j) != c]
      peak = N.sum(im[c] > goodlist) == len(goodlist)
      if peak:
        iniposn.append(c); inipeak.append(im[c])
        im1 = mclean(im1, c, beam)

    return inipeak, iniposn, im1

def watershed(image, mask=None, markers=None, beam=None, thr=None):
      import numpy as N
      from copy import deepcopy as cp
      import scipy.ndimage as nd
      #import matplotlib.pyplot as pl
      #import pylab as pl

      if thr==None: thr = -1e9
      if mask==None: mask = N.zeros(image.shape, bool)
      if beam==None: beam = (2.0, 2.0, 0.0)
      if markers==None:
        inipeak, iniposn, im1 = get_maxima(image, mask, thr, image.shape, beam)
        ng = len(iniposn); markers = N.zeros(image.shape, int)
        for i in range(ng): markers[iniposn[i]] = i+2
        markers[N.unravel_index(N.argmin(image), image.shape)] = 1

      im1 = cp(image)
      if im1.min() < 0.: im1 = im1-im1.min()
      im1 = 255 - im1/im1.max()*255
      opw = nd.watershed_ift(N.array(im1, N.uint8), markers)

      return opw, markers

def get_kwargs(kwargs, key, typ, default):

    obj = True
    if kwargs.has_key(key):
      obj = kwargs[key]
    if not isinstance(obj, typ):
      obj = default

    return obj

def read_image_from_file(filename, img, indir, quiet=False):
    """ Reads data and header from indir/filename.

    We can use either pyfits or pyrap depending on the value
    of img.use_io = 'fits'/'rap'

    PyFITS is required, as it is used to standardize the header format. pyrap
    is optional.
    """
    import mylogger
    import os
    import numpy as N

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Readfile")
    if indir == None or indir == './':
        prefix = ''
    else:
        prefix = indir + '/'
    image_file = prefix + filename

    # Check that file exists
    if not os.path.exists(image_file):
        img._reason = 'File does not exist'
        return None

    # If img.use_io is set, then use appropriate io module
    if img.use_io != '':
        if img.use_io == 'fits':
            import pyfits
            try:
                fits = pyfits.open(image_file, mode="readonly", ignore_missing_end=True)
            except IOError, err:
                img._reason = 'Problem reading file.\nOriginal error: {0}'.format(str(err))
                return None
        if img.use_io == 'rap':
            import pyrap.images as pim
            try:
                inputimage = pim.image(image_file)
            except IOError, err:
                img._reason = 'Problem reading file.\nOriginal error: {0}'.format(str(err))
                return None
    else:
        # Simple check of whether pyrap and pyfits are available
        # We need pyfits version 2.2 or greater to use the
        # "ignore_missing_end" argument to pyfits.open().
        try:
            from distutils.version import StrictVersion
            import pyfits
            if StrictVersion(pyfits.__version__) > StrictVersion('2.2'):
                has_pyfits = True
            else:
                raise RuntimeError("PyFITS (version 2.2 or greater) is required.")
        except ImportError, err:
            raise RuntimeError("PyFITS is required.")
        try:
            import pyrap.images as pim
            has_pyrap = True
        except ImportError, err:
            has_pyrap = False
            e_pyrap = str(err)

        # First assume image is a fits file, and use pyfits to open it (if
        # available). If that fails, try to use pyrap if available.
        failed_read = False
        reason = 0
        try:
            fits = pyfits.open(image_file, mode="readonly", ignore_missing_end=True)
            img.use_io = 'fits'
        except IOError, err:
            e_pyfits = str(err)
            if has_pyrap:
                try:
                    inputimage = pim.image(image_file)
                    img.use_io = 'rap'
                except IOError, err:
                    e_pyrap = str(err)
                    failed_read = True
                    img._reason = 'File is not a valid FITS, CASA, or HDF5 image.'
            else:
                failed_read = True
                e_pyrap = "Pyrap unavailable"
                img._reason = 'Problem reading file.'
        if failed_read:
            img._reason += '\nOriginal errors: {0}\n {1}'.format(e_pyfits, e_pyrap)
            return None

    # Now that image has been read in successfully, get data and header
    if not quiet:
        mylogger.userinfo(mylog, "Opened '"+image_file+"'")
    if img.use_io == 'rap':
        data = inputimage.getdata()
        hdr = convert_pyrap_header(inputimage)
    if img.use_io == 'fits':
        data = fits[0].data
        hdr = fits[0].header
        fits.close()

    # Make sure data is in proper order. Final order is [pol, chan, x (RA), y (DEC)],
    # so we need to rearrange dimensions if they are not in this order. Use the
    # ctype FITS keywords to determine order of dimensions.
    mylog.info("Original data shape of " + image_file +': ' +str(data.shape))
    ctype_in = []
    for i in range(len(data.shape)):
        key_val_raw = hdr['CTYPE' + str(i+1)]
        key_val = key_val_raw.split('-')[0]
        ctype_in.append(key_val.strip())
    ctype_in.reverse() # Need to reverse order, as pyfits does this

    if 'RA' not in ctype_in or 'DEC' not in ctype_in:
        if 'GLON' not in ctype_in or 'GLAT' not in ctype_in:
            raise RuntimeError("Image data not found")
        else:
            lat_lon = True
    else:
        lat_lon = False
    if len(ctype_in) > 2 and 'FREQ' not in ctype_in:
        from pywcs import WCS
        t = WCS(hdr)
        t.wcs.fix()
        spec_indx = t.wcs.spec
        if spec_indx != -1:
            ctype_in.reverse()
            ctype_in[spec_indx] = 'FREQ'
            ctype_in.reverse()

    if lat_lon:
        ctype_out = ['STOKES', 'FREQ', 'GLON', 'GLAT']
    else:
        ctype_out = ['STOKES', 'FREQ', 'RA', 'DEC']
    indx_out = [-1, -1, -1, -1]
    indx_in = range(len(data.shape))
    for i in indx_in:
        for j in range(4):
            if ctype_in[i] == ctype_out[j]:
                indx_out[j] = i
    shape_out = [1, 1, data.shape[indx_out[2]], data.shape[indx_out[3]]]
    if indx_out[0] != -1:
        shape_out[0] = data.shape[indx_out[0]]
    if indx_out[1] != -1:
        shape_out[1] = data.shape[indx_out[1]]

    ### now we need to transpose columns to get the right order
    axes = range(len(data.shape))
    indx_out.reverse()
    for indx in indx_out:
        if indx != -1:
            axes.remove(indx)
            axes.insert(0, indx)
    data = data.transpose(*axes)
    data.shape = data.shape[0:4] # trim unused dimensions (if any)
    data = data.reshape(shape_out)
    mylog.info("Final data shape (npol, nchan, x, y): " + str(data.shape))

    ### and make a copy of it to get proper layout & byteorder
    data = N.array(data, order='C',
                   dtype=data.dtype.newbyteorder('='))

    ### trim image if trim_box is specified
    if img.opts.trim_box != None:
        img.trim_box = img.opts.trim_box
        xmin, xmax, ymin, ymax = img.trim_box
        if xmin < 0: xmin = 0
        if ymin < 0: ymin = 0
        if xmax > data.shape[2]: xmax = data.shape[2]
        if ymax > data.shape[3]: ymax = data.shape[3]
        if xmin >= xmax or ymin >= ymax:
            raise RuntimeError("The trim_box option does not specify a valid part of the image.")
        data = data[:, :, xmin:xmax, ymin:ymax]

        # Adjust WCS keywords
        hdr['crpix1'] -= xmin
        hdr['crpix2'] -= ymin
    else:
        img.trim_box = None

    return data, hdr


def convert_pyrap_header(pyrap_image):
    """Converts a pyrap header to a PyFITS header."""
    import tempfile
    import pyfits

    tfile = tempfile.NamedTemporaryFile(delete=False)
    pyrap_image.tofits(tfile.name)
    hdr = pyfits.getheader(tfile.name)
    return hdr


def write_image_to_file(use, filename, image, img, outdir=None,
                                           clobber=True):
    """ Writes image array to dir/filename using pyfits"""
    import numpy as N
    import os
    import mylogger

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Writefile")

    if filename == 'SAMP':
        import tempfile
        if not hasattr(img,'samp_client'):
            s, private_key = start_samp_proxy()
            img.samp_client = s
            img.samp_key = private_key

        # Broadcast image to SAMP Hub
        temp_im = make_fits_image(N.transpose(image), img.wcs_obj, img.beam, img.frequency)
        tfile = tempfile.NamedTemporaryFile(delete=False)
        temp_im.writeto(tfile.name,  clobber=clobber)
        send_fits_image(img.samp_client, img.samp_key, 'PyBDSM image', tfile.name)
    else:
        # Write image to FITS file
        import pyfits
        if outdir == None:
            outdir = img.indir
        if not os.path.exists(outdir) and outdir != '':
            os.makedirs(outdir)
        if os.path.exists(outdir + filename):
            if clobber:
                os.remove(outdir + filename)
            else:
                return
        temp_im = make_fits_image(N.transpose(image), img.wcs_obj, img.beam, img.frequency)
        temp_im.writeto(outdir + filename,  clobber=clobber)

def make_fits_image(imagedata, wcsobj, beam, freq):
    """Makes a simple FITS hdulist appropriate for single-channel images"""
    import pyfits
    shape_out = [1, 1, imagedata.shape[0], imagedata.shape[1]]
    hdu = pyfits.PrimaryHDU(imagedata.reshape(shape_out))
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header

    # Add WCS info
#     wcs_header = wcsobj.to_header()
#     for key in wcs_header.keys():
#         header.update(key, wcs_header[key])
    header.update('CRVAL1', wcsobj.wcs.crval[0])
    header.update('CDELT1', wcsobj.wcs.cdelt[0])
    header.update('CRPIX1', wcsobj.wcs.crpix[0])
    header.update('CUNIT1', wcsobj.wcs.cunit[0])
    header.update('CTYPE1', wcsobj.wcs.ctype[0])
    header.update('CRVAL2', wcsobj.wcs.crval[1])
    header.update('CDELT2', wcsobj.wcs.cdelt[1])
    header.update('CRPIX2', wcsobj.wcs.crpix[1])
    header.update('CUNIT2', wcsobj.wcs.cunit[1])
    header.update('CTYPE2', wcsobj.wcs.ctype[1])

    # Add STOKES info
    header.update('CRVAL3', 1)
    header.update('CDELT3', 1)
    header.update('CRPIX3', 1)
    header.update('CUNIT3', '')
    header.update('CTYPE3', 'STOKES')

    # Add or alter frequency info if needed
    header.update('CRVAL4', freq)
    header.update('CDELT4', 0.0)
    header.update('CRPIX4', 1)
    header.update('CUNIT4', 'Hz')
    header.update('CTYPE4', 'FREQ')

    # Add beam info
    header.update('BMAJ', beam[0])
    header.update('BMIN', beam[1])
    header.update('BPA', beam[2])


    # Add STOKES info
    header.update('CRVAL3', 1)
    header.update('CDELT3', 1)
    header.update('CRPIX3', 1)
    header.update('CUNIT3', '')
    header.update('CTYPE3', 'STOKES')

    # Add or alter frequency info if needed
    if wcsobj.wcs.spec != -1:
        header.update('CRVAL' + str(wcsobj.wcs.spec + 1), freq)
    else:
        header.update('CRVAL4', freq)
        header.update('CDELT4', 0.0)
        header.update('CRPIX4', 1)
        header.update('CUNIT4', 'Hz')
        header.update('CTYPE4', 'FREQ')

    # Add beam info
    header.update('BMAJ', beam[0])
    header.update('BMIN', beam[1])
    header.update('BPA', beam[2])

    hdulist[0].header = header
    return hdulist

def connect(mask):
    """ Find if a mask is singly or multiply connected """

    import scipy.ndimage as nd

    connectivity = nd.generate_binary_structure(2,2)
    labels, count = nd.label(mask, connectivity)
    if count > 1 :
      connected = 'multiple'
    else:
      connected = 'single'

    return connected, count

def area_polygon(points):
    """ Given an ANGLE ORDERED array points of [[x], [y]], find the total area by summing each successsive
    triangle with the centre """
    import numpy as N

    x, y = points
    n_tri = len(x)-1
    cenx, ceny = N.mean(x), N.mean(y)

    area = 0.0
    for i in range(n_tri):
      p1, p2, p3 = N.array([cenx, ceny]), N.array([x[i], y[i]]), N.array([x[i+1], y[i+1]])
      t_area= N.linalg.norm(N.cross((p2 - p1), (p3 - p1)))/2.
      area += t_area

    return area

def convexhull_deficiency(isl):
    """ Finds the convex hull for the island and returns the deficiency.
    Code taken from http://code.google.com/p/milo-lab/source/browse/trunk/src/toolbox/convexhull.py?spec=svn140&r=140
    """

    import random
    import time
    import numpy as N
    import scipy.ndimage as nd

    def _angle_to_point(point, centre):
        """calculate angle in 2-D between points and x axis"""
        delta = point - centre
        if delta[0] == 0.0:
            res = N.pi/2.0
        else:
            res = N.arctan(delta[1] / delta[0])
        if delta[0] < 0:
            res += N.pi
        return res

    def area_of_triangle(p1, p2, p3):
        """calculate area of any triangle given co-ordinates of the corners"""
        return N.linalg.norm(N.cross((p2 - p1), (p3 - p1)))/2.

    def convex_hull(points):
        """Calculate subset of points that make a convex hull around points
        Recursively eliminates points that lie inside two neighbouring points until only convex hull is remaining.
        points : ndarray (2 x m) array of points for which to find hull
        Returns: hull_points : ndarray (2 x n), convex hull surrounding points """

        n_pts = points.shape[1]
        #assert(n_pts > 5)
        centre = points.mean(1)
        angles = N.apply_along_axis(_angle_to_point, 0, points, centre)
        pts_ord = points[:,angles.argsort()]
        pts = [x[0] for x in zip(pts_ord.transpose())]
        prev_pts = len(pts) + 1
        k = 0
        while prev_pts > n_pts:
            prev_pts = n_pts
            n_pts = len(pts)
            i = -2
            while i < (n_pts - 2):
                Aij = area_of_triangle(centre, pts[i], pts[(i + 1) % n_pts])
                Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], \
                                       pts[(i + 2) % n_pts])
                Aik = area_of_triangle(centre, pts[i], pts[(i + 2) % n_pts])
                if Aij + Ajk < Aik:
                    del pts[i+1]
                i += 1
                n_pts = len(pts)
            k += 1
        return N.asarray(pts)

    mask = ~isl.mask_active
    points = N.asarray(N.where(mask - nd.binary_erosion(mask)))
    hull_pts = list(convex_hull(points))   # these are already in angle-sorted order

    hull_pts.append(hull_pts[0])
    hull_pts = N.transpose(hull_pts)

    isl_area = isl.size_active
    hull_area = area_polygon(hull_pts)
    ratio1 = hull_area/(isl_area - 0.5*len(hull_pts[0]))

    return ratio1


def open_isl(mask, index):
    """ Do an opening on a mask, divide left over pixels among opened sub islands. Mask = True => masked pixel """
    import scipy.ndimage as nd
    import numpy as N

    connectivity = nd.generate_binary_structure(2,2)
    ft = N.ones((index,index), int)

    open = nd.binary_opening(~mask, ft)
    open = check_1pixcontacts(open)  # check if by removing one pixel from labels, you can split a sub-island
    labels, n_subisl = nd.label(open, connectivity)  # get label/rank image for open. label = 0 for masked pixels
    labels, mask = assign_leftovers(mask, open, n_subisl, labels)  # add the leftover pixels to some island

    if labels != None:
        isl_pixs = [len(N.where(labels==i)[0]) for i in range(1,n_subisl+1)]
        isl_pixs = N.array(isl_pixs)/float(N.sum(isl_pixs))
    else:
        isl_pixs = None

    return n_subisl, labels, isl_pixs

def check_1pixcontacts(open):
    import scipy.ndimage as nd
    import numpy as N
    from copy import deepcopy as cp

    connectivity = nd.generate_binary_structure(2,2)
    ind = N.transpose(N.where(open[1:-1,1:-1] > 0)) + [1,1]   # exclude boundary to make it easier
    for pixel in ind:
      x, y = pixel
      grid = cp(open[x-1:x+2, y-1:y+2]); grid[1,1] = 0
      grid = N.where(grid == open[tuple(pixel)], 1, 0)
      ll, nn = nd.label(grid, connectivity)
      if nn > 1:
        open[tuple(pixel)] = 0

    return open

def assign_leftovers(mask, open, nisl, labels):
    """
    Given isl and the image of the mask after opening (open) and the number of new independent islands n,
    connect up the left over pixels to the new islands if they connect to only one island and not more.
    Assign the remaining to an island. We need to assign the leftout pixels to either of many sub islands.
    Easiest is to assign to the sub island with least size.
    """
    import scipy.ndimage as nd
    import numpy as N
    from copy import deepcopy as cp

    n, m = mask.shape
    leftout = ~mask - open

    connectivity = nd.generate_binary_structure(2,2)
    mlabels, count = nd.label(leftout, connectivity)
    npix = [len(N.where(labels==b)[0]) for b in range(1,nisl+1)]

    for i_subisl in range(count):
      c_list = []    # is list of all bordering pixels of the sub island
      ii = i_subisl+1
      coords = N.transpose(N.where(mlabels==ii))  # the coordinates of island i of left-out pixels
      for co in coords:
        co8 = [[x,y] for x in range(co[0]-1,co[0]+2) for y in range(co[1]-1,co[1]+2) if x >=0 and y >=0 and x <n and y<m]
        c_list.extend([tuple(cc) for cc in co8 if mlabels[tuple(cc)] == 0])
      c_list = list(set(c_list))     # to avoid duplicates
      vals = N.array([labels[c] for c in c_list])
      belongs = list(set(vals[N.nonzero(vals)]))
      if len(belongs) == 0:
        # No suitable islands found => mask pixels
        for cc in coords:
            mask = (mlabels == ii)
#             mask[cc] = True
            return None, mask
      if len(belongs) == 1:
        for cc in coords:
          labels[tuple(cc)] = belongs[0]
      else:                             # get the border pixels of the islands
        nn = [npix[b-1] for b in belongs]
        addto = belongs[N.argmin(nn)]
        for cc in coords:
          labels[tuple(cc)] = addto

    return labels, mask


def _float_approx_equal(x, y, tol=1e-18, rel=1e-7):
    if tol is rel is None:
        raise TypeError('cannot specify both absolute and relative errors are None')
    tests = []
    if tol is not None: tests.append(tol)
    if rel is not None: tests.append(rel*abs(x))
    assert tests
    return abs(x - y) <= max(tests)


def approx_equal(x, y, *args, **kwargs):
    """approx_equal(float1, float2[, tol=1e-18, rel=1e-7]) -> True|False
    approx_equal(obj1, obj2[, *args, **kwargs]) -> True|False

    Return True if x and y are approximately equal, otherwise False.

    If x and y are floats, return True if y is within either absolute error
    tol or relative error rel of x. You can disable either the absolute or
    relative check by passing None as tol or rel (but not both).

    For any other objects, x and y are checked in that order for a method
    __approx_equal__, and the result of that is returned as a bool. Any
    optional arguments are passed to the __approx_equal__ method.

    __approx_equal__ can return NotImplemented to signal that it doesn't know
    how to perform that specific comparison, in which case the other object is
    checked instead. If neither object have the method, or both defer by
    returning NotImplemented, approx_equal falls back on the same numeric
    comparison used for floats.

    >>> almost_equal(1.2345678, 1.2345677)
    True
    >>> almost_equal(1.234, 1.235)
    False

    """
    if not (type(x) is type(y) is float):
        # Skip checking for __approx_equal__ in the common case of two floats.
        methodname = '__approx_equal__'
        # Allow the objects to specify what they consider "approximately equal",
        # giving precedence to x. If either object has the appropriate method, we
        # pass on any optional arguments untouched.
        for a,b in ((x, y), (y, x)):
            try:
                method = getattr(a, methodname)
            except AttributeError:
                continue
            else:
                result = method(b, *args, **kwargs)
                if result is NotImplemented:
                    continue
                return bool(result)
    # If we get here without returning, then neither x nor y knows how to do an
    # approximate equal comparison (or are both floats). Fall back to a numeric
    # comparison.
    return _float_approx_equal(x, y, *args, **kwargs)

def isl_tosplit(isl, opts):
    """ Splits an island and sends back parameters """
    import numpy as N

    size_extra5 = opts.splitisl_size_extra5
    frac_bigisl3 = opts.splitisl_frac_bigisl3

    connected, count = connect(isl.mask_active)
    index = 0
    n_subisl3, labels3, isl_pixs3 = open_isl(isl.mask_active, 3)
    n_subisl5, labels5, isl_pixs5 = open_isl(isl.mask_active, 5)
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
    convex_def =  convexhull_deficiency(isl)
    #print 'CONVEX = ',convex_def

    if opts.plot_islands:
        try:
            import matplotlib.pyplot as pl
            pl.figure()
            pl.suptitle('Island '+str(isl.island_id))
            pl.subplot(2,2,1); pl.imshow(N.transpose(isl.image*~isl.mask_active), origin='lower', interpolation='nearest'); pl.title('Image')
            pl.subplot(2,2,2); pl.imshow(N.transpose(labels3), origin='lower', interpolation='nearest'); pl.title('labels3')
            pl.subplot(2,2,3); pl.imshow(N.transpose(labels5), origin='lower', interpolation='nearest'); pl.title('labels5')
        except ImportError:
            print "\033[31;1mWARNING\033[0m: Matplotlib not found. Plotting disabled."
    if index == 0: return [index, n_subisl5, labels5]
    else: return [index, n_subisl, labels]


class NullDevice():
    """Null device to suppress stdout, etc."""
    def write(self, s):
        pass

def ch0_aperture_flux(img, posn_pix, aperture_pix):
    """Measure ch0 flux inside radius aperture_pix pixels centered on posn_pix.

    Returns [flux, fluxE]
    """
    import numpy as N

    if aperture_pix == None:
        return [0.0, 0.0]

    # Make ch0 and rms subimages
    xlo = posn_pix[0]-int(aperture_pix)-1
    if xlo < 0:
        xlo = 0
    xhi = posn_pix[0]+int(aperture_pix)+1
    if xhi > img.ch0.shape[0]:
        xhi = img.ch0.shape[0]
    ylo = posn_pix[1]-int(aperture_pix)-1
    if ylo < 0:
        ylo = 0
    yhi = posn_pix[1]+int(aperture_pix)+1
    if yhi > img.ch0.shape[1]:
        yhi = img.ch0.shape[1]

    aper_im = img.ch0[xlo:xhi, ylo:yhi]
    aper_rms = img.rms[xlo:xhi, ylo:yhi]
    posn_pix_new = [posn_pix[0]-xlo, posn_pix[1]-ylo]
    aper_flux = aperture_flux(aperture_pix, posn_pix_new, aper_im, aper_rms, img.pixel_beamarea)
    return aper_flux

def aperture_flux(aperture_pix, posn_pix, aper_im, aper_rms, beamarea):
    """Returns aperture flux and error"""
    import numpy as N

    dist_mask = generate_aperture(aper_im.shape[0], aper_im.shape[1], posn_pix[1], posn_pix[0], aperture_pix)
    aper_mask = N.where(dist_mask)
    if N.size(aper_mask) == 0:
        return [0.0, 0.0]
    aper_flux = N.nansum(aper_im[aper_mask])/beamarea # Jy
    pixels_in_source = N.sum(~N.isnan(aper_im[aper_mask])) # number of unmasked pixels assigned to current source
    aper_fluxE = nanmean(aper_rms[aper_mask]) * N.sqrt(pixels_in_source/beamarea) # Jy
    return [aper_flux, aper_fluxE]

def generate_aperture(ysize, xsize, ycenter, xcenter, radius):
    """Makes a mask for a circular aperture"""
    import numpy

    x, y = numpy.mgrid[0:ysize,0:xsize]
    return ((x - ycenter)**2 + (y - xcenter)**2 <= radius**2) * 1

def getTerminalSize():
    """
    returns (lines:int, cols:int)
    """
    import os, struct
    def ioctl_GWINSZ(fd):
        import fcntl, termios
        return struct.unpack("hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "1234"))
    # try stdin, stdout, stderr
    for fd in (0, 1, 2):
        try:
            return ioctl_GWINSZ(fd)
        except:
            pass
    # try os.ctermid()
    try:
        fd = os.open(os.ctermid(), os.O_RDONLY)
        try:
            return ioctl_GWINSZ(fd)
        finally:
            os.close(fd)
    except:
        pass
    # try `stty size`
    try:
        return tuple(int(x) for x in os.popen("stty size", "r").read().split())
    except:
        pass
    # try environment variables
    try:
        return tuple(int(os.getenv(var)) for var in ("LINES", "COLUMNS"))
    except:
        pass
    # Give up. return 0.
    return (0, 0)

def eval_func_tuple(f_args):
    """Takes a tuple of a function and args, evaluates and returns result

    This function (in addition to itertools) gets around limitation that
    multiple-argument sequences are not supported by multiprocessing.
    """
    return f_args[0](*f_args[1:])


def start_samp_proxy():
    """Starts (registers) and returns a SAMP proxy"""
    import os
    import xmlrpclib

    lockfile = os.path.expanduser('~/.samp')
    if not os.path.exists(lockfile):
        raise RuntimeError("A running SAMP hub was not found.")
    else:
        HUB_PARAMS = {}
        for line in open(lockfile):
            if not line.startswith('#'):
                key, value = line.split('=', 1)
                HUB_PARAMS[key] = value.strip()

    # Set up proxy
    s = xmlrpclib.ServerProxy(HUB_PARAMS['samp.hub.xmlrpc.url'])

    # Register with Hub
    metadata = {"samp.name": 'PyBDSM', "samp.description.text": 'PyBDSM: the Python Blob Detection and Source Measurement software'}
    result = s.samp.hub.register(HUB_PARAMS['samp.secret'])
    private_key = result['samp.private-key']
    s.samp.hub.declareMetadata(private_key, metadata)
    return s, private_key


def stop_samp_proxy(img):
    """Stops (unregisters) a SAMP proxy"""
    import os

    if hasattr(img, 'samp_client'):
        lockfile = os.path.expanduser('~/.samp')
        if os.path.exists(lockfile):
            img.samp_client.samp.hub.unregister(img.samp_key)


def send_fits_image(s, private_key, name, file_path):
    """Send a SAMP notification to load a fits image."""
    import os

    message = {}
    message['samp.mtype'] = "image.load.fits"
    message['samp.params'] = {}
    message['samp.params']['url'] = 'file://' + os.path.abspath(file_path)
    message['samp.params']['name'] = name
    lockfile = os.path.expanduser('~/.samp')
    if not os.path.exists(lockfile):
        raise RuntimeError("A running SAMP hub was not found.")
    else:
        s.samp.hub.notifyAll(private_key, message)

def send_fits_table(s, private_key, name, file_path):
    """Send a SAMP notification to load a fits table."""
    import os

    message = {}
    message['samp.mtype'] = "table.load.fits"
    message['samp.params'] = {}
    message['samp.params']['url'] = 'file://' + os.path.abspath(file_path)
    message['samp.params']['name'] = name
    lockfile = os.path.expanduser('~/.samp')
    if not os.path.exists(lockfile):
        raise RuntimeError("A running SAMP hub was not found.")
    else:
        s.samp.hub.notifyAll(private_key, message)

def send_highlight_row(s, private_key, url, row_id):
    """Send a SAMP notification to highlight a row in a table."""
    import os

    message = {}
    message['samp.mtype'] = "table.highlight.row"
    message['samp.params'] = {}
    message['samp.params']['row'] = str(row_id)
    message['samp.params']['url'] = url
    lockfile = os.path.expanduser('~/.samp')
    if not os.path.exists(lockfile):
        raise RuntimeError("A running SAMP hub was not found.")
    else:
        s.samp.hub.notifyAll(private_key, message)

def send_coords(s, private_key, coords):
    """Send a SAMP notification to point at given coordinates."""
    import os

    message = {}
    message['samp.mtype'] = "coord.pointAt.sky"
    message['samp.params'] = {}
    message['samp.params']['ra'] = str(coords[0])
    message['samp.params']['dec'] = str(coords[1])
    lockfile = os.path.expanduser('~/.samp')
    if not os.path.exists(lockfile):
        raise RuntimeError("A running SAMP hub was not found.")
    else:
        s.samp.hub.notifyAll(private_key, message)
