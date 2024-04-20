# some functions
from __future__ import print_function
from __future__ import absolute_import

try:
    # For Python 2
    basestring = basestring
except NameError:
    basestring = str

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
            print('Not yet implemented')

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
        print(" We do not trust polynomial fits > 3 ")
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
    from .const import fwsig, pi
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
        print(" Something wrong with the parameters passed - need multiples of 6 !")
    else:
        ngaus = int(len(indx)/6)
        params = N.zeros(6*ngaus)
        params[N.where(indx==1)[0]] = c
        params[N.where(indx==0)[0]] = p_tofix
        for i in range(ngaus):
            gau = params[i*6:i*6+6]
            val = val + gaus_2d(gau, x, y)

    return val

def g2param(g, adj=False):
    """Convert gaussian object g to param list [amp, cenx, ceny, sigx, sigy, theta] """
    from .const import fwsig
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
    from .const import fwsig
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

    from .const import fwsig

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
    import numpy as N
    from .gausfit import Gaussian

    rad = 180.0/N.pi
    if isinstance(g, Gaussian):
        param = g2param(g)
    else:
        if isinstance(g, list) and len(g)>=6:
            param = g
        else:
            raise RuntimeError("Input to drawellipse neither Gaussian nor list")

    size = [param[3], param[4], param[5]]
    size_fwhm = corrected_size(size)
    th=N.arange(0, 370, 10)
    x1=size_fwhm[0]*N.cos(th/rad)
    y1=size_fwhm[1]*N.sin(th/rad)
    x2=x1*N.cos(param[5]/rad)-y1*N.sin(param[5]/rad)+param[1]
    y2=x1*N.sin(param[5]/rad)+y1*N.cos(param[5]/rad)+param[2]

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
    from .const import fwsig

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

    if mask is None:
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
    if N.all(m3/m1 > m2*m2):
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

        if p0 is None:
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
            if cov is not None:
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


def angsep(ra1, dec1, ra2, dec2):
    """Returns angular separation between two coordinates (all in degrees)"""
    import math

    const = math.pi/180.
    ra1 = ra1*const
    rb1 = dec1*const
    ra2 = ra2*const
    rb2 = dec2*const

    v1_1 = math.cos(ra1)*math.cos(rb1)
    v1_2 = math.sin(ra1)*math.cos(rb1)
    v1_3 = math.sin(rb1)

    v2_1 = math.cos(ra2)*math.cos(rb2)
    v2_2 = math.sin(ra2)*math.cos(rb2)
    v2_3 = math.sin(rb2)

    w = ( (v1_1-v2_1)**2 + (v1_2-v2_2)**2 + (v1_3-v2_3)**2 )/4.0

    x = math.sqrt(w)
    y = math.sqrt(max(0.0, 1.0-w))
    angle = 2.0*math.atan2(x, y)/const
    return angle


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
    import scipy.fft
    from scipy import ndimage

    shape=image.shape

    f1=scipy.fft.fft(image, shape[0], axis=0)
    f2=scipy.fft.fft(f1, shape[1], axis=1)

    s=ndimage.fourier_shift(f2,shift, axis=0)

    y1=scipy.fft.ifft(s, shape[1], axis=1)
    y2=scipy.fft.ifft(y1, shape[0], axis=0)

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
    from .const import fwsig
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

    if mask is not None and mask.shape != data.shape:
        print('Data and mask array dont have the same shape, ignoring mask')
        mask = None
    if err is not None and err.shape != data.shape:
        print('Data and error array dont have the same shape, ignoring error')
        err = None

    if mask is None: mask = N.zeros(data.shape, bool)
    g_ind = N.where(~N.ravel(mask))[0]

    if err is None:
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


def get_errors(img, p, stdav, bm_pix=None, fixed_to_beam=False):
    """ Returns errors on the fitted Gaussian parameters, using the
    equations from Condon 1997 (PASP, 109, 166) and Condon et al.
    1998 (ApJ, 115, 1693)

    Parameters:
    img: Image object (needed for pixel beam info)
    p: list of Gaussian parameters: [peak, x0, y0, maj, min, pa, tot]
    stdav: estimate of the image noise at the Gaussian's position
    bm_pix: optional pixel beam to be used instead of that in img
    fixed_to_beam: True if the fits were done with the
        size fixed to that of the beam, False otherwise

    Returned list includes errors on:
        peak flux [Jy/beam]
        x_0 [pix]
        y_0 [pix]
        e_maj [pix]
        e_min [pix]
        e_pa [deg]
        e_tot [Jy]

    """
    from .const import fwsig
    from math import sqrt, log, pow, pi
    from . import mylogger
    import numpy as N

    mylog = mylogger.logging.getLogger("PyBDSM.Compute")

    if len(p) % 7 > 0:
        mylog.error("Gaussian parameters passed have to have 7n numbers")
    ngaus = int(len(p)/7)
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
            if bm_pix is None:
                bm_pix = N.array([img.pixel_beam()[0]*fwsig, img.pixel_beam()[1]*fwsig, img.pixel_beam()[2]])
            dumr = sqrt(abs(size[0] * size[1] / (4.0 * bm_pix[0] * bm_pix[1])))  # first term of Eq. 26 of Condon+ (1998)
            dumrr1 = 1.0 + bm_pix[0] * bm_pix[1] / (size[0] * size[0])  # second term of Eq. 26 of Condon+ (1998)
            dumrr2 = 1.0 + bm_pix[0] * bm_pix[1] / (size[1] * size[1])  # third term of Eq. 26 of Condon+ (1998)
            dumrr3 = dumr * pp[0] / stdav  # product of first and fourth terms of Eq. 26 of Condon+ (1998)
            d1 = sqrt(8.0 * log(2.0))
            d2 = (size[0] * size[0] - size[1] * size[1]) / (size[0] * size[1])  # last term of Eq. 30 of Condon+ (1998)
            try:
                # The following three errors are calculated using Eq. 21 of Condon (1997),
                # using Eq. 26 of Condon+ (1998) for rho
                e_peak = pp[0] * sq2 / (dumrr3 * pow(dumrr1, 0.75) * pow(dumrr2, 0.75))
                e_maj = size[0] * sq2 / (dumrr3 * pow(dumrr1, 1.25) * pow(dumrr2, 0.25))
                e_min = size[1] * sq2 / (dumrr3 * pow(dumrr1, 0.25) * pow(dumrr2, 1.25))

                # The following two errors are calculated using Eq. 27 of Condon+ (1998)
                pa_rad = size[2] * pi / 180.0
                e_x0 = sqrt( (e_maj * N.sin(pa_rad))**2 + (e_min * N.cos(pa_rad))**2 ) / d1
                e_y0 = sqrt( (e_maj * N.cos(pa_rad))**2 + (e_min * N.sin(pa_rad))**2 ) / d1

                # The following error is calculated using Eq. 30 of Condon+ (1998)
                e_pa = 2.0 / (d2 * dumrr3 * pow(dumrr1, 0.25) * pow(dumrr2, 1.25))
                e_pa = e_pa * 180.0/pi

                # The following error is calculated using Eq. 36 of Condon+ (1998)
                e_tot = pp[6] * sqrt(e_peak * e_peak / (pp[0] * pp[0]) + (0.25 / dumr / dumr) *
                        (e_maj * e_maj / (size[0] * size[0]) + e_min * e_min / (size[1] * size[1])))
            except:
                e_peak = 0.0
                e_x0 = 0.0
                e_y0 = 0.0
                e_maj = 0.0
                e_min = 0.0
                e_pa = 0.0
                e_tot = 0.0
            if abs(e_pa) > 180.0:
                e_pa = 180.0
            if fixed_to_beam:
                # When the size was fixed to that of the beam during the fit, set
                # uncertainties on the size to zero and reduce the error in the fluxes
                # by sqrt(2) (see Eq. 25 of Condon 1997)
                e_maj = 0.0
                e_min = 0.0
                e_pa = 0.0
                e_peak /= sq2
                e_tot /= sq2
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

    if mask is not None and mask.shape != image.shape:
        print('Data and mask array dont have the same shape, ignoring mask')
        mask = None
    if err is not None and err.shape != image.shape:
        print('Data and error array dont have the same shape, ignoring error')
        err = None
    if mask is None: mask = N.zeros(image.shape, bool)

    g_ind = N.where(~N.ravel(mask))[0]

    ngaus = len(gaus)
    if ngaus > 0:
        p_ini = []
        for g in gaus:
            p_ini = p_ini + g2param(g, adj)
        p_ini = N.array(p_ini)

        if fitfix is None: fitfix = [0]*ngaus
        ind = N.ones(6*ngaus)                                     # 1 => fit ; 0 => fix
        for i in range(ngaus):
            if fitfix[i] == 1: ind[i*6+1:i*6+6] = 0
            if fitfix[i] == 2: ind[i*6+3:i*6+6] = 0
            if fitfix[i] == 3: ind[i*6+1:i*6+3] = 0
        ind = N.array(ind)
        p_tofit = p_ini[N.where(ind==1)[0]]
        p_tofix = p_ini[N.where(ind==0)[0]]
        if err is None: err = N.ones(image.shape)

        errorfunction = lambda p, x, y, p_tofix, ind, image, err, g_ind: \
                       N.ravel((gaus_2d_itscomplicated(p, x, y, p_tofix, ind)-image)/err)[g_ind]
        try:
            p, success = leastsq(errorfunction, p_tofit, args=(x, y, p_tofix, ind, image, err, g_ind))
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
    from .const import fwsig
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

def get_maxima(im, mask, thr, shape, beam, im_pos=None):
    """ Gets the peaks in an image """
    from copy import deepcopy as cp
    import numpy as N

    if im_pos is None:
        im_pos = im
    im1 = cp(im)
    ind = N.array(N.where(~mask)).transpose()
    ind = [tuple(coord) for coord in ind if im_pos[tuple(coord)] > thr]
    n, m = shape
    iniposn = []
    inipeak = []
    for c in ind:
        goodlist = [im_pos[i,j] for i in range(c[0]-1,c[0]+2) for j in range(c[1]-1,c[1]+2) \
                     if i>=0 and i<n and j>=0 and j<m and (i,j) != c]
        peak = N.sum(im_pos[c] > goodlist) == len(goodlist)
        if peak:
            iniposn.append(c)
            inipeak.append(im[c])
            im1 = mclean(im1, c, beam)

    return inipeak, iniposn, im1

def watershed(image, mask=None, markers=None, beam=None, thr=None):
    import numpy as N
    from copy import deepcopy as cp
    import scipy.ndimage as nd
    #import matplotlib.pyplot as pl
    #import pylab as pl

    if thr is None: thr = -1e9
    if mask is None: mask = N.zeros(image.shape, bool)
    if beam is None: beam = (2.0, 2.0, 0.0)
    if markers is None:
        inipeak, iniposn, im1 = get_maxima(image, mask, thr, image.shape, beam)
        ng = len(iniposn); markers = N.zeros(image.shape, int)
        for i in range(ng): markers[iniposn[i]] = i+2
        markers[N.unravel_index(N.argmin(image), image.shape)] = 1

    im1 = cp(image)
    if im1.min() < 0.: im1 = im1-im1.min()
    im1 = 255 - im1/im1.max()*255
    opw = nd.watershed_ift(N.array(im1, N.uint16), markers)

    return opw, markers

def get_kwargs(kwargs, key, typ, default):
    obj = True
    if key in kwargs:
        obj = kwargs[key]
    if not isinstance(obj, typ):
        obj = default

    return obj

def read_image_from_file(filename, img, indir, quiet=False):
    """ Reads data and header from indir/filename.

    We can use either pyfits or python-casacore depending on the value
    of img.use_io = 'fits'/'rap'

    PyFITS is required, as it is used to standardize the header format. python-casacore
    is optional.
    """
    from . import mylogger
    import os
    import numpy as N
    from astropy.io import fits as pyfits
    from astropy.wcs import WCS
    from copy import deepcopy as cp
    import warnings

    # Check if casacore is available, which is needed for 'rap' type of image I/O.
    try:
        import casacore.images as pim
        has_casacore = True
    except ImportError as err:
        has_casacore = False

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Readfile")
    if indir is None or indir == './':
        prefix = ''
    else:
        prefix = indir + '/'
    image_file = prefix + filename

    # Check that file exists
    if not os.path.exists(image_file):
        img._reason = f'File {image_file} does not exist'
        return None

    # If img.use_io is set, then use appropriate io module
    if img.use_io:
        # Sanity check: only 'fits' and 'rap' are supported image I/O types
        if img.use_io not in ('fits', 'rap'):
            raise ValueError(f"Invalid image I/O type '{img.use_io}'. "
                             "Supported types are: 'fits' and 'rap'")
        if img.use_io == 'fits':
            try:
                fits = pyfits.open(image_file, mode="readonly", ignore_missing_end=True)
            except IOError as err:
                img._reason = f'Problem reading {image_file}.\nOriginal error: {err}'
                return None
        if img.use_io == 'rap':
            if not has_casacore:
                img._reason = f'Problem reading {image_file}.\nCasacore is unavailable'
                return None
            try:
                inputimage = pim.image(image_file)
            except IOError as err:
                img._reason = f'Problem reading {image_file}.\nOriginal error: {err}'
                return None
    else:
        # First assume image is a fits file, and use pyfits to open it.
        # If that fails, try to use casacore if available.
        failed_read = False
        try:
            fits = pyfits.open(image_file, mode="readonly", ignore_missing_end=True)
            img.use_io = 'fits'
        except IOError as err:
            e_pyfits = str(err)
            if has_casacore:
                try:
                    inputimage = pim.image(image_file)
                    img.use_io = 'rap'
                except IOError as err:
                    e_casacore = str(err)
                    failed_read = True
                    img._reason = 'File is not a valid FITS, CASA, or HDF5 image.'
            else:
                failed_read = True
                e_casacore = "Casacore unavailable"
                img._reason = f'Problem reading {image_file}.'
        if failed_read:
            img._reason += f'\nOriginal error: {e_pyfits}\n {e_casacore}'
            return None

    # Now that image has been read in successfully, get header (data is loaded
    # later to take advantage of sectioning if trim_box is specified).
    if not quiet:
        mylogger.userinfo(mylog, "Opened '"+image_file+"'")
    if img.use_io == 'rap':
        tmpdir = os.path.join(img.outdir, img.parentname+'_tmp')
        hdr = convert_casacore_header(inputimage, tmpdir)
        coords = inputimage.coordinates()
        img.coords_dict = coords.dict()
        if 'telescope' in img.coords_dict:
            img._telescope = img.coords_dict['telescope']
        else:
            img._telescope = None
    if img.use_io == 'fits':
        hdr = fits[0].header
        img.coords_dict = None
        if 'TELESCOP' in hdr:
            img._telescope = hdr['TELESCOP']
        else:
            img._telescope = None

    # Make sure data is in proper order. Final order is [pol, chan, x (RA), y (DEC)],
    # so we need to rearrange dimensions if they are not in this order. Use the
    # ctype FITS keywords to determine order of dimensions. Note that both PyFITS
    # and casacore reverse the order of the axes relative to NAXIS, so we must too.
    naxis = hdr['NAXIS']
    data_shape = []
    for i in range(naxis):
        data_shape.append(hdr['NAXIS'+str(i+1)])
    data_shape.reverse()
    data_shape = tuple(data_shape)
    mylog.info("Original data shape of " + image_file +': ' +str(data_shape))
    ctype_in = []
    for i in range(naxis):
        key_val_raw = hdr['CTYPE' + str(i+1)]
        key_val = key_val_raw.split('-')[0]
        ctype_in.append(key_val.strip())
    if 'RA' not in ctype_in or 'DEC' not in ctype_in:
        if 'GLON' not in ctype_in or 'GLAT' not in ctype_in:
            raise RuntimeError("Image data not found")
        else:
            lat_lon = True
    else:
        lat_lon = False

    # Check for incorrect spectral units. For example, "M/S" is not
    # recognized by PyWCS as velocity ("S" is actually Siemens, not
    # seconds). Note that we check CUNIT3 and CUNIT4 even if the
    # image has only 2 axes, as the header may still have these
    # entries.
    for i in range(4):
        key_val_raw = hdr.get('CUNIT' + str(i+1))
        if key_val_raw is not None:
            if 'M/S' in key_val_raw or 'm/S' in key_val_raw or 'M/s' in key_val_raw:
                hdr['CUNIT' + str(i+1)] = 'm/s'
            if 'HZ' in key_val_raw or 'hZ' in key_val_raw or 'hz' in key_val_raw:
                hdr['CUNIT' + str(i+1)] = 'Hz'
            if 'DEG' in key_val_raw or 'Deg' in key_val_raw:
                hdr['CUNIT' + str(i+1)] = 'deg'

    # Make sure that the spectral axis has been identified properly
    if len(ctype_in) > 2 and 'FREQ' not in ctype_in:
        from astropy.wcs import FITSFixedWarning
        # TODO: Is it still needed and/or desirable to filter warnings?
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=DeprecationWarning)
            warnings.filterwarnings("ignore",category=FITSFixedWarning)
            t = WCS(hdr)
            t.wcs.fix()
        spec_indx = t.wcs.spec
        if spec_indx != -1:
            ctype_in[spec_indx] = 'FREQ'

    # Now reverse the axes order to match PyFITS/casacore order and define the
    # final desired order (cytpe_out) and shape (shape_out).
    ctype_in.reverse()
    if lat_lon:
        ctype_out = ['STOKES', 'FREQ', 'GLON', 'GLAT']
    else:
        ctype_out = ['STOKES', 'FREQ', 'RA', 'DEC']
    indx_out = [-1, -1, -1, -1]
    indx_in = range(naxis)
    for i in indx_in:
        for j in range(4):
            if ctype_in[i] == ctype_out[j]:
                indx_out[j] = i
    shape_out = [1, 1, data_shape[indx_out[2]], data_shape[indx_out[3]]]
    if indx_out[0] != -1:
        shape_out[0] = data_shape[indx_out[0]]
    if indx_out[1] != -1:
        shape_out[1] = data_shape[indx_out[1]]
    indx_out = [a for a in indx_out if a >= 0] # trim unused axes

    # Read in data. If only a subsection of the image is desired (as defined
    # by the trim_box option), we can try to use PyFITS to read only that section.
    img._original_naxis = data_shape
    img._original_shape = (shape_out[2], shape_out[3])
    img._xy_hdr_shift = (0, 0)
    if img.opts.trim_box is not None:
        img.trim_box = [int(b) for b in img.opts.trim_box]
        xmin, xmax, ymin, ymax = img.trim_box
        if xmin < 0: xmin = 0
        if ymin < 0: ymin = 0
        if xmax > shape_out[2]: xmax = shape_out[2]
        if ymax > shape_out[3]: ymax = shape_out[3]
        if xmin >= xmax or ymin >= ymax:
            raise RuntimeError("The trim_box option does not specify a valid part of the image.")
        shape_out_untrimmed = shape_out[:]
        shape_out[2] = xmax-xmin
        shape_out[3] = ymax-ymin

        if img.use_io == 'fits':
            sx = slice(int(xmin),int(xmax))
            sy = slice(int(ymin),int(ymax))
            sn = slice(None)
            s_array = [sx, sy]
            for i in range(naxis-2):
                s_array.append(sn)
            s_array.reverse() # to match ordering of data array returned by PyFITS
            if naxis == 2:
                data = fits[0].section[s_array[0], s_array[1]]
            elif naxis == 3:
                data = fits[0].section[s_array[0], s_array[1], s_array[2]]
            elif naxis == 4:
                data = fits[0].section[s_array[0], s_array[1], s_array[2], s_array[3]]
            else:
                # If more than 4 axes, just read in the whole image and
                # do the trimming after reordering.
                data = fits[0].data
            fits.close()
            data = data.transpose(*indx_out) # transpose axes to final order
            data.shape = data.shape[0:4] # trim unused dimensions (if any)
            if naxis > 4:
                data = data.reshape(shape_out_untrimmed) # Add axes if needed
                data = data[:, :, xmin:xmax, ymin:ymax] # trim to trim_box
            else:
                data = data.reshape(shape_out) # Add axes if needed
        else:
            # With casacore, just read in the whole image and then trim
            data = inputimage.getdata()
            data = data.transpose(*indx_out) # transpose axes to final order
            data.shape = data.shape[0:4] # trim unused dimensions (if any)
            data = data.reshape(shape_out_untrimmed) # Add axes if needed
            data = data[:, :, xmin:xmax, ymin:ymax] # trim to trim_box

        # Adjust WCS keywords for trim_box starting x and y.
        hdr['crpix1'] -= xmin
        hdr['crpix2'] -= ymin
        img._xy_hdr_shift = (xmin, ymin)
    else:
        if img.use_io == 'fits':
            data = fits[0].data
            fits.close()
        else:
            data = inputimage.getdata()
        data = data.transpose(*indx_out) # transpose axes to final order
        data.shape = data.shape[0:4] # trim unused dimensions (if any)
        data = data.reshape(shape_out) # Add axes if needed

    mylog.info("Final data shape (npol, nchan, x, y): " + str(data.shape))

    return data, hdr


def convert_casacore_header(casacore_image, tmpdir):
    """Converts a casacore header to a PyFITS header."""
    import tempfile
    import os
    import atexit
    import shutil
    try:
        from astropy.io import fits as pyfits
    except ImportError as err:
        import pyfits

    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    tfile = tempfile.NamedTemporaryFile(delete=False, dir=tmpdir)
    casacore_image.tofits(tfile.name)
    hdr = pyfits.getheader(tfile.name)
    if os.path.isfile(tfile.name):
        os.remove(tfile.name)

    # Register deletion of temp directory at exit to be sure it is deleted
    atexit.register(shutil.rmtree, tmpdir, ignore_errors=True)

    return hdr


def write_image_to_file(use, filename, image, img, outdir=None,
                        pad_image=False, clobber=True, is_mask=False):
    """ Writes image array to outdir/filename"""
    import numpy as N
    import os
    from . import mylogger

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Writefile")

    wcs_obj = img.wcs_obj
    if pad_image and img.opts.trim_box is not None:
        # Pad image to original size
        xsize, ysize = img._original_shape
        xmin, ymin = img._xy_hdr_shift
        image_pad = N.zeros((xsize, ysize), dtype=N.float32)
        image_pad[xmin:xmin+image.shape[0], ymin:ymin+image.shape[1]] = image
        image = image_pad
    else:
        xmin = 0
        ymin = 0

    if not hasattr(img, '_telescope'):
        telescope = None
    else:
        telescope = img._telescope

    if filename == 'SAMP':
        import tempfile
        if not hasattr(img,'samp_client'):
            s, private_key = start_samp_proxy()
            img.samp_client = s
            img.samp_key = private_key

        # Broadcast image to SAMP Hub
        temp_im = make_fits_image(N.transpose(image), wcs_obj, img.beam,
            img.frequency, img.equinox, telescope, xmin=xmin, ymin=ymin,
            is_mask=is_mask)
        tfile = tempfile.NamedTemporaryFile(delete=False)
        try:
            temp_im.writeto(tfile.name, overwrite=clobber)
        except TypeError:
            # The "overwrite" argument was added in astropy v1.3, so fall back to "clobber"
            # if it doesn't work
            temp_im.writeto(tfile.name, clobber=clobber)
        send_fits_image(img.samp_client, img.samp_key, 'PyBDSM image', tfile.name)
    else:
        # Write image to FITS file
        if outdir is None:
            outdir = img.indir
        if not os.path.exists(outdir) and outdir != '':
            os.makedirs(outdir)
        outfilename = os.path.join(outdir, filename)
        if os.path.isfile(outfilename):
            if clobber:
                os.remove(outfilename)
            else:
                return
        if os.path.isdir(outfilename):
            if clobber:
                os.system("rm -rf "+outfilename)
            else:
                return
        temp_im = make_fits_image(N.transpose(image), wcs_obj, img.beam,
            img.frequency, img.equinox, telescope, xmin=xmin, ymin=ymin,
            is_mask=is_mask, shape=(img.shape[1], img.shape[0], image.shape[1],
            image.shape[0]))
        if use == 'rap':
            outfile = outfilename + '.fits'
        else:
            outfile = outfilename
        try:
            temp_im.writeto(outfile,  overwrite=clobber)
        except TypeError:
            # The "overwrite" argument was added in astropy v1.3, so fall back to "clobber"
            # if it doesn't work
            temp_im.writeto(outfile,  clobber=clobber)
        temp_im.close()

        if use == 'rap':
            # For CASA images, read in FITS image and convert
            try:
                import casacore.images as pim
                import casacore.tables as pt
                import os
                outimage = pim.image(outfile)
                outimage.saveas(outfilename, overwrite=clobber)

                # For masks, use the coordinates dictionary from the input
                # image, as this is needed in order for the
                # image to work as a clean mask in CASA.
                if is_mask:
                    if img.coords_dict is None:
                        mylog.warning('Mask header information may be incomplete.')
                    else:
                        outtable = pt.table(outfilename, readonly=False, ack=False)
                        outtable.putkeywords({'coords': img.coords_dict})
                        outtable.done()

            except ImportError as err:
                import os
                os.remove(outfile)
                raise RuntimeError("Error importing python-casacore. CASA image could not "
                                   "be writen. Use img_format = 'fits' instead.")


def make_fits_image(imagedata, wcsobj, beam, freq, equinox, telescope, xmin=0, ymin=0,
                    is_mask=False, shape=None):
    """Makes a simple FITS hdulist appropriate for single-channel images"""
    import numpy as np
    from astropy.io import fits as pyfits

    # If mask, expand to all channels and Stokes for compatibility with casa
    if is_mask and shape is not None:
        shape_out = shape
    else:
        shape_out = [1, 1, imagedata.shape[0], imagedata.shape[1]]
    hdu = pyfits.PrimaryHDU(np.resize(imagedata, shape_out))
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header

    # Add WCS info
    header['CRVAL1'] = wcsobj.wcs.crval[0]
    header['CDELT1'] = wcsobj.wcs.cdelt[0]
    header['CRPIX1'] = wcsobj.wcs.crpix[0] + xmin
    header['CUNIT1'] = str(wcsobj.wcs.cunit[0]).strip().lower() # needed due to bug in pywcs/astropy
    header['CTYPE1'] = wcsobj.wcs.ctype[0]
    header['CRVAL2'] = wcsobj.wcs.crval[1]
    header['CDELT2'] = wcsobj.wcs.cdelt[1]
    header['CRPIX2'] = wcsobj.wcs.crpix[1] + ymin
    header['CUNIT2'] = str(wcsobj.wcs.cunit[1]).strip().lower()
    header['CTYPE2'] = wcsobj.wcs.ctype[1]

    # Add STOKES info
    header['CRVAL3'] = 1.0
    header['CDELT3'] = 1.0
    header['CRPIX3'] = 1.0
    header['CUNIT3'] = ''
    header['CTYPE3'] = 'STOKES'

    # Add frequency info
    header['RESTFRQ'] = freq
    header['CRVAL4'] = freq
    header['CDELT4'] = 3e8
    header['CRPIX4'] = 1.0
    header['CUNIT4'] = 'HZ'
    header['CTYPE4'] = 'FREQ'
    header['SPECSYS'] = 'TOPOCENT'

    # Add beam info
    if not is_mask:
        header['BMAJ'] = beam[0]
        header['BMIN'] = beam[1]
        header['BPA'] = beam[2]

    # Add equinox
    header['EQUINOX'] = equinox

    # Add telescope
    if telescope is not None:
        header['TELESCOP'] = telescope

    hdulist[0].header = header
    return hdulist

def retrieve_map(img, map_name):
    """Returns a map cached on disk."""
    import numpy as N
    import os

    filename = get_name(img, map_name)
    if not os.path.isfile(filename):
        return None
    infile = open(filename, 'rb')
    data = N.load(infile)
    infile.close()
    return data

def store_map(img, map_name, map_data):
    """Caches a map to disk."""
    import numpy as N

    filename = get_name(img, map_name)
    outfile = open(filename, 'wb')
    N.save(outfile, map_data)
    outfile.close()

def del_map(img, map_name):
    """Deletes a cached map."""
    import os

    filename = get_name(img, map_name)
    if os.path.isfile(filename):
        os.remove(filename)

def get_name(img, map_name):
    """Returns name of cache file."""
    import os

    if img._pi:
        pi_text = 'pi'
    else:
        pi_text = 'I'
    suffix = '/w%i_%s/' % (img.j, pi_text)
    dir = img.tempdir + suffix
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir + map_name + '.bin'

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
    points = N.asarray(N.where(mask ^ nd.binary_erosion(mask)))
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

    if labels is not None:
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
    leftout = ~mask ^ open

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
    if n_subisl3 > 0 and isl_pixs3 is not None:                                 # open 3 breaks up island
        max_sub3 = N.max(isl_pixs3)
        if max_sub3 < frac_bigisl3 : open3 = True       # if biggest sub island isnt too big
    if n_subisl5 > 0 and isl_pixs5 is not None:                                 # open 5 breaks up island
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
            print("\033[31;1mWARNING\033[0m: Matplotlib not found. Plotting disabled.")
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

    if aperture_pix is None:
        return [0.0, 0.0]

    # Make ch0 and rms subimages
    ch0 = img.ch0_arr
    shape = ch0.shape
    xlo = int(posn_pix[0]) - int(aperture_pix) - 1
    if xlo < 0:
        xlo = 0
    xhi = int(posn_pix[0]) + int(aperture_pix) + 1
    if xhi > shape[0]:
        xhi = shape[0]
    ylo = int(posn_pix[1]) - int(aperture_pix) - 1
    if ylo < 0:
        ylo = 0
    yhi = int(posn_pix[1]) + int(aperture_pix) + 1
    if yhi > shape[1]:
        yhi = shape[1]

    mean = img.mean_arr
    rms = img.rms_arr
    aper_im = ch0[int(xlo):int(xhi), int(ylo):int(yhi)] - mean[int(xlo):int(xhi), int(ylo):int(yhi)]
    aper_rms = rms[int(xlo):int(xhi), int(ylo):int(yhi)]
    posn_pix_new = [int(posn_pix[0])-xlo, int(posn_pix[1])-ylo]
    pixel_beamarea = img.pixel_beamarea()
    aper_flux = aperture_flux(aperture_pix, posn_pix_new, aper_im, aper_rms, pixel_beamarea)
    return aper_flux

def aperture_flux(aperture_pix, posn_pix, aper_im, aper_rms, beamarea):
    """Returns aperture flux and error"""
    import numpy as N

    dist_mask = generate_aperture(aper_im.shape[0], aper_im.shape[1], posn_pix[0], posn_pix[1], aperture_pix)
    aper_mask = N.where(dist_mask.astype(bool))
    if N.size(aper_mask) == 0:
        return [0.0, 0.0]
    aper_flux = N.nansum(aper_im[aper_mask])/beamarea # Jy
    pixels_in_source = N.sum(~N.isnan(aper_im[aper_mask])) # number of unmasked pixels assigned to current source
    aper_fluxE = nanmean(aper_rms[aper_mask]) * N.sqrt(pixels_in_source/beamarea) # Jy
    return [aper_flux, aper_fluxE]

def generate_aperture(xsize, ysize, xcenter, ycenter, radius):
    """Makes a mask (1 = inside aperture) for a circular aperture"""
    import numpy

    x, y = numpy.mgrid[0.5:xsize, 0.5:ysize]
    mask = ((x - xcenter)**2 + (y - ycenter)**2 <= radius**2) * 1
    return mask

def make_src_mask(mask_size, posn_pix, aperture_pix):
    """Makes an island mask (1 = inside aperture) for a given source position.
    """
    import numpy as N

    xsize, ysize = mask_size
    if aperture_pix is None:
        return N.zeros((xsize, ysize), dtype=int)

    # Make subimages
    xlo = int(posn_pix[0]-int(aperture_pix)-1)
    if xlo < 0:
        xlo = 0
    xhi = int(posn_pix[0]+int(aperture_pix)+1)
    if xhi > xsize:
        xhi = xsize
    ylo = int(posn_pix[1]-int(aperture_pix)-1)
    if ylo < 0:
        ylo = 0
    yhi = int(posn_pix[1]+int(aperture_pix)+1)
    if yhi > ysize:
        yhi = ysize

    mask = N.zeros((xsize, ysize), dtype=int)
    posn_pix_new = [posn_pix[0]-xlo, posn_pix[1]-ylo]
    submask_xsize = xhi - xlo
    submask_ysize = yhi - ylo
    submask = generate_aperture(submask_xsize, submask_ysize, posn_pix_new[0], posn_pix_new[1], aperture_pix)
    submask_slice = [slice(int(xlo), int(xhi)), slice(int(ylo), int(yhi))]
    mask[tuple(submask_slice)] = submask
    return mask

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
    try:
        # Python 3
        from xmlrpc.client import ServerProxy
    except ImportError:
        # Python 2
        from xmlrpclib import ServerProxy

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
    s = ServerProxy(HUB_PARAMS['samp.hub.xmlrpc.url'])

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

def make_curvature_map(subim):
    """Makes a curvature map with the Aegean curvature algorithm
    (Hancock et al. 2012)

    The Aegean algorithm uses a curvature map to identify regions of negative
    curvature. These regions then define distinct sources.
    """
    import scipy.signal as sg
    import numpy as N
    import sys

    # Make average curavature map:
    curv_kernal = N.array([[1, 1, 1],[1, -8, 1],[1, 1, 1]])
    # The next step prints meaningless warnings, so suppress them
    original_stdout = sys.stdout  # keep a reference to STDOUT
    sys.stdout = NullDevice()  # redirect the real STDOUT
    curv_map = sg.convolve2d(subim, curv_kernal)
    sys.stdout = original_stdout  # turn STDOUT back on

    return curv_map


def bstat(indata, mask, kappa_npixbeam):
    """Numpy version of the c++ bstat routine

    Uses the PySE method for calculating the clipped mean and rms of an array.
    This method is superior to the c++ bstat routine (see section 2.7.3 of
    http://dare.uva.nl/document/174052 for details) and, since the Numpy
    functions used here are written in c, there should be no big computational
    penalty in using Python code.
    """
    import numpy
    from scipy.special import erf, erfcinv

    # Flatten array
    skpix = indata.flatten()
    if mask is not None:
        msk_flat = mask.flatten()
        unmasked = numpy.where(~msk_flat)
        skpix = skpix[unmasked]

    ct = skpix.size
    iter = 0
    c1 = 1.0
    c2 = 0.0
    maxiter = 200
    converge_num = 1e-6
    m_raw = numpy.mean(skpix)
    r_raw = numpy.std(skpix, ddof=1)

    while (c1 >= c2) and (iter < maxiter):
        npix = skpix.size
        if kappa_npixbeam > 0.0:
            kappa = kappa_npixbeam
        else:
            npixbeam = abs(kappa_npixbeam)
            kappa = numpy.sqrt(2.0)*erfcinv(1.0 / (2.0*npix/npixbeam))
            if kappa < 3.0:
                kappa = 3.0
        lastct = ct
        medval = numpy.median(skpix)
        sig = numpy.std(skpix)
        wsm = numpy.where(abs(skpix-medval) < kappa*sig)
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]

        c1 = abs(ct - lastct)
        c2 = converge_num * lastct
        iter += 1

    mean  = numpy.mean(skpix)
    median = numpy.median(skpix)
    sigma = numpy.std(skpix, ddof=1)
    mode = 2.5*median - 1.5*mean
    if sigma > 0.0:
        skew_par = abs(mean - median)/sigma
    else:
        raise RuntimeError("A region with an unphysical rms value has been found. "
            "Please check the input image.")

    if skew_par <= 0.3:
        m = mode
    else:
        m = median

    r1 = numpy.sqrt(2.0*numpy.pi)*erf(kappa/numpy.sqrt(2.0))
    r = numpy.sqrt(sigma**2 * (r1 / (r1 - 2.0*kappa*numpy.exp(-kappa**2/2.0))))

    return m_raw, r_raw, m, r, iter


def centered(arr, newshape):
    """Return the center newshape portion of the array

    This function is a copy of the private _centered() function in
    scipy.signal.signaltools
    """
    import numpy as np

    newshape = np.asarray(newshape)
    currshape = np.array(arr.shape)
    startind = (currshape - newshape) // 2
    endind = startind + newshape
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


def set_up_output_paths(opts):
    """Returns various paths and filenames related to output

    The opts input is either an instance of <class 'bdsf.opts.Opts'> or a
    dict generated by that class.

    The outputs are:
        - parentname: the name of the image, with the path and extension removed
          (if it is a common image extension)
        - output_basedir: the output directory, where the log file and
          other optional outputs of the process_image task are placed
    """
    import os

    # Get filename and outdir from opts
    if type(opts) is dict:
        filename = opts['filename']
        outdir = opts['outdir']
    else:
        # opts is type <class 'bdsf.opts.Opts'>, so options are stored
        # as attributes
        filename = opts.filename
        outdir = opts.outdir

    # Try to trim common extensions from filename to make the parent filename,
    # used for various output purposes
    root, ext = os.path.splitext(filename)
    if ext in ['.fits', '.FITS', '.image']:
        fname = root
    elif ext in ['.gz', '.GZ']:
        root2, ext2 = os.path.splitext(root)
        if ext2 in ['.fits', '.FITS', '.image']:
            fname = root2
        else:
            fname = root
    else:
        fname = filename
    parentname = os.path.basename(fname)

    # Determine the base output directory
    if outdir is None:
        output_basedir = os.path.abspath(os.path.dirname(filename))
    else:
        output_basedir = os.path.abspath(outdir)

    # Make the output directory if needed
    if not os.path.exists(output_basedir):
        os.makedirs(output_basedir)

    # Check that we have write permission to the base directory
    if not os.access(output_basedir, os.W_OK):
        raise RuntimeError("Cannot write to the output directory '{0}' (permission denied). "
                           "Please specify an output directory to which you have "
                           "write permission using the 'outdir' option.".format(output_basedir))

    return parentname, output_basedir

def fix_gaussian_axes(major, minor, pa):
    """Check a Gaussian for switched axes and fix if found

    Returns corrected (major, minor, pa)
    """
    if major < minor:
        major, minor = minor, major
        pa += 90.0
    pa = divmod(pa, 180)[1]  # restrict to range [0, 180)

    return (major, minor, pa)
