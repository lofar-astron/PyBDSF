
import numpy as N
from math import *
import matplotlib.pyplot as pl
#import pylab as pl

def readinfile(fname, skip=0):
    """ Read the rows in fname and returns a column-wise list. If the column is a 
    number then it is a numpy array, else a list of strings. It ignores skip number 
    of lines at the beginning. Default skip=0"""

    f=open(fname, 'r')
    dlist = []
    for i,line in enumerate(f):
        if i >= skip:
          dlist.append(line.split())
    f.close()

    new= [[r[col] for r in dlist] for col in range(len(dlist[0]))]

    dlist = []
    for col in new:
        ele = col[0]
        try:
            x=float(ele)
            dlist.append(N.asarray(col, dtype=float))
        except:
            dlist.append(col)
             
    return dlist

def justdist(ra1,ra2,dec1,dec2):
    """ Computes distance on surface of a sphere in arcsec 
    between (ra1,dec1) and (ra2,dec2) in degrees. """

    rad=pi/180.0

    r1=(ra1-180.0)*rad
    r2=(ra2-180.0)*rad
    d1=dec1*rad
    d2=dec2*rad
    dist=acos(sin(d1)*sin(d2)+cos(d1)*cos(d2)*cos(r1-r2))*3600.0/rad  

    return dist

# Fake an ellipse using an N-sided polygon
def Ellipse((cenx,ceny), (radx, rady), alpha, nbin=50, **kwargs):

    theta = 2*pi/nbin*N.arange(nbin) 
    xs = cenx + radx * N.cos(theta)
    ys = ceny + rady * N.sin(theta)

    return pl.Polygon(zip(xs, ys), closed=True, alpha=alpha, **kwargs)

def convertcoords(params):
    """ Convert ra, dec in deg to hh,mm etc and vice versa. Input is a list."""

    if len(params) not in [2, 6]:
      print "convertcoords: Coordinates dont make sense"
    else:
      if len(params) == 2:
        ra, dec = params 
        dumr = ra/15.0
        hh = int(dumr); dum = (dumr - hh)*60.0
        mm = int(dum); ss = (dum - mm)*60.0
        sgn = 1
        if dec < 0: sgn = -1
        dumr = abs(dec)
        dd = int(dumr); dum = (dumr - dd)*60.0
        ma = int(dum); sa = (dum - ma)*60.0
        return hh, mm, ss, sgn, dd, ma, sa
      if len(params) == 6:
        hh, mm, ss, dd, ma, sa = params
        sgn = 1
        if dd < 0 or ma < 0 or sa < 0: sgn = -1
        ra = (hh+mm/60.0+ss/3600.0)*15.0
        dec = sgn*(abs(dd)+abs(ma)/60.0+abs(sa)/3600.0)
        return ra, dec

def B1950toJ2000(Bcoord):
    """ Precess using Aoki et al. 1983. Same results as NED to ~0.2asec """
    from math import sin, cos, pi, sqrt, asin, acos
    import numpy as N

    rad = 180.0/pi
    ra, dec = Bcoord

    A = N.array([-1.62557e-6, -0.31919e-6, -0.13843e-6])
    M = N.array([[0.9999256782, 0.0111820609, 0.00485794], [-0.0111820610, 0.9999374784, -0.0000271474], \
                 [-0.0048579477, -0.0000271765, 0.9999881997]])

    r0=N.zeros(3)
    r0[0]=cos(dec/rad)*cos(ra/rad)
    r0[1]=cos(dec/rad)*sin(ra/rad)
    r0[2]=sin(dec/rad)

    r0A=N.sum(r0*A)
    r1=r0-A+r0A*r0
    r = N.sum(M.transpose()*r1, axis = 1)

    rscal = sqrt(N.sum(r*r))
    decj=asin(r[2]/rscal)*rad 

    d1=r[0]/rscal/cos(decj/rad)
    d2=r[1]/rscal/cos(decj/rad)
    raj=acos(d1)*rad 
    if d2 < 0.0: raj = 360.0 - raj

    Jcoord = [raj, decj]
    return Jcoord

