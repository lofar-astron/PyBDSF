
def shapelet_coeff(nmax=20,basis='cartesian'):
    """ Computes shapelet coefficient matrix for cartesian and polar 
      hc=shapelet_coeff(nmax=10, basis='cartesian') or
      hc=shapelet_coeff(10) or hc=shapelet_coeff().
      hc(nmax) will be a nmax+1 X nmax+1 matrix."""
    import numpy as N

    order=nmax+1
    if basis == 'polar':
       raise NotImplementedError, "Polar shapelets not yet implemented."

    hc=N.zeros([order,order])
    hnm1=N.zeros(order); hn=N.zeros(order)

    hnm1[0]=1.0; hn[0]=0.0; hn[1]=2.0
    hc[0]=hnm1
    hc[1]=hn
    for ind in range(3,order+1):
        n=ind-2
        hnp1=-2.0*n*hnm1
        hnp1[1:] += 2.0*hn[:order-1]
        hc[ind-1]=hnp1
        hnm1=hn
        hn=hnp1

    return hc


