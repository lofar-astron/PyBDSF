"""
This is for pybdsm for calculating spectral index. We assume a linear spectral index
in log(freq) and then each channel has a flux which is bit wrong because of the colour
correction problem within that band. 

Now we average n such channels. There will be another error made, partly because of the
colour correction now for channels (including discretisation) and the colour correction 
of the earlier 2nd order colour correction. 

This is to see how much they differ. Refer notebook for forumlae.

"""

import numpy as N
import pylab as pl
import math

nchan = N.array([9, 17])
alpha_arr = N.arange(-1.3, -0.3, 0.1)
deltanu = N.array([0.05e6, 0.1e6, 0.2e6])

freq = N.arange(40.0e6, 200.0e6, 10.0e6)

pl1 = pl.figure()
pl2 = pl.figure()
pl3 = pl.figure()
k = 0
for inchan, n in enumerate(nchan):
    for ibw, bw in enumerate(deltanu):
        k += 1
        for ia, alpha in enumerate(alpha_arr):
            f_diff1 = N.zeros(len(freq))
            f_diff2 = N.zeros(len(freq))
            for ifreq, f in enumerate(freq):
                f_arr = N.arange(f-(n-1)/2*bw, f+(n+1)/2*bw, bw)
                f_naive = N.mean(f_arr)
                f1 = N.power(f_arr, alpha)
                f2 = N.power(f_arr, alpha-2.0)
            
                f1 = 1.0/n*N.sum(f1)
                f2 = 1.0/n*N.sum(f2)*bw*bw*alpha*(alpha-1.0)/24.0

                f_eff1 = N.power(f1, 1.0/alpha)
                f_eff2 = N.power(f1+f2, 1.0/alpha)

                f_diff1[ifreq] = f_naive - f_eff2
                f_diff2[ifreq] = f_eff1 - f_eff2

            fig = pl.figure(pl1.number)
            adjustprops = dict(wspace=0.5, hspace=0.5) 
            fig.subplots_adjust(**adjustprops)                        
            ax = pl.subplot(2,3,k)
            pl.plot(freq/1e6, f_diff1/1e3)
            pl.title('n='+str(n)+'; bw='+str(bw/1e6)+' MHz')
            pl.xlabel('Freq(MHz)')
            pl.ylabel('Diff in freq (kHz)')
            pl.setp(ax.get_xticklabels(), rotation='vertical', fontsize=12)

            fig = pl.figure(pl2.number)
            adjustprops = dict(wspace=0.5, hspace=0.5) 
            fig.subplots_adjust(**adjustprops)                        
            ax = pl.subplot(2,3,k)
            pl.plot(freq/1e6, f_diff2)
            pl.title('n='+str(n)+'; bw='+str(bw/1e6)+' MHz')
            pl.xlabel('Freq(MHz)')
            pl.ylabel('Diff due to 2nd order (Hz)')
            pl.setp(ax.get_xticklabels(), rotation='vertical', fontsize=12)
        
            fig = pl.figure(pl3.number)
            adjustprops = dict(wspace=0.9, hspace=0.5) 
            fig.subplots_adjust(**adjustprops)                        
            ax = pl.subplot(2,3,k)
            f2 = f_naive+5e6
            y = f_diff1*alpha/f_naive/math.log(f_naive/(f2))
            pl.plot(freq/1e6, y)
            pl.title('n='+str(n)+'; bw='+str(bw/1e6)+' MHz')
            pl.xlabel('Freq(MHz)')
            pl.ylabel('Error in sp.in. for f2=f1+10MHz')
            pl.setp(ax.get_xticklabels(), rotation='vertical', fontsize=12)

pl.figure(pl1.number)
pl.savefig('colourcorr_full.png')
pl.figure(pl2.number)
pl.savefig('colourcorr_order1-2.png')
pl.figure(pl3.number)
pl.savefig('colourcorr_delta_spin.png')


