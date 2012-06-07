
"""
This is to write the output of pyBDSM in the format of fBDSM files such that the results
can be plotted using Anaamika, with _fbdsm_plot.so.

This writes the islandlist files and gaul and srl files.

"""

import numpy as N
import mylogger
#import _py2fbdsm as fbp


def write_fbdsm_islands(img):
    """ Writes islandlist and islandlist.asc. """

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Op_fBDSM  ")
    nisl = len(img.islands)
    n, m = img.ch0.shape[0:2]
    maxmem = 0
    for isl in img.islands: maxmem = max(maxmem, isl.size_active) 
    imagename = img.imagename
    runcode = 'aq'
    scratch = img.opts.fbdsm_scratch
    arr = img.ch0
                                        # construct the rank array
    island = N.zeros((n,m), int)
    for isl in img.islands: 
        ori = isl.origin
        x, y = N.where(~isl.mask_active)
        x += ori[0]; y += ori[1]
        for i in range(len(x)):
            island[x[i],y[i]] = isl.island_id + 1
 
    fbp.sub_iland_mat2list(imagename,nisl,maxmem,runcode,scratch,island,arr)

    fn = imagename + '.rank'
    str1 = ' '+str(nisl)+' '+str(maxmem)
    fbp.writearray_bin_int(island,n,m,fn,str1,runcode)
  
    mylog.info('Wrote '+imagename+'.rank.nmg and '+imagename+'.islandlist(.asc) in '+scratch)

##################################################################################################

def write_fbdsm_gaul(img):
    """ Writes gaul and gaul.bin files. Write the gaul file first and then call a fortran program
        to write gaul.bin. Also writes the header file."""
    import functions as func
    from const import fwsig

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Op_fBDSM  ")
    scratch = img.opts.fbdsm_scratch
    imagename = img.imagename
    n, m = img.ch0.shape[0:2]
    maxgaus = 0
    for isl in img.islands:
        maxgaus = max(maxgaus, isl.ngaus)
    rmsmap = {True : 'map', False: 'const'}
    fa = open(scratch+imagename+'.gaul', 'w')

    fa.write(' Gaussian list made by PyBDSM'+'\n')
    fa.write(' Image_name '+img.filename+'\n')
    fa.write(' Image_size_x '+str(n)+'\n')
    fa.write(' Image_size_y '+str(m)+'\n')
    fa.write(' Island_list_name '+imagename+'\n')
    fa.write(' Number_of_islands '+str(len(img.islands))+'\n')
    fa.write(' Number_of_gaussians '+str(img.ngaus)+'\n')
    fa.write(' Max_gaussians_per_island '+str(maxgaus)+'\n')
    fa.write(' RMS_map '+rmsmap[img.opts.rms_map]+'\n')
    fa.write(' Sigma '+str(img.clipped_rms)+'\n')
    fa.write(' Detect_threshold '+str(img.opts.thresh_pix)+'\n')
    fa.write(' gaul# island# flag tot_Jy err peak_Jy err   RA err DEC err  xpos_pix err ypos_pix '+\
             'err bmaj_asec_fw err bmin_asec_fw err bpa_deg err deconv_bmaj_bmin_bpa_asec_fw &errors '+\
             'src_rms src_av isl_rms isl_av spin e_spin src#  blc1  blc2  trc1 trc2  dumr1 dumr2 '+\
             'dumr3 dumr4 dumr5 dumr6'+'\n')
    fa.write(' fmt 113 "(3(i7,1x),4(1Pe11.3,1x),4(0Pf13.9,1x),16(0Pf11.5,1x),4(1Pe11.3,1x),'+\
             '0Pf7.3,1x,0Pf7.3,1x,i7,4(1x,i5),6(1x,1Pe11.3))"\n')

    for g in img.gaussians:#():
        gidx = g.gaus_num
        iidx = g.island_id+1
        A = g.peak_flux
        ra, dec = g.centre_sky
        x, y = g.centre_pix
        x += 1; y += 1
        shape = g.size_sky
        shape[0] = shape[0]*3600; shape[1] = shape[1]*3600
        eA = g.peak_fluxE
        era, edec = g.centre_skyE
        ex, ey = g.centre_pixE
        eshape = g.size_skyE
        eshape[0] = eshape[0]*3600; eshape[1] = eshape[1]*3600

        str1 = 3*"%7i "+4*"%11.3e "+4*"%13.9f "+16*"%11.5f "+4*"%11.3e "+ \
               "%7.3f %7.3f %7i"+4*" %5i"+6*" %11.3e" +"\n"
        str2 = str1 % (gidx, iidx, 0,    0., 0., A, eA, \
               ra, era, dec, edec, x, ex,y, ey, \
               shape[0], eshape[0], shape[1], eshape[1], shape[2], eshape[2], \
               0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,  0.,0.,  g.source_id, \
               0,0,0,0, 0.,0.,0.,0.,0.,0.)
        fa.write(str2)
    fbp.gaul2gaulbin(fa.name)
    fa.close()

    cdelt = N.array(img.wcs_obj.acdelt)
    crpix = N.array(img.wcs_obj.crpix)
    crval = N.array(img.wcs_obj.crval)
    ctype = N.array(img.wcs_obj.ctype)
    for i in range(3): ctype[i] += ' '*(8-len(ctype[i]))
    bm_pix = N.array(img.pixel_beam)*fwsig
    if img.wcs_obj.crota == None:
      crota = N.array(img.wcs_obj.crota)
    else:
      crota = N.zeros(3)
    fbp.writefitshead(imagename,crpix,ctype[0],ctype[1],ctype[2],cdelt,crval,crota,bm_pix,scratch)
    fbp.writefitshead(img.filename,crpix,ctype[0],ctype[1],ctype[2],cdelt,crval,crota,bm_pix,scratch)

    mylog.info('Wrote '+imagename+'.gaul, '+imagename+'.gaul.bin and '+imagename+'.header in '+scratch)

##################################################################################################



