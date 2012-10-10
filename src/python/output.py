"""Module output.

Defines functions that write the results of source detection in a
variety of formats. These are then used as methods of Image objects
and/or are called by the outlist operation if output_all is True.
"""
from image import Op

class Op_outlist(Op):
    """Write out list of Gaussians

    Currently 6 output formats are supported:
    - FIST list
    - BBS list
    - star list
    - kvis annotations
    - ascii list
    - ds9 region list

    All output lists are generated atm.
    """
    def __call__(self, img):
        if img.opts.output_all:
            import os
            if len(img.gaussians) > 0:
                dir = img.basedir + '/catalogues/'
                if not os.path.exists(dir):
                    os.mkdir(dir)
                self.write_bbs(img, dir)
                self.write_gaul(img, dir)
                self.write_srl(img, dir)
                self.write_aips(img, dir)
                self.write_kvis(img, dir)
                self.write_ds9(img, dir)
                self.write_gaul_FITS(img, dir)
                self.write_srl_FITS(img, dir)
            if not os.path.exists(img.basedir + '/misc/'):
                os.mkdir(img.basedir + '/misc/')
            self.write_opts(img, img.basedir + '/misc/')
            self.save_opts(img, img.basedir + '/misc/')
            img.completed_Ops.append('outlist')

    def write_bbs(self, img, dir):
        """ Writes the gaussian list as a bbs-readable file"""
        prefix = ''
        if img.extraparams.has_key('bbsprefix'):
            prefix = img.extraparams['bbsprefix']+'_'
        if img.extraparams.has_key('bbsname'):
            name = img.extraparams['bbsname']
        else:
            name = img.imagename
        fname = dir + name + '.sky_in'

        # Write Gaussian list
        write_bbs_gaul(img, filename=fname, srcroot=img.opts.srcroot,
                       patch=img.opts.bbs_patches,
                       sort_by='flux', clobber=True)


    def write_gaul(self, img, dir):
        """ Writes the gaussian list as an ASCII file"""
        fname = dir + img.imagename + '.gaul'
        write_ascii_list(img, filename=fname, sort_by='indx',
                         clobber=True, objtype='gaul')

    def write_srl(self, img, dir):
        """ Writes the source list as an ASCII file"""
        fname = dir + img.imagename + '.srl'
        write_ascii_list(img, filename=fname, sort_by='indx',
                         clobber=True, objtype='srl')

    def write_aips(self, img, dir):
        """ Writes the gaussian list an AIPS STAR file"""
        fname = dir + img.imagename + '.star'
        write_star(img, filename=fname, sort_by='indx',
                   clobber=True)

    def write_kvis(self, img, dir):
        """ Writes the gaussian list as a kvis file"""
        fname = dir + img.imagename + '.kvis.ann'
        write_kvis_ann(img, filename=fname, sort_by='indx',
                       clobber=True)

    def write_ds9(self, img, dir):
        """ Writes the gaussian list as a ds9 region file"""
        fname = dir + img.imagename + '.ds9.reg'
        write_ds9_list(img, filename=fname, srcroot=img.opts.srcroot,
                       clobber=True, deconvolve=False)

    def write_gaul_FITS(self, img, dir):
        """ Writes the gaussian list as FITS binary table"""
        fname = dir + img.imagename+'.gaul.FITS'
        write_fits_list(img, filename=fname, sort_by='indx',
                        clobber=True, objtype='gaul')

    def write_srl_FITS(self, img, dir):
        """ Writes the source list as FITS binary table"""
        fname = dir + img.imagename+'.srl.FITS'
        write_fits_list(img, filename=fname, sort_by='indx',
                        clobber=True, objtype='srl')

    def write_shap_FITS(self, img, dir):
        """ Writes the shapelet list as a FITS file"""
        fname = dir + img.imagename + '.shap.FITS'
        write_fits_list(img, filename=fname, sort_by='indx',
                        clobber=True, objtype='shap')

    def write_opts(self, img, dir):
        """ Writes input parameters to a text file."""
        import inspect
        import types
        import mylogger

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
        fname = 'parameters_used'
        f = open(dir+fname, 'w')
        mylog.info('Writing '+dir+fname)
        for attr in inspect.getmembers(img.opts):
          if attr[0][0] != '_':
            if isinstance(attr[1], (int, str, bool, float, types.NoneType, tuple, list)):
              f.write('%-40s' % attr[0])
              f.write(repr(attr[1])+'\n')

              # Also print the values derived internally. They are all stored
              # in img with the same name (e.g., img.opts.beam --> img.beam)
              if hasattr(img, attr[0]):
                  used = img.__getattribute__(attr[0])
                  if used != attr[1] and isinstance(used, (int, str, bool, float,
                                                           types.NoneType, tuple,
                                                           list)):
                      f.write('%-40s' % '    Value used')
                      f.write(repr(used)+'\n')
        f.close()

    def save_opts(self, img, dir):
        """ Saves input parameters to a PyBDSM save file."""
        import interface
        import mylogger

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
        fname = 'parameters.sav'
        mylog.info('Writing '+dir+fname)
        interface.save_pars(img, dir+fname, quiet=True)


def ra2hhmmss(deg):
    """Convert RA coordinate (in degrees) to HH MM SS"""

    from math import modf
    if deg < 0:
        deg += 360.0
        #raise RuntimeError("Negative RA")
    x, hh = modf(deg/15.)
    x, mm = modf(x*60)
    ss = x*60

    return (int(hh), int(mm), ss)

def dec2ddmmss(deg):
    """Convert DEC coordinate (in degrees) to DD MM SS"""

    from math import modf
    sign = (-1 if deg < 0 else 1)
    x, dd = modf(abs(deg))
    x, ma = modf(x*60)
    sa = x*60

    return (int(dd), int(ma), sa, sign)

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

def write_bbs_gaul(img, filename=None, srcroot=None, patch=None,
                   incl_primary=True, sort_by='flux',
                   clobber=False):
    """Writes Gaussian list to a BBS sky model"""
    import numpy as N
    from const import fwsig
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM.write_gaul")
    if int(img.equinox) != 2000 and int(img.equinox) != 1950:
        mylog.warning('Equinox of input image is not J2000 or B1950. '\
                          'Sky model may not be appropriate for BBS.')
    if int(img.equinox) == 1950:
        mylog.warning('Equinox of input image is B1950. Coordinates '\
                          'will be precessed to J2000.')

    outl, outn, patl = list_and_sort_gaussians(img, patch=patch,
                                               root=srcroot, sort_by=sort_by)
    outstr_list = make_bbs_str(img, outl, outn, patl)

    if filename == None:
        filename = img.imagename + '.sky_in'
    if os.path.exists(filename) and clobber == False:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, 'w')
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_lsm_gaul(img, filename=None, srcroot=None, patch=None,
                   incl_primary=True, sort_by='flux',
                   clobber=False):
    """Writes Gaussian list to a Sagecal lsm sky model"""
    import numpy as N
    from const import fwsig
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM.write_gaul")
    if int(img.equinox) != 2000 and int(img.equinox) != 1950:
        mylog.warning('Equinox of input image is not J2000 or B1950. '\
                          'Sky model may not be appropriate for Sagecal.')
    if int(img.equinox) == 1950:
        mylog.warning('Equinox of input image is B1950. Coordinates '\
                          'will be precessed to J2000.')

    outl, outn, patl = list_and_sort_gaussians(img, patch=patch,
                                               root=srcroot, sort_by=sort_by)
    outstr_list = make_lsm_str(img, outl, outn)

    if filename == None:
        filename = img.imagename + '.lsm'
    if os.path.exists(filename) and clobber == False:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, 'w')
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_ds9_list(img, filename=None, srcroot=None, deconvolve=False,
                   clobber=False, objtype='gaul'):
    """Writes Gaussian list to a ds9 region file"""
    import numpy as N
    from const import fwsig
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if objtype == 'gaul':
        outl, outn, patl = list_and_sort_gaussians(img, patch=None)
    elif objtype == 'srl':
        root = img.parentname
        outl = [img.sources]
        outn = []
        for src in img.sources:
            outn.append(root + '_i' + str(src.island_id) + '_s' +
                            str(src.source_id))
        outn = [outn]
    outstr_list = make_ds9_str(img, outl, outn, deconvolve=deconvolve)
    if filename == None:
        filename = img.imagename + '.' + objtype + '.reg'
    if os.path.exists(filename) and clobber == False:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, "w")
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_ascii_list(img, filename=None, sort_by='indx',
                     incl_chan=False, incl_empty=False, clobber=False, objtype='gaul'):
    """Writes Gaussian list to an ASCII file"""
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if objtype == 'gaul':
        outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)
    elif objtype == 'srl':
        outl = [img.sources]
    outstr_list = make_ascii_str(img, outl, objtype=objtype)
    if filename == None:
        if objtype == 'gaul':
            filename = img.imagename + '.gaul'
        elif objtype == 'srl':
            filename = img.imagename + '.srl'
    if os.path.exists(filename) and clobber == False:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, "w")
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_casa_gaul(img, filename=None, clobber=False):
    """Writes a clean box file for use in casapy"""
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    outl, outn, patl = list_and_sort_gaussians(img, patch=None)
    outstr_list = make_casa_str(img, outl)
    if filename == None:
        filename = img.imagename + '.box'
    if os.path.exists(filename) and clobber == False:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, "w")
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_fits_list(img, filename=None, sort_by='index', objtype='gaul',
                    incl_chan=False, incl_empty=False, clobber=False):
    """ Write as FITS binary table.
    """
    import mylogger
    import pyfits
    import os
    import numpy as N
    from _version import __version__, __revision__

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if objtype == 'gaul':
        outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)
    elif objtype == 'srl':
        outl = [img.sources]
    elif objtype == 'shap':
        outl = [img.islands]

    nmax = 0
    if objtype == 'shap':
        # loop over shapelets and get maximum size of coefficient matrix
        for isl in outl[0]:
            if isl.shapelet_nmax > nmax:
                nmax = isl.shapelet_nmax
        nmax += 1

    if img.opts.aperture != None:
        incl_aper = True
    else:
        incl_aper = False
    cvals, cnames, cformats, cunits = make_output_columns(outl[0][0], fits=True,
                                                          objtype=objtype,
                                                          incl_spin=img.opts.spectralindex_do,
                                                          incl_chan=img.opts.incl_chan,
                                                          incl_pol=img.opts.polarisation_do,
                                                          incl_aper=incl_aper,
                                                          incl_empty=incl_empty,
                                                          nmax=nmax, nchan=img.nchan)
    out_list = make_fits_list(img, outl, objtype=objtype, nmax=nmax)
    col_list = []
    for ind, col in enumerate(out_list):
      list1 = pyfits.Column(name=cnames[ind], format=cformats[ind],
                            unit=cunits[ind], array=N.array(out_list[ind]))
      col_list.append(list1)
    tbhdu = pyfits.new_table(col_list)
    if objtype == 'gaul':
        tbhdu.header.add_comment('Gaussian list for '+img.filename)
    elif objtype == 'srl':
        tbhdu.header.add_comment('Source list for '+img.filename)
    elif objtype == 'shap':
        tbhdu.header.add_comment('Shapelet list for '+img.filename)
    tbhdu.header.add_comment('Generated by PyBDSM version %s (LUS revision %s)'
                             % (__version__, __revision__))
    freq = "%.5e" % img.frequency
    tbhdu.header.add_comment('Reference frequency of the detection ("ch0") image: %s Hz' % freq)
    tbhdu.header.add_comment('Equinox : %s' % img.equinox)
    tbhdu.header.update('INIMAGE', img.filename, 'Filename of image')
    tbhdu.header.update('FREQ0', float(freq), 'Reference frequency')
    tbhdu.header.update('EQUINOX', img.equinox, 'Equinox')
    if filename == None:
        filename = img.imagename + '.' + objtype + '.fits'
    if os.path.exists(filename) and clobber == False:
        return None
    mylog.info('Writing ' + filename)
    tbhdu.writeto(filename, clobber=True)
    return filename


def write_kvis_ann(img, filename=None, sort_by='indx',
                   clobber=False):
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if filename == None:
        filename = img.imagename + '.kvis.ann'
    if os.path.exists(filename) and clobber == False:
        return None
    f = open(filename, 'w')
    mylog.info('Writing '+filename)
    f.write("### KVis annotation file\n\n")
    f.write("color green\n\n")

    outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)
    for g in outl[0]:
        iidx = g.island_id
        ra, dec = g.centre_sky
        shape = g.size_sky

        str = 'text   %10.5f %10.5f   %d\n' % \
            (ra, dec, iidx)
        f.write(str)
        str = 'ellipse %10.5f %10.5f   %10.7f %10.7f %10.4f\n' % \
            (ra, dec, shape[0], shape[1], shape[2])
        f.write(str)
    f.close()
    return filename


def write_star(img, filename=None, sort_by='indx',
               clobber=False):
    from output import ra2hhmmss, dec2ddmmss
    import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if filename == None:
        filename = img.imagename + '.star'
    if os.path.exists(filename) and clobber == False:
        return None
    f = open(filename, 'w')
    mylog.info('Writing '+filename)

    outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)

    for g in outl[0]:
        A = g.peak_flux
        ra, dec = g.centre_sky
        shape = g.size_sky
        ### convert to canonical representation
        ra = ra2hhmmss(ra)
        dec= dec2ddmmss(dec)
        decsign = ('-' if dec[3] < 0 else '+')

        str = '%2i %2i %6.3f ' \
              '%c%2i %2i %6.3f ' \
              '%9.4f %9.4f %7.2f ' \
              '%2i %13.7f %10s\n' % \
              (ra[0], ra[1], ra[2],
               decsign, dec[0], dec[1], dec[2],
               shape[0]*3600, shape[1]*3600, shape[2],
               4, A, '')

        f.write(str)
    f.close()
    return filename


def make_bbs_str(img, glist, gnames, patchnames):
    """Makes a list of string entries for a BBS sky model."""
    from output import ra2hhmmss
    from output import dec2ddmmss
    import numpy as N

    outstr_list = []
    freq = "%.5e" % img.frequency
    if patchnames[0] == None:
        outstr_list.append("format = Name, Type, Ra, Dec, I, Q, U, V, "\
                               "MajorAxis, MinorAxis, Orientation, "\
                               "ReferenceFrequency='"+freq+"', "\
                               "SpectralIndex='[]'\n\n")
    else:
        outstr_list.append("format = Name, Type, Patch, Ra, Dec, I, Q, U, V, "\
                               "MajorAxis, MinorAxis, Orientation, "\
                               "ReferenceFrequency='"+freq+"', "\
                               "SpectralIndex='[]'\n\n")
    patchname_last = ''
    for pindx, patch_name in enumerate(patchnames): # loop over patches
      if patch_name != None and patch_name != patchname_last:
          outstr_list.append(', , ' + patch_name + ', 00:00:00, +00.00.00\n')
          patchname_last = patch_name
      gaussians_in_patch = glist[pindx]
      names_in_patch = gnames[pindx]
      for gindx, g in enumerate(gaussians_in_patch):
          if g.gaus_num >= 0 or (g.gaus_num < 0 and img.opts.incl_empty):
              src_name = names_in_patch[gindx]
              ra, dec = g.centre_sky
              if img.equinox == 1950:
                  ra, dec = B1950toJ2000([ra, dec])
              ra = ra2hhmmss(ra)
              sra = str(ra[0]).zfill(2)+':'+str(ra[1]).zfill(2)+':'+str("%.3f" % (ra[2])).zfill(6)
              dec = dec2ddmmss(dec)
              decsign = ('-' if dec[3] < 0 else '+')
              sdec = decsign+str(dec[0]).zfill(2)+'.'+str(dec[1]).zfill(2)+'.'+str("%.3f" % (dec[2])).zfill(6)
              total = str("%.3e" % (g.total_flux))
              deconv = g.deconv_size_sky
              if deconv[0] == 0.0  and deconv[1] == 0.0:
                  stype = 'POINT'
                  deconv[2] = 0.0
              else:
                  stype = 'GAUSSIAN'
              deconv1 = str("%.5e" % (deconv[0]*3600.0))
              deconv2 = str("%.5e" % (deconv[1]*3600.0))
              deconv3 = str("%.5e" % (deconv[2]))
              deconvstr = deconv1 + ', ' + deconv2 + ', ' + deconv3
              specin = '-0.8'
              if hasattr(g, 'spec_indx'):
                  if g.spec_indx != None and N.isfinite(g.spec_indx):
                      specin = str("%.3e" % (g.spec_indx))
              sep = ', '
              if img.opts.polarisation_do:
                  Q_flux = str("%.3e" % (g.total_flux_Q))
                  U_flux = str("%.3e" % (g.total_flux_U))
                  V_flux = str("%.3e" % (g.total_flux_V))
              else:
                  Q_flux = '0.0'
                  U_flux = '0.0'
                  V_flux = '0.0'
              if patch_name == None:
                  outstr_list.append(src_name + sep + stype + sep + sra + sep +
                                     sdec + sep + total + sep + Q_flux + sep +
                                     U_flux + sep + V_flux + sep +
                                     deconvstr + sep + freq + sep +
                                     '[' + specin + ']\n')
              else:
                  outstr_list.append(src_name + sep + stype + sep + patch_name +
                                     sep + sra + sep + sdec + sep + total + sep +
                                     Q_flux + sep + U_flux + sep + V_flux + sep +
                                     deconvstr + sep + freq + sep +
                                     '[' + specin + ']\n')
          else:
            outstr_list.pop()
    return outstr_list

def make_lsm_str(img, glist, gnames):
    """Makes a list of string entries for a BBS sky model."""
    from output import ra2hhmmss
    from output import dec2ddmmss
    from const import fwsig
    import numpy as N

    outstr_list = []
    freq = "%.5e" % img.frequency
    outstr_list.append("## LSM file\n### Name  | RA (hr,min,sec) | DEC (deg,min,sec) | I | Q | U | V | SI | RM | eX | eY | eP | freq0\n\n")
    for gindx, g in enumerate(glist[0]):
        if g.gaus_num >= 0 or (g.gaus_num < 0 and img.opts.incl_empty):
            src_name = gnames[0][gindx]
            ra, dec = g.centre_sky
            if img.equinox == 1950:
                ra, dec = B1950toJ2000([ra, dec])
            ra = ra2hhmmss(ra)
            sra = str(ra[0]).zfill(2)+' '+str(ra[1]).zfill(2)+' '+str("%.3f" % (ra[2])).zfill(6)
            dec = dec2ddmmss(dec)
            decsign = ('-' if dec[3] < 0 else '+')
            sdec = decsign+str(dec[0]).zfill(2)+' '+str(dec[1]).zfill(2)+' '+str("%.3f" % (dec[2])).zfill(6)
            size = g.size_sky #  degrees, in terms of FWHM
            src_area = 1.1331*size[0]*size[1]*fwsig*fwsig*3600.0**2 # area of source in arcsec**2
            peak = str("%.3e" % (g.total_flux/src_area)) # peak flux in Jy/arcsec**2
            deconv = g.deconv_size_sky
            if deconv[0] == 0.0  and deconv[1] == 0.0:
                sname = 'P' + src_name
                deconv[2] = 0.0
            else:
                sname = 'G' + src_name
            deconv1 = str("%.5e" % (deconv[0]*N.pi/180.0/2.0))
            deconv2 = str("%.5e" % (deconv[1]*N.pi/180.0/2.0))
            deconv3 = str("%.5e" % (deconv[2]*N.pi/180.0/2.0))
            deconvstr = deconv1 + ' ' + deconv2 + ' ' + deconv3
            specin = '-0.8'
            if hasattr(g, 'spec_indx'):
                if g.spec_indx != None and N.isfinite(g.spec_indx):
                    specin = str("%.3e" % (g.spec_indx))
            sep = ' '
            if img.opts.polarisation_do:
                Q_flux = str("%.3e" % (g.total_flux_Q/src_area))
                U_flux = str("%.3e" % (g.total_flux_U/src_area))
                V_flux = str("%.3e" % (g.total_flux_V/src_area))
            else:
                Q_flux = '0.0'
                U_flux = '0.0'
                V_flux = '0.0'
            outstr_list.append(sname + sep + sra + sep +
                                   sdec + sep + peak + sep + Q_flux + sep +
                                   U_flux + sep + V_flux + sep +
                                   specin + sep + '0' + sep + deconvstr + sep +
                                   freq + sep + '\n')
    return outstr_list


def make_ds9_str(img, glist, gnames, deconvolve=False):
    """Makes a list of string entries for a ds9 region file."""
    outstr_list = []
    freq = "%.5e" % img.frequency
    if img.equinox == None:
        equinox = 'fk5'
    else:
        if int(img.equinox) == 2000:
            equinox = 'fk5'
        elif int(img.equinox) == 1950:
            equinox = 'fk4'
        else:
            mylog.warning('Equinox of input image is not J2000 or B1950. '\
                                  'Regions may not be correct.')
            equinox = 'fk5'

    outstr_list.append('# Region file format: DS9 version 4.0\nglobal color=green '\
                           'font="helvetica 10 normal" select=1 highlite=1 edit=1 '\
                           'move=1 delete=1 include=1 fixed=0 source\n'+equinox+'\n')

    for gindx, g in enumerate(glist[0]):
        if g.gaus_num >= 0 or (g.gaus_num < 0 and img.opts.incl_empty):
            src_name = gnames[0][gindx]
            try:
                ra, dec = g.centre_sky
            except AttributeError:
                ra, dec = g.posn_sky_centroid
            if deconvolve:
                deconv = g.deconv_size_sky
            else:
                deconv = g.size_sky
            if deconv[0] == 0.0 and deconv[1] == 0.0:
                stype = 'POINT'
                deconv[2] = 0.0
                region = 'point(' + str(ra) + ',' + str(dec) + \
                    ') # point=cross width=2 text={' + src_name + '}\n'
            else:
                # ds9 can't handle 1-D Gaussians, so make sure they are 2-D
                if deconv[0] < 1.0/3600.0: deconv[0] = 1.0/3600.0
                if deconv[1] < 1.0/3600.0: deconv[1] = 1.0/3600.0
                stype = 'GAUSSIAN'
                region = 'ellipse(' + str(ra) + ',' + str(dec) + ',' + \
                    str(deconv[0]*3600.0) + '",' + str(deconv[1]*3600.0) + \
                    '",' + str(deconv[2]+90.0) + ') # text={' + src_name + '}\n'
            outstr_list.append(region)
    return outstr_list


def make_ascii_str(img, glist, objtype='gaul'):
    """Makes a list of string entries for an ascii region file."""
    from _version import __version__, __revision__
    outstr_list = []
    freq = "%.5e" % img.frequency

    if objtype == 'gaul':
        outstr_list.append('# Gaussian list for '+img.filename+'\n')
    elif objtype == 'srl':
        outstr_list.append('# Source list for '+img.filename+'\n')
    outstr_list.append('# Generated by PyBDSM version %s (LUS revision %s)\n'
                       % (__version__, __revision__))
    outstr_list.append('# Reference frequency of the detection ("ch0") image: %s Hz\n' % freq)
    outstr_list.append('# Equinox : %s \n\n' % img.equinox)
    val_list = []
    if img.opts.aperture != None:
        incl_aper = True
    else:
        incl_aper = False

    for i, g in enumerate(glist[0]):
        cvals, cnames, cformats, cunits = make_output_columns(g, fits=False,
                                                              objtype=objtype,
                                                              incl_spin=img.opts.spectralindex_do,
                                                              incl_chan=img.opts.incl_chan,
                                                              incl_pol=img.opts.polarisation_do,
                                                              incl_aper=incl_aper,
                                                              incl_empty = img.opts.incl_empty,
                                                              nchan=img.nchan)
        if cvals != None:
            cformats[-1] += "\n"
            if i == 0:
                outstr_list.append("# " + " ".join(cnames) + "\n")
            outstr_list.append(" ".join(cformats) % tuple(cvals))
    return outstr_list


def make_fits_list(img, glist, objtype='gaul', nmax=30):
    import functions as func

    out_list = []
    if img.opts.aperture != None:
        incl_aper = True
    else:
        incl_aper = False
    for g in glist[0]:
        cvals, ext1, ext2, ext3 = make_output_columns(g, fits=True, objtype=objtype,
                                                      incl_spin=img.opts.spectralindex_do,
                                                      incl_chan=img.opts.incl_chan,
                                                      incl_pol=img.opts.polarisation_do,
                                                      incl_aper=incl_aper,
                                                      incl_empty=img.opts.incl_empty,
                                                      nmax=nmax, nchan=img.nchan)
        if cvals != None:
            out_list.append(cvals)
    out_list = func.trans_gaul(out_list)
    return out_list


def make_casa_str(img, glist):
    """Makes a list of string entries for a casa clean box file."""
    import functions as func
    outstr_list = []
    sep = ' '
    scale = 2.0
    for gindx, g in enumerate(glist[0]):
        x, y = g.centre_pix
        xsize, ysize, ang = g.size_pix # FWHM
        ellx, elly = func.drawellipse(g)
        blc = [int(min(ellx)), int(min(elly))]
        trc = [int(max(ellx)), int(max(elly))]

        blc[0] -= (x - blc[0]) * scale
        blc[1] -= (y - blc[1]) * scale
        trc[0] += (trc[0] - x) * scale
        trc[1] += (trc[1] - y) * scale
        # Format is: <id> <blcx> <blcy> <trcx> <trcy>
        # Note that we use gindx rather than g.gaus_num so that
        # all Gaussians will have a unique id, even if wavelet
        # Gaussians are included.
        outstr_list.append(str(gindx+1) + sep + str(blc[0]) + sep +
                           str(blc[1]) + sep + str(trc[0]) + sep +
                           str(trc[1]) +'\n')
    return outstr_list


def write_islands(img):
    import numpy as N
    import os

    ### write out island properties for reference since achaar doesnt work.
    filename = img.basedir + '/misc/'
    if not os.path.exists(filename): os.mkdir(filename)
    filename = filename + 'island_file'

    if img.j == 0:
      f = open(filename, 'w')
      f.write('Wavelet# Island_id  bbox origin shape mask_active mask_noisy size_active mean rms max_value ngaul gresid_mean '+\
              'gresid_rms resid_rms resid_mean nsource \n')
    else:
      f = open(filename, 'a')

    for isl in img.islands:
      f.write('%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %10i %10i %10i %.3e %.3e %.3e %5i %.3e %.3e %5i \n' \
              % (img.j, isl.island_id, isl.bbox[0].start, isl.bbox[0].stop, isl.bbox[1].start, isl.bbox[1].stop, \
              isl.origin[0], isl.origin[1], isl.shape[0], isl.shape[1], N.sum(~isl.mask_active), N.sum(~isl.mask_noisy), \
              isl.size_active, isl.mean, isl.rms, isl.max_value, len(isl.gaul), isl.gresid_mean, isl.gresid_rms, \
              len(isl.sources)))

    f.close()

def get_src(src_list, srcid):
    """Returns the source for srcid or None if not found"""
    for src in src_list:
        if src.source_id == srcid:
            return src
    return None

def list_and_sort_gaussians(img, patch=None, root=None,
                            sort_by='index'):
    """Returns sorted lists of Gaussians and their names and patch names.

    patch - can be "single", "gaussian", "source", or None

    Returns (outlist, outnames, patchnames)
    outlist is [[g1, g2, g3], [g4], ...]
    outnames is [['root_i2_s1_g1', 'root_i2_s1_g2', 'root_i2_s1_g3'], ...]
    patchnames is ['root_patch_s1', 'root_patch_s2', ...]

    The names are root_iXX_sXX_gXX (or wXX_iXX_sXX_gXX for wavelet Gaussians)
    """
    import numpy as N
    import functions as func

    # Define lists
    if root == None:
        root = img.parentname
    gauslist = []
    gausname = []
    outlist = []
    outnames = []
    patchnames = []
    patchnames_sorted = []
    gausflux = [] # fluxes of Gaussians
    gausindx = [] # indices of Gaussians
    patchflux = [] # total flux of each patch
    patchindx = [] # indices of sources

    # If a mask image is to be used to define patches, read it in and
    # make a rank image from it
#     if patch == 'mask':
#         patches_mask = func.readimage(mask_file)
#         act_pixels = patches_mask
#         rank = len(patches_mask.shape)
#         connectivity = nd.generate_binary_structure(rank, rank)
#         labels, count = nd.label(act_pixels, connectivity)
#         mask_labels = labels


    src_list = img.sources
    for src in src_list:
        for g in src.gaussians:
            gauslist.append(g)
            gausflux.append(g.total_flux)
            gausindx.append(g.gaus_num)
            jstr = '_w' + str(g.jlevel)
            gausname.append(root + jstr + '_i' + str(src.island_id) + '_s' +
                            str(src.source_id) + '_g' + str(g.gaus_num))
            if patch == 'gaussian':
                outlist.append(gauslist)
                outnames.append(gausname)
                patchnames.append(root + '_patch' + jstr + '_g' + str(g.gaus_num))
                patchflux.append(N.sum(gausflux))
                patchindx.append(g.gaus_num)
                gauslist = [] # reset for next Gaussian
                gausname = []
                gausflux = []
                gausindx = []
#                 if patch == 'mask':
#                     patchnum = mask_labels[g.centre_pix]

        if patch == 'source':
            sorted_gauslist = list(gauslist)
            sorted_gausname = list(gausname)
            if sort_by == 'flux':
                # Sort Gaussians by flux within each source
                indx = range(len(gausflux))
                indx.sort(lambda x,y: cmp(gausflux[x],gausflux[y]), reverse=True)
            elif sort_by == 'index':
                # Sort Gaussians by index within each source
                indx = range(len(gausindx))
                indx.sort(lambda x,y: cmp(gausindx[x],gausindx[y]), reverse=False)
            else:
                # Unrecognized property --> Don't sort
                indx = range(len(gausindx))
            for i, si in enumerate(indx):
                sorted_gauslist[i] = gauslist[si]
                sorted_gausname[i] = gausname[si]

            outlist.append(sorted_gauslist)
            outnames.append(sorted_gausname)
            patchnames.append(root + '_patch' + '_s' + str(src.source_id))
            patchflux.append(N.sum(gausflux))
            patchindx.append(src.source_id)
            gauslist = [] # reset for next source
            gausname = []
            gausflux = []

    # Sort
    if patch == 'single' or patch == None:
        outlist = [list(gauslist)]
        outlist_sorted = [list(gauslist)]
        outnames = [list(gausname)]
        outnames_sorted = [list(gausname)]
        if patch == 'single':
            patchnames = [root + '_patch']
        else:
            patchnames = [None]
        if sort_by == 'flux':
            # Sort by Gaussian flux
            indx = range(len(gauslist))
            indx.sort(lambda x,y: cmp(gausflux[x],gausflux[y]), reverse=True)
        elif sort_by == 'index':
            # Sort by Gaussian index
            indx = range(len(gausindx))
            indx.sort(lambda x,y: cmp(gausindx[x],gausindx[y]), reverse=False)
        else:
            # Unrecognized property --> Don't sort
            indx = range(len(gausindx))
        for i, si in enumerate(indx):
                outlist_sorted[0][i] = outlist[0][si]
                outnames_sorted[0][i] = outnames[0][si]
                patchnames_sorted = list(patchnames)
    else:
        outlist_sorted = list(outlist)
        outnames_sorted = list(outnames)
        patchnames_sorted = list(patchnames)
        if sort_by == 'flux':
            # Sort by patch flux
            indx = range(len(patchflux))
            indx.sort(lambda x,y: cmp(patchflux[x],patchflux[y]), reverse=True)
        elif sort_by == 'index':
            # Sort by source index
            indx = range(len(patchindx))
            indx.sort(lambda x,y: cmp(patchindx[x],patchindx[y]), reverse=False)
        else:
            # Unrecognized property --> Don't sort
            indx = range(len(gausindx))

        for i, si in enumerate(indx):
            outlist_sorted[i] = outlist[si]
            outnames_sorted[i] = outnames[si]
            patchnames_sorted[i] = patchnames[si]

    return (outlist_sorted, outnames_sorted, patchnames_sorted)

def make_output_columns(obj, fits=False, objtype='gaul', incl_spin=False,
                        incl_chan=False, incl_pol=False, incl_aper=False,
                        incl_empty=False, nmax=30, nchan=1):
    """Returns a list of column names, formats, and units for Gaussian, Source, or Shapelet"""
    import numpy as N

    # First, define a list of columns in order desired, using the names of
    # the attributes of the object
    if objtype == 'gaul':
        names = ['gaus_num', 'island_id', 'source_id', 'jlevel',
                 'centre_sky', 'centre_skyE', 'total_flux',
                 'total_fluxE', 'peak_flux', 'peak_fluxE',
                 'centre_pix', 'centre_pixE', 'size_sky', 'size_skyE',
                 'deconv_size_sky',
                 'deconv_size_skyE', 'total_flux_isl', 'total_flux_islE', 'rms',
                 'mean', 'gresid_rms', 'gresid_mean',
                 'code']
    elif objtype == 'srl':
        if incl_aper:
            infix = ['aperture_flux', 'aperture_fluxE']
        else:
            infix = []
        names = ['source_id', 'island_id', 'posn_sky_centroid',
                 'posn_sky_centroidE', 'total_flux',
                 'total_fluxE',
                 'peak_flux_max', 'peak_flux_maxE'] + infix + \
                 ['posn_sky_max', 'posn_sky_maxE',
                 'posn_pix_centroid', 'posn_pix_centroidE', 'posn_pix_max',
                 'posn_pix_maxE',
                 'size_sky', 'size_skyE', 'deconv_size_sky',
                 'deconv_size_skyE', 'total_flux_isl', 'total_flux_islE',
                 'rms_isl', 'mean_isl', 'gresid_rms',
                 'gresid_mean', 'code']
    elif objtype == 'shap':
        names = ['island_id', 'posn_sky_centroid',
                 'posn_sky_centroidE', 'total_flux',
                 'total_fluxE',
                 'peak_flux_max', 'peak_flux_maxE', 'posn_sky_max', 'posn_sky_maxE',
                 'posn_pix_centroid', 'posn_pix_centroidE', 'posn_pix_max',
                 'posn_pix_maxE', 'rms_isl', 'mean_isl', 'shapelet_basis' ,
                 'shapelet_beta', 'shapelet_nmax', 'shapelet_cf']
    else:
        print 'Object type unrecongnized.'
        return (None, None, None, None)
    if incl_spin:
        names += ['spec_indx', 'e_spec_indx']
    if incl_chan:
        names += ['specin_flux', 'specin_fluxE', 'specin_freq']
    if incl_pol:
        names += ['total_flux_Q', 'total_fluxE_Q', 'total_flux_U', 'total_fluxE_U',
                  'total_flux_V', 'total_fluxE_V', 'lpol_fraction', 'lpol_fraction_loerr',
                  'lpol_fraction_hierr', 'cpol_fraction', 'cpol_fraction_loerr',
                  'cpol_fraction_hierr', 'tpol_fraction',  'tpol_fraction_loerr',
                  'tpol_fraction_hierr', 'lpol_angle', 'lpol_angle_err']
    cnames = []
    cformats = []
    cunits = []
    cvals = []
    skip_next = False
    for n, name in enumerate(names):
        if hasattr(obj, name):
            if name in ['specin_flux', 'specin_fluxE', 'specin_freq']:
                # As these are variable length lists, they must
                # (unfortunately) be treated differently.
                    val = obj.__getattribute__(name)
                    colname = obj.__class__.__dict__[name]._colname
                    units = obj.__class__.__dict__[name]._units
                    for i in range(nchan):
                        if i < len(val):
                            cvals.append(val[i])
                            cnames.append(colname[0]+'_ch'+str(i+1))
                            cunits.append(units[0])
                        else:
                            cvals.append(N.NaN)
                            cnames.append(colname[0]+'_ch'+str(i+1))
                            cunits.append(units[0])
            else:
                if not skip_next:
                    val = obj.__getattribute__(name)
                    colname = obj.__class__.__dict__[name]._colname
                    units = obj.__class__.__dict__[name]._units
                    if units == None:
                        units = ' '
                    if isinstance(val, list):
                        # This is a list, so handle it differently. We assume the next
                        # entry will have the errors, and they are interleaved to be
                        # in the order (val, error).
                        next_name = names[n+1]
                        val_next = obj.__getattribute__(next_name)
                        colname_next = obj.__class__.__dict__[next_name]._colname
                        units_next = obj.__class__.__dict__[next_name]._units
                        if units_next == None:
                            units_next = ' '
                        for i in range(len(val)):
                            cvals.append(val[i])
                            cvals.append(val_next[i])
                            cnames.append(colname[i])
                            cnames.append(colname_next[i])
                            cunits.append(units[i])
                            cunits.append(units_next[i])
                        skip_next = True
                    elif isinstance(val, N.ndarray):
                        # This is a numpy array, so flatten it
                        tarr = val.flatten()
                        tarr2 = N.resize(tarr, nmax**2)
                        tarr2[tarr.shape[0]:] = N.NaN
                        cvals.append(tarr2)
                        cnames.append(colname)
                        cunits.append(units)
                    else:
                        cvals.append(val)
                        cnames.append(colname)
                        cunits.append(units)
                else:
                    skip_next = False

    for i, v in enumerate(cvals):
        if fits:
            if isinstance(v, int):
                cformats.append('1J')
            if isinstance(v, float):
                cformats.append('1D')
            if isinstance(v, str):
                cformats.append('1A')
            if isinstance(v, N.ndarray):
                cformats.append('%iD' % (nmax**2,))
        else:
            if isinstance(v, int):
                cformats.append('%4d')
            if isinstance(v, float):
                cformats.append('%10f')
            if isinstance(v, str):
                cformats.append('%4s')

    if objtype == 'gaul':
        if obj.gaus_num < 0 and not incl_empty:
            return (None, cnames, cformats, cunits)
    if objtype == 'srl':
        if obj.source_id < 0 and not incl_empty:
            return (None, cnames, cformats, cunits)
    return (cvals, cnames, cformats, cunits)
