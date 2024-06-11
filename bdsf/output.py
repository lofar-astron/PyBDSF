"""Module output.

Defines functions that write the results of source detection in a
variety of formats. These are then used as methods of Image objects
and/or are called by the outlist operation if output_all is True.
"""
from __future__ import print_function
from __future__ import absolute_import

from .image import Op


class Op_outlist(Op):
    """Write out list of Gaussians

    All available output lists are generated atm.
    """
    def __call__(self, img):
        if img.opts.output_all:
            import os
            if len(img.gaussians) > 0:
                dir = img.basedir + '/catalogues/'
                if not os.path.exists(dir):
                    os.makedirs(dir)
                self.write_bbs(img, dir)
                self.write_lsm(img, dir)
                self.write_gaul(img, dir)
                self.write_srl(img, dir)
                self.write_aips(img, dir)
                self.write_kvis(img, dir)
                self.write_ds9(img, dir, objtype='gaul')
                self.write_ds9(img, dir, objtype='srl')
                self.write_gaul_FITS(img, dir)
                self.write_srl_FITS(img, dir)
            if not os.path.exists(img.basedir + '/misc/'):
                os.makedirs(img.basedir + '/misc/')
            self.write_opts(img, img.basedir + '/misc/')
            self.save_opts(img, img.basedir + '/misc/')
            img.completed_Ops.append('outlist')

    def write_bbs(self, img, dir):
        """ Writes the gaussian list as a bbs-readable file"""
        if 'bbsname' in img.extraparams:
            name = img.extraparams['bbsname']
        else:
            name = img.imagename
        fname = dir + name + '.sky_in'

        # Write Gaussian list
        write_bbs_gaul(img, filename=fname, srcroot=img.opts.srcroot,
                       patch=img.opts.bbs_patches, sort_by='flux',
                       clobber=True, incl_empty=img.opts.incl_empty,
                       correct_proj=img.opts.correct_proj)

    def write_lsm(self, img, dir):
        """ Writes the gaussian list as an SAGECAL file"""
        fname = dir + img.imagename + '.lsm'
        write_lsm_gaul(img, filename=fname, sort_by='indx',
                       clobber=True,
                       incl_empty=img.opts.incl_empty)

    def write_gaul(self, img, dir):
        """ Writes the gaussian list as an ASCII file"""
        fname = dir + img.imagename + '.gaul'
        write_ascii_list(img, filename=fname, sort_by='indx',
                         clobber=True, objtype='gaul',
                         incl_empty=img.opts.incl_empty)

    def write_srl(self, img, dir):
        """ Writes the source list as an ASCII file"""
        fname = dir + img.imagename + '.srl'
        write_ascii_list(img, filename=fname, sort_by='indx',
                         clobber=True, objtype='srl',
                         incl_empty=img.opts.incl_empty)

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

    def write_ds9(self, img, dir, objtype='gaul'):
        """ Writes the gaussian list as a ds9 region file"""
        fname = dir + img.imagename + '.' + objtype + '.ds9.reg'
        write_ds9_list(img, filename=fname, srcroot=img.opts.srcroot,
                       clobber=True, deconvolve=False, objtype=objtype,
                       incl_empty=img.opts.incl_empty,)

    def write_gaul_FITS(self, img, dir):
        """ Writes the gaussian list as FITS binary table"""
        fname = dir + img.imagename+'.gaul.FITS'
        write_fits_list(img, filename=fname, sort_by='indx',
                        clobber=True, objtype='gaul',
                        incl_empty=img.opts.incl_empty,)

    def write_srl_FITS(self, img, dir):
        """ Writes the source list as FITS binary table"""
        fname = dir + img.imagename+'.srl.FITS'
        write_fits_list(img, filename=fname, sort_by='indx',
                        clobber=True, objtype='srl',
                        incl_empty=img.opts.incl_empty, incl_chan=img.opts.incl_chan)

    def write_shap_FITS(self, img, dir):
        """ Writes the shapelet list as a FITS file"""
        fname = dir + img.imagename + '.shap.FITS'
        write_fits_list(img, filename=fname, sort_by='indx',
                        clobber=True, objtype='shap')

    def write_opts(self, img, dir):
        """ Writes input parameters to a text file."""
        import inspect
        from . import mylogger

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
        fname = 'parameters_used'
        f = open(dir+fname, 'w')
        mylog.info('Writing '+dir+fname)
        for attr in inspect.getmembers(img.opts):
            if attr[0][0] != '_':
                if isinstance(attr[1], (int, str, bool, float, type(None), tuple, list)):
                    f.write('%-40s' % attr[0])
                    f.write(repr(attr[1])+'\n')

                    # Also print the values derived internally. They are all stored
                    # in img with the same name (e.g., img.opts.beam --> img.beam)
                    if hasattr(img, attr[0]):
                        used = img.__getattribute__(attr[0])
                        if used != attr[1] and isinstance(used, (int, str, bool, float,
                                                                 type(None), tuple,
                                                                 list)):
                            f.write('%-40s' % '    Value used')
                            f.write(repr(used)+'\n')
        f.close()

    def save_opts(self, img, dir):
        """ Saves input parameters to a PyBDSM save file."""
        from . import interface
        from . import mylogger

        mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
        fname = 'parameters.sav'
        mylog.info('Writing '+dir+fname)
        interface.save_pars(img, dir+fname, quiet=True)


def ra2hhmmss(deg):
    """Convert RA coordinate (in degrees) to HH MM SS"""

    from math import modf
    if deg < 0:
        deg += 360.0
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
    M = N.array([[0.9999256782, 0.0111820609, 0.00485794], [-0.0111820610, 0.9999374784, -0.0000271474],
                 [-0.0048579477, -0.0000271765, 0.9999881997]])

    r0 = N.zeros(3)
    r0[0] = cos(dec/rad) * cos(ra/rad)
    r0[1] = cos(dec/rad) * sin(ra/rad)
    r0[2] = sin(dec/rad)

    r0A = N.sum(r0*A)
    r1 = r0 - A + r0A*r0
    r = N.sum(M.transpose()*r1, axis=1)

    rscal = sqrt(N.sum(r*r))
    decj = asin(r[2]/rscal)*rad

    d1 = r[0] / rscal / cos(decj/rad)
    d2 = r[1] / rscal / cos(decj/rad)
    raj = acos(d1)*rad
    if d2 < 0.0:
        raj = 360.0 - raj

    Jcoord = [raj, decj]
    return Jcoord


def write_bbs_gaul(img, filename=None, srcroot=None, patch=None,
                   incl_primary=True, sort_by='flux',
                   clobber=False, incl_empty=False, correct_proj=True):
    """Writes Gaussian list to a BBS sky model"""
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM.write_gaul")
    if int(img.equinox) != 2000 and int(img.equinox) != 1950:
        mylog.warning('Equinox of input image is not J2000 or B1950. '
                      'Sky model may not be appropriate for BBS.')
    if int(img.equinox) == 1950:
        mylog.warning('Equinox of input image is B1950. Coordinates '
                      'will be precessed to J2000.')

    outl, outn, patl = list_and_sort_gaussians(img, patch=patch,
                                               root=srcroot, sort_by=sort_by)
    outstr_list = make_bbs_str(img, outl, outn, patl, incl_empty=incl_empty,
                               correct_proj=correct_proj)

    if filename is None:
        filename = img.imagename + '.sky_in'
    if os.path.exists(filename) and not clobber:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, 'w')
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_lsm_gaul(img, filename=None, srcroot=None, patch=None,
                   incl_primary=True, sort_by='flux',
                   clobber=False, incl_empty=False):
    """Writes Gaussian list to a SAGECAL lsm sky model"""
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM.write_gaul")
    if int(img.equinox) != 2000 and int(img.equinox) != 1950:
        mylog.warning('Equinox of input image is not J2000 or B1950. '
                      'Sky model may not be appropriate for Sagecal.')
    if int(img.equinox) == 1950:
        mylog.warning('Equinox of input image is B1950. Coordinates '
                      'will be precessed to J2000.')

    outl, outn, patl = list_and_sort_gaussians(img, patch=patch,
                                               root=srcroot, sort_by=sort_by)
    outstr_list = make_lsm_str(img, outl, outn, incl_empty=incl_empty)

    if filename is None:
        filename = img.imagename + '.lsm'
    if os.path.exists(filename) and not clobber:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, 'w')
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_ds9_list(img, filename=None, srcroot=None, deconvolve=False,
                   clobber=False, incl_empty=False, objtype='gaul'):
    """Writes Gaussian list to a ds9 region file"""
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if objtype == 'gaul':
        outl, outn, patl = list_and_sort_gaussians(img, patch=None)
    elif objtype == 'srl':
        root = img.parentname
        outl = [img.sources]
        if incl_empty:
            # Append the dummy sources for islands without any unflagged Gaussians
            outl[0] += img.dsources
        outn = []
        for src in img.sources:
            outn.append(root + '_i' + str(src.island_id) + '_s' +
                        str(src.source_id))
        if incl_empty:
            # Append the dummy sources for islands without any unflagged Gaussians
            for dsrc in img.dsources:
                outn.append(root + '_i' + str(dsrc.island_id) + '_s' +
                            str(dsrc.source_id))
        outn = [outn]
    outstr_list = make_ds9_str(img, outl, outn, deconvolve=deconvolve,
                               objtype=objtype, incl_empty=incl_empty)
    if filename is None:
        filename = img.imagename + '.' + objtype + '.reg'
    if os.path.exists(filename) and not clobber:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, "w")
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_ascii_list(img, filename=None, sort_by='indx', format='ascii',
                     incl_chan=False, incl_empty=False, clobber=False, objtype='gaul'):
    """Writes Gaussian list to an ASCII file"""
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if objtype == 'gaul':
        outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)
    elif objtype == 'srl':
        outl = [img.sources]
        if incl_empty:
            # Append the dummy sources for islands without any unflagged Gaussians
            outl[0] += img.dsources
    outstr_list = make_ascii_str(img, outl, objtype=objtype, incl_chan=incl_chan,
                                 incl_empty=incl_empty, format=format)
    if filename is None:
        if objtype == 'gaul':
            filename = img.imagename + '.gaul'
        elif objtype == 'srl':
            filename = img.imagename + '.srl'
    if os.path.exists(filename) and not clobber:
        return None
    mylog.info('Writing ' + filename)
    f = open(filename, "w")
    for s in outstr_list:
        f.write(s)
    f.close()
    return filename


def write_casa_gaul(img, filename=None,  incl_empty=False, clobber=False):
    """Writes a clean box file for use in casapy"""
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    outl, outn, patl = list_and_sort_gaussians(img, patch=None)
    outstr_list = make_casa_str(img, outl)
    if filename is None:
        filename = img.imagename + '.box'
    if os.path.exists(filename) and not clobber:
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
    from . import mylogger
    import os
    import numpy as N
    from astropy.io import fits as pyfits
    from ._version import __version__

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if objtype == 'gaul':
        outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)
    elif objtype == 'srl':
        outl = [img.sources]
        if incl_empty:
            # Append the dummy sources for islands without any unflagged Gaussians
            outl[0] += img.dsources
    elif objtype == 'shap':
        outl = [[isl for isl in img.islands if hasattr(isl, 'shapelet_nmax')]]

    nmax = 0
    if objtype == 'shap':
        # loop over shapelets and get maximum size of coefficient matrix
        for isl in outl[0]:
            if hasattr(isl, 'shapelet_nmax'):
                if isl.shapelet_nmax > nmax:
                    nmax = isl.shapelet_nmax
        nmax += 1

    if img.opts.aperture is not None:
        incl_aper = True
    else:
        incl_aper = False
    if len(outl[0]) > 0:
        cvals, cnames, cformats, cunits = make_output_columns(outl[0][0], fits=True,
                                                              objtype=objtype,
                                                              incl_spin=img.opts.spectralindex_do,
                                                              incl_chan=incl_chan,
                                                              incl_pol=img.opts.polarisation_do,
                                                              incl_aper=incl_aper,
                                                              incl_empty=incl_empty,
                                                              nmax=nmax, nchan=img.nchan)
    out_list = make_fits_list(img, outl, objtype=objtype, nmax=nmax, incl_empty=incl_empty, incl_chan=incl_chan)
    col_list = []
    for ind, col in enumerate(out_list):
        list1 = pyfits.Column(name=cnames[ind], format=cformats[ind],
                              unit=cunits[ind], array=N.array(out_list[ind]))
        col_list.append(list1)
    if len(col_list) == 0:
        col_list = [pyfits.Column(name='Blank', format='1J')]

    tbhdu = pyfits.BinTableHDU.from_columns(col_list)

    if objtype == 'gaul':
        tbhdu.header.add_comment('Gaussian list for '+img.filename)
    elif objtype == 'srl':
        tbhdu.header.add_comment('Source list for '+img.filename)
    elif objtype == 'shap':
        tbhdu.header.add_comment('Shapelet list for '+img.filename)
    tbhdu.header.add_comment('Generated by PyBDSM version %s'
                             % (__version__, ))
    freq = "%.5e" % img.frequency
    tbhdu.header.add_comment('Reference frequency of the detection ("ch0") image: %s Hz' % freq)
    tbhdu.header.add_comment('Equinox : %s' % img.equinox)
    tbhdu.header['INIMAGE'] = (img.filename, 'Filename of image')
    tbhdu.header['FREQ0'] = (float(freq), 'Reference frequency')
    tbhdu.header['EQUINOX'] = (img.equinox, 'Equinox')

    for key in img.header.keys():
        if key in ['HISTORY', 'COMMENT', '']:
            continue
        tbhdu.header.add_comment('%s = %s' % (key, repr(img.header[key])))

    if filename is None:
        filename = img.imagename + '.' + objtype + '.fits'
    if os.path.exists(filename) and not clobber:
        return None
    mylog.info('Writing ' + filename)
    try:
        tbhdu.writeto(filename, overwrite=True)
    except TypeError:
        # The "overwrite" argument was added in astropy v1.3, so fall back to "clobber"
        # if it doesn't work
        tbhdu.writeto(filename, clobber=True)
    return filename


def write_kvis_ann(img, filename=None, sort_by='indx',
                   clobber=False):
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if filename is None:
        filename = img.imagename + '.kvis.ann'
    if os.path.exists(filename) and not clobber:
        return None
    f = open(filename, 'w')
    mylog.info('Writing '+filename)
    f.write("### KVis annotation file\n\n")
    f.write("color green\n\n")

    outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)
    for g in outl[0]:
        iidx = g.island_id
        # kvis does not correct for postion-dependent angle or pixel scale
        # for region files, so we must use the uncorrected values
        ra, dec = g.centre_sky
        shape = g.size_sky_uncorr

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
    from .output import ra2hhmmss, dec2ddmmss
    from . import mylogger
    import os

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"Output")
    if filename is None:
        filename = img.imagename + '.star'
    if os.path.exists(filename) and not clobber:
        return None
    f = open(filename, 'w')
    mylog.info('Writing '+filename)

    outl, outn, patl = list_and_sort_gaussians(img, patch=None, sort_by=sort_by)

    for g in outl[0]:
        A = g.peak_flux
        ra, dec = g.centre_sky
        shape = g.size_sky_uncorr
        # convert to canonical representation
        ra = ra2hhmmss(ra)
        dec = dec2ddmmss(dec)
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


def make_bbs_str(img, glist, gnames, patchnames, objtype='gaul',
                 incl_empty=False, correct_proj=True):
    """Makes a list of string entries for a BBS sky model."""
    from .output import ra2hhmmss
    from .output import dec2ddmmss
    import numpy as N

    outstr_list = []
    freq = "%.5e" % img.frequency

    if len(patchnames) == 0:
        # Handle empty list: just write default header
        outstr_list.append("format = Name, Type, Ra, Dec, I, Q, U, V, "
                           "MajorAxis, MinorAxis, Orientation, "
                           "ReferenceFrequency='"+freq+"', "
                           "SpectralIndex='[]'\n\n")
    elif patchnames[0] is None:
        outstr_list.append("format = Name, Type, Ra, Dec, I, Q, U, V, "
                           "MajorAxis, MinorAxis, Orientation, "
                           "ReferenceFrequency='"+freq+"', "
                           "SpectralIndex='[]'\n\n")
    else:
        outstr_list.append("format = Name, Type, Patch, Ra, Dec, I, Q, U, V, "
                           "MajorAxis, MinorAxis, Orientation, "
                           "ReferenceFrequency='"+freq+"', "
                           "SpectralIndex='[]'\n\n")
    if objtype == 'shap':
        raise RuntimeError("Shapelets not yet supported in the BBS format.")
    else:
        patchname_last = ''
        for pindx, patch_name in enumerate(patchnames):  # loop over patches
            if patch_name is not None and patch_name != patchname_last:
                outstr_list.append(', , ' + patch_name + ', 00:00:00, +00.00.00\n')
                patchname_last = patch_name
            gaussians_in_patch = glist[pindx]
            names_in_patch = gnames[pindx]
            for gindx, g in enumerate(gaussians_in_patch):
                if g.gaus_num >= 0 or (g.gaus_num < 0 and incl_empty):
                    src_name = names_in_patch[gindx]
                    ra, dec = g.centre_sky
                    if img.equinox == 1950:
                        ra, dec = B1950toJ2000([ra, dec])
                    ra = ra2hhmmss(ra)
                    sra = str(ra[0]).zfill(2)+':'+str(ra[1]).zfill(2)+':'+str("%.6f" % (ra[2])).zfill(6)
                    dec = dec2ddmmss(dec)
                    decsign = ('-' if dec[3] < 0 else '+')
                    sdec = decsign+str(dec[0]).zfill(2)+'.'+str(dec[1]).zfill(2)+'.'+str("%.6f" % (dec[2])).zfill(6)
                    total = str("%.3e" % (g.total_flux))
                    if correct_proj:
                        deconv = list(g.deconv_size_sky)
                    else:
                        deconv = list(g.deconv_size_sky_uncorr)
                    if deconv[0] == 0.0 and deconv[1] == 0.0:
                        stype = 'POINT'
                        deconv[2] = 0.0
                    else:
                        stype = 'GAUSSIAN'
                    deconv1 = str("%.5e" % (deconv[0]*3600.0))
                    deconv2 = str("%.5e" % (deconv[1]*3600.0))
                    deconv3 = str("%.5e" % (deconv[2]))
                    deconvstr = deconv1 + ', ' + deconv2 + ', ' + deconv3
                    specin = '-0.8'
                    if 'spectralindex' in img.completed_Ops:
                        if g.spec_indx is not None and N.isfinite(g.spec_indx):
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
                    if patch_name is None:
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


def make_lsm_str(img, glist, gnames, incl_empty=False):
    """Makes a list of string entries for a SAGECAL sky model."""
    from .output import ra2hhmmss
    from .output import dec2ddmmss
    import numpy as N
    from ._version import __version__

    outstr_list = ["# SAGECAL sky model\n"]
    freq = "%.5e" % img.frequency
    outstr_list.append('# Generated by PyBDSM version %s\n'
                       % (__version__, ))
    outstr_list.append("# Name  | RA (hr,min,sec) | DEC (deg,min,sec) | I | Q | U | V | SI | RM | eX | eY | eP | freq0\n\n")
    for gindx, g in enumerate(glist[0]):
        if g.gaus_num >= 0 or (g.gaus_num < 0 and incl_empty):
            src_name = gnames[0][gindx]
            ra, dec = g.centre_sky
            if img.equinox == 1950:
                ra, dec = B1950toJ2000([ra, dec])
            ra = ra2hhmmss(ra)
            sra = str(ra[0]).zfill(2)+' '+str(ra[1]).zfill(2)+' '+str("%.6f" % (ra[2])).zfill(6)
            dec = dec2ddmmss(dec)
            decsign = ('-' if dec[3] < 0 else '+')
            sdec = decsign+str(dec[0]).zfill(2)+' '+str(dec[1]).zfill(2)+' '+str("%.6f" % (dec[2])).zfill(6)
            total = str("%.3e" % (g.total_flux))
            deconv = list(g.deconv_size_sky)
            if deconv[0] == 0.0 and deconv[1] == 0.0:
                sname = 'P' + src_name
                deconv[2] = 0.0
            else:
                sname = 'G' + src_name
                # Make sure Gaussian is not 1-D, as SAGECAL cannot handle these
                if deconv[0] < 1e-5:
                    deconv[0] = 1e-5
                if deconv[1] < 1e-5:
                    deconv[1] = 1e-5
            # The following conversions taken from the SABECAL script "convert_skymodel.py"
            deconv1 = str("%.5e" % (deconv[0]*N.pi/180.0/2.0))
            deconv2 = str("%.5e" % (deconv[1]*N.pi/180.0/2.0))
            deconv3 = str("%.5e" % (N.pi/2-(N.pi-deconv[2]/180.0*N.pi)))
            deconvstr = deconv1 + ' ' + deconv2 + ' ' + deconv3
            specin = '-0.8'
            if 'spectralindex' in img.completed_Ops:
                if g.spec_indx is not None and N.isfinite(g.spec_indx):
                    specin = str("%.3e" % (g.spec_indx))
            sep = ' '
            if img.opts.polarisation_do:
                Q_flux = str("%.3e" % g.total_flux_Q)
                U_flux = str("%.3e" % g.total_flux_U)
                V_flux = str("%.3e" % g.total_flux_V)
            else:
                Q_flux = '0.0'
                U_flux = '0.0'
                V_flux = '0.0'
            outstr_list.append(sname + sep + sra + sep +
                               sdec + sep + total + sep + Q_flux + sep +
                               U_flux + sep + V_flux + sep +
                               specin + sep + '0' + sep + deconvstr + sep +
                               freq + sep + '\n')
    return outstr_list


def make_ds9_str(img, glist, gnames, deconvolve=False, objtype='gaul', incl_empty=False):
    """Makes a list of string entries for a ds9 region file."""
    from . import mylogger

    outstr_list = []
    if img.equinox is None:
        equinox = 'fk5'
    else:
        if int(img.equinox) == 2000:
            equinox = 'fk5'
        elif int(img.equinox) == 1950:
            equinox = 'fk4'
        else:
            mylog = mylogger.logging.getLogger("PyBDSM.write_ds9")
            mylog.warning('Equinox of input image is not J2000 or B1950. '
                          'Regions may not be correct.')
            equinox = 'fk5'

    outstr_list.append('# Region file format: DS9 version 4.0\nglobal color=green '
                       'font="helvetica 10 normal" select=1 highlite=1 edit=1 '
                       'move=1 delete=1 include=1 fixed=0 source\n'+equinox+'\n')

    for gindx, g in enumerate(glist[0]):
        if objtype == 'gaul':
            objid = g.gaus_num
        else:
            objid = g.source_id
        if objid >= 0 or (objid < 0 and incl_empty):
            src_name = gnames[0][gindx]
            if objtype == 'gaul':
                ra, dec = g.centre_sky
            else:
                ra, dec = g.posn_sky_centroid

            # ds9 does not correct for postion-dependent angle or pixel scale
            # for region files, so we must use the uncorrected values
            if deconvolve:
                deconv = g.deconv_size_sky_uncorr
            else:
                deconv = g.size_sky_uncorr
            if deconv[0] == 0.0 and deconv[1] == 0.0:
                deconv[2] = 0.0
                region = 'point(' + str(ra) + ',' + str(dec) + \
                    ') # point=cross width=2 text={' + src_name + '}\n'
            else:
                # ds9 can't handle 1-D Gaussians, so make sure they are 2-D
                if deconv[0] < 1.0/3600.0:
                    deconv[0] = 1.0/3600.0
                if deconv[1] < 1.0/3600.0:
                    deconv[1] = 1.0/3600.0
                region = 'ellipse(' + str(ra) + ',' + str(dec) + ',' + \
                    str(deconv[0]*3600.0) + '",' + str(deconv[1]*3600.0) + \
                    '",' + str(deconv[2]+90.0) + ') # text={' + src_name + '}\n'
            outstr_list.append(region)
    return outstr_list


def make_ascii_str(img, glist, objtype='gaul', format='ascii', incl_empty=False,
                   incl_chan=False):
    """Makes a list of string entries for an ascii region file."""
    from ._version import __version__
    outstr_list = []
    freq = "%.5e" % img.frequency

    if objtype == 'gaul':
        outstr_list.append('# Gaussian list for '+img.filename+'\n')
    elif objtype == 'srl':
        outstr_list.append('# Source list for '+img.filename+'\n')
    outstr_list.append('# Generated by PyBDSM version %s\n'
                       % (__version__, ))
    outstr_list.append('# Reference frequency of the detection ("ch0") image: %s Hz\n' % freq)
    outstr_list.append('# Equinox : %s \n\n' % img.equinox)
    if img.opts.aperture is not None:
        incl_aper = True
    else:
        incl_aper = False

    for i, g in enumerate(glist[0]):
        cvals, cnames, cformats, cunits = make_output_columns(g, fits=False,
                                                              objtype=objtype,
                                                              incl_spin=img.opts.spectralindex_do,
                                                              incl_chan=incl_chan,
                                                              incl_pol=img.opts.polarisation_do,
                                                              incl_aper=incl_aper,
                                                              incl_empty=incl_empty,
                                                              nchan=img.nchan)
        if cvals is not None:
            cformats[-1] += "\n"
            if format == 'ascii':
                if i == 0:
                    outstr_list.append("# " + " ".join(cnames) + "\n")
                outstr_list.append(" ".join(cformats).format(*cvals))
            else:
                if i == 0:
                    outstr_list.append("# " + ", ".join(cnames) + "\n")
                outstr_list.append(", ".join(cformats).format(*cvals))
    return outstr_list


def make_fits_list(img, glist, objtype='gaul', nmax=30, incl_empty=False,
                   incl_chan=False):
    from . import functions as func

    out_list = []
    if img.opts.aperture is not None:
        incl_aper = True
    else:
        incl_aper = False
    for g in glist[0]:
        cvals, ext1, ext2, ext3 = make_output_columns(g, fits=True, objtype=objtype,
                                                      incl_spin=img.opts.spectralindex_do,
                                                      incl_chan=incl_chan,
                                                      incl_pol=img.opts.polarisation_do,
                                                      incl_aper=incl_aper,
                                                      incl_empty=incl_empty,
                                                      nmax=nmax, nchan=img.nchan)
        if cvals is not None:
            out_list.append(cvals)
    out_list = func.trans_gaul(out_list)
    return out_list


def make_casa_str(img, glist):
    """Makes a list of string entries for a casa region file."""
    from . import functions as func
    outstr_list = ['#CRTFv0 CASA Region Text Format version 0\n']
    scale = 2.0  # scale box to 2 times FWHM of Gaussian
    for gindx, g in enumerate(glist[0]):
        x, y = g.centre_pix
        ellx, elly = func.drawellipse(g)
        blc = [min(ellx), min(elly)]
        trc = [max(ellx), max(elly)]

        blc[0] -= (x - blc[0]) * scale
        blc[1] -= (y - blc[1]) * scale
        trc[0] += (trc[0] - x) * scale
        trc[1] += (trc[1] - y) * scale

        blc_sky = img.pix2sky(blc)
        trc_sky = img.pix2sky(trc)

        blc_sky_str = convert_radec_str(blc_sky[0], blc_sky[1])
        trc_sky_str = convert_radec_str(trc_sky[0], trc_sky[1])

        # Format is: box [ [<blcx>, <blcy>], [<trcx>, <trcy>] ]
        # Note that we use gindx rather than g.gaus_num so that
        # all Gaussians will have a unique id, even if wavelet
        # Gaussians are included.
        outstr_list.append('box [[' + ', '.join(blc_sky_str) + '], [' +
                           ', '.join(trc_sky_str) + ']] coord=J2000\n')
    return outstr_list


def write_islands(img):
    import numpy as N
    import os

    # write out island properties for reference since achaar doesnt work.
    filename = img.basedir + '/misc/'
    if not os.path.exists(filename):
        os.makedirs(filename)
    filename = filename + 'island_file'

    if img.j == 0:
        f = open(filename, 'w')
        f.write('Wavelet# Island_id  bbox origin shape mask_active mask_noisy size_active mean rms max_value ngaul gresid_mean ' +
                'gresid_rms resid_rms resid_mean nsource \n')
    else:
        f = open(filename, 'a')

    for isl in img.islands:
        f.write('%5i %5i %5i %5i %5i %5i %5i %5i %5i %5i %10i %10i %10i %.3e %.3e %.3e %5i %.3e %.3e %5i \n'
                % (img.j, isl.island_id, isl.bbox[0].start, isl.bbox[0].stop, isl.bbox[1].start, isl.bbox[1].stop,
                   isl.origin[0], isl.origin[1], isl.shape[0], isl.shape[1], N.sum(~isl.mask_active), N.sum(~isl.mask_noisy),
                   isl.size_active, isl.mean, isl.rms, isl.max_value, len(isl.gaul), isl.gresid_mean, isl.gresid_rms,
                   len(isl.sources)))

    f.close()


def get_src(src_list, srcid):
    """Returns the source for srcid or None if not found"""
    for src in src_list:
        if src.source_id == srcid:
            return src
    return None


def convert_radec_str(ra, dec):
    """Takes ra, dec in degrees and returns BBS/CASA strings"""
    ra = ra2hhmmss(ra)
    sra = str(ra[0]).zfill(2)+':'+str(ra[1]).zfill(2)+':'+str("%.3f" % (ra[2])).zfill(6)
    dec = dec2ddmmss(dec)
    decsign = ('-' if dec[3] < 0 else '+')
    sdec = decsign+str(dec[0]).zfill(2)+'.'+str(dec[1]).zfill(2)+'.'+str("%.3f" % (dec[2])).zfill(6)
    return sra, sdec


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
    from . import functions as func

    # Define lists
    if root is None:
        root = img.parentname
    gauslist = []
    gausname = []
    outlist = []
    outnames = []
    patchnames = []
    patchnames_sorted = []
    gausflux = []  # fluxes of Gaussians
    gausindx = []  # indices of Gaussians
    patchflux = []  # total flux of each patch
    patchindx = []  # indices of sources
    patchnums = []  # number of patch from mask

    # If a mask image is to be used to define patches, read it in and
    # make a rank image from it
    use_mask = False
    if patch not in ['single', 'gaussian', 'source', None]:
        mask_file = img.opts.bbs_patches_mask
        patches_mask, hdr = func.read_image_from_file(mask_file, img, img.indir)
        use_mask = True
        act_pixels = patches_mask[0, 0]
        rank = len(act_pixels.shape)
        import scipy.ndimage as nd
        connectivity = nd.generate_binary_structure(rank, rank)
        mask_labels, count = nd.label(act_pixels, connectivity)

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
                gauslist = []  # reset for next Gaussian
                gausname = []
                gausflux = []
                gausindx = []
            if use_mask:
                patchnums.append(mask_labels[g.centre_pix[0], g.centre_pix[1]])

        if patch == 'source':
            sorted_gauslist = list(gauslist)
            sorted_gausname = list(gausname)
            if sort_by == 'flux':
                # Sort Gaussians by flux within each source
                indx = N.argsort(N.array(gausflux)).tolist()
                indx.reverse()
            elif sort_by == 'index':
                # Sort Gaussians by index within each source
                indx = N.argsort(N.array(gausindx)).tolist()
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
            gauslist = []  # reset for next source
            gausname = []
            gausflux = []

    if use_mask:
        unique_patch_ids = set(patchnums)

        # Check if there is a patch with id = 0. If so, this means there were
        # some Gaussians that fell outside of the regions in the patch
        # mask file.
        if 0 in unique_patch_ids:
            from . import mylogger
            mylog = mylogger.logging.getLogger("PyBDSM.write_gaul")
            mylog.warning('Some sources fall outside of the regions '
                          'defined in the mask file. These sources are not '
                          'included in the output sky model.')
        for p in unique_patch_ids:
            if p != 0:
                in_patch = N.where(patchnums == p)
                outlist.append(N.array(gauslist)[in_patch].tolist())
                outnames.append(N.array(gausname)[in_patch].tolist())
                patchnames.append('patch_'+str(p))
                patchflux.append(N.sum(N.array(gausflux)[in_patch]))
                patchindx.append(p)

    # Sort
    if patch == 'single' or patch is None:
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
            indx = N.argsort(N.array(gausflux)).tolist()
            indx.reverse()
        elif sort_by == 'index':
            # Sort by Gaussian index
            indx = N.argsort(N.array(gausindx)).tolist()
        else:
            # Unrecognized property --> Don't sort
            indx = list(range(len(gausindx)))
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
            indx = N.argsort(N.array(patchflux)).tolist()
            indx.reverse()
        elif sort_by == 'index':
            # Sort by source index
            indx = N.argsort(N.array(patchindx)).tolist()
        else:
            # Unrecognized property --> Don't sort
            indx = list(range(len(gausindx)))

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
                 'size_sky_uncorr', 'size_skyE_uncorr',
                 'deconv_size_sky', 'deconv_size_skyE',
                 'deconv_size_sky_uncorr', 'deconv_size_skyE_uncorr',
                 'total_flux_isl', 'total_flux_islE', 'rms',
                 'mean', 'gresid_rms', 'gresid_mean', 'wave_rms', 'wave_mean',
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
                 'size_sky', 'size_skyE',
                 'size_sky_uncorr', 'size_skyE_uncorr',
                 'deconv_size_sky', 'deconv_size_skyE',
                 'deconv_size_sky_uncorr', 'deconv_size_skyE_uncorr',
                 'total_flux_isl', 'total_flux_islE',
                 'rms_isl', 'mean_isl', 'gresid_rms',
                 'gresid_mean', 'code']
    elif objtype == 'shap':
        names = ['island_id', 'shapelet_posn_sky', 'shapelet_posn_skyE',
                 'shapelet_basis', 'shapelet_beta', 'shapelet_nmax', 'shapelet_cf']
    else:
        print('Object type unrecongnized.')
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
                colname = obj.__dict__[name+'_def']._colname
                units = obj.__dict__[name+'_def']._units
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
                    colname = obj.__dict__[name+'_def']._colname
                    units = obj.__dict__[name+'_def']._units
                    if units is None:
                        units = ' '
                    if isinstance(val, list) or isinstance(val, tuple):
                        # This is a list, so handle it differently. We assume the next
                        # entry will have the errors, and they are interleaved to be
                        # in the order (val, error).
                        next_name = names[n+1]
                        val_next = obj.__getattribute__(next_name)
                        colname_next = obj.__dict__[next_name+'_def']._colname
                        units_next = obj.__dict__[next_name+'_def']._units
                        if units_next is None:
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
                cformats.append('J')
            elif isinstance(v, float) or isinstance(v, N.float32) or isinstance(v, N.float64):
                cformats.append('D')
            elif isinstance(v, str):
                cformats.append('A')
            elif isinstance(v, N.ndarray):
                cformats.append('%iD' % (nmax**2,))
            else:
                raise RuntimeError("Format not supported.")
        else:
            if isinstance(v, int):
                cformats.append('{'+str(i)+':4d}')
            elif isinstance(v, float) or isinstance(v, N.float32) or isinstance(v, N.float64):
                cformats.append('{'+str(i)+':.14f}')
            elif isinstance(v, str):
                cformats.append('{'+str(i)+':4s}')
            else:
                raise RuntimeError("Format not supported.")

    if objtype == 'gaul':
        if obj.gaus_num < 0 and not incl_empty:
            return (None, cnames, cformats, cunits)
    if objtype == 'srl':
        if obj.source_id < 0 and not incl_empty:
            return (None, cnames, cformats, cunits)
    return (cvals, cnames, cformats, cunits)
