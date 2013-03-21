"""Module readimage.

Defines operation Op_readimage which initializes image and WCS

The current implementation tries to reduce input file to 2D if
possible, as this makes more sence atm. One more important thing
to note -- in its default configuration pyfits will read data
in non-native format, so we have to convert it before usage. See
the read_image_from_file in functions.py for details.

Lastly, wcs and spectal information are stored in the PyWCS
object img.wcs_obj.
"""

import numpy as N
from image import *
from functions import read_image_from_file
import mylogger
import sys

Image.imagename = String(doc="Identifier name for output files")
Image.filename = String(doc="Name of input file without extension")
Image.bbspatchnum = Int(doc="To keep track of patch number for bbs file "\
                            "for seperate patches per source")
Image.frequency = Float(doc="Frequency in the header")
Image.use_io = String(doc="pyfits or pyrap")
Image.j = Int(doc="Wavelet order j, 0 for normal run")
Image.freq_pars = Tuple((0.0, 0.0, 0.0),
                        doc="Frequency prarmeters from the header: (crval, cdelt, crpix)")
Image.waveletimage = Bool(doc="Whether a wavelet transform image of not")
Image.equinox = Float(2000.0, doc='Equinox of input image from header')

class Op_readimage(Op):
    """Image file loader

    Loads image and configures wcslib machinery for it.
    """
    def __call__(self, img):
        import time, os
        mylog = mylogger.logging.getLogger("PyBDSM." + img.log + "Readimage")

        if img.opts.filename == '':
            raise RuntimeError('Image file name not specified.')

        # Check for trailing "/" in file name (since CASA images are directories).
        # Although the general rule is to not alter the values in opts (only the
        # user should be able to alter these), in this case there is no harm in
        # replacing the file name in opts with the '/' trimmed off.
        if img.opts.filename[-1] == '/':
            img.opts.filename = img.opts.filename[:-1]

        # Determine indir if not explicitly given by user (in img.opts.indir)
        if img.opts.indir == None:
            indir = os.path.dirname(img.opts.filename)
            if indir == '':
                indir = './'
            img.indir = indir
        else:
            img.indir = img.opts.indir

        image_file = os.path.basename(img.opts.filename)
        result = read_image_from_file(image_file, img, img.indir)
        if result == None:
            raise RuntimeError("Cannot open file " + repr(image_file) + ". " + img._reason)
        else:
            data, hdr = result

        # Store data and header in img. If polarisation_do = False, only store pol == 'I'
        img.nchan = data.shape[1]
        img.nstokes = data.shape[0]
        mylogger.userinfo(mylog, 'Image size',
                          str(data.shape[-2:]) + ' pixels')
        mylogger.userinfo(mylog, 'Number of channels',
                          '%i' % data.shape[1])
        mylogger.userinfo(mylog, 'Number of Stokes parameters',
                          '%i' % data.shape[0])
        if img.opts.polarisation_do and data.shape[0] == 1:
            img.opts.polarisation_do = False
            mylog.warning('Image has Stokes I only. Polarisation module disabled.')
        if img.opts.polarisation_do or data.shape[0] == 1:
            img.image = data
        else:
            img.image = data[0, :].reshape(1, data.shape[1], data.shape[2], data.shape[3])
        img.header = hdr
        img.j = 0

        ### initialize wcs conversion routines
        self.init_wcs(img)
        self.init_beam(img)
        self.init_freq(img)
        year, code = self.get_equinox(img)
        if year == None:
            mylog.info('Equinox not found in image header. Assuming J2000.')
            img.equinox = 2000.0
        else:
            mylog.info('Equinox of image is %f.' % year)
            img.equinox = year

        # Try to trim common extensions from filename
        root, ext = os.path.splitext(img.opts.filename)
        if ext in ['.fits', '.FITS', '.image']:
            fname = root
        elif ext in ['.gz', '.GZ']:
            root2, ext2 = os.path.splitext(root)
            if ext2 in ['.fits', '.FITS', '.image']:
                fname = root2
            else:
                fname = root
        else:
            fname = img.opts.filename
        img.filename = img.opts.filename
        img.parentname = fname
        img.imagename = fname + '.pybdsm'
        img.basedir = './' + fname + '_pybdsm/'
        if img.opts.output_all:
            # Set up directory to write output to
            opdir = img.opts.opdir_overwrite
            if opdir not in ['overwrite', 'append']:
                img.opts.opdir_overwrite = 'append'
            if opdir == 'append':
                mylog.info('Appending output files to directory ' + img.basedir)
            else:
                mylog.info('Overwriting output files (if any) in directory ' + img.basedir)
                if os.path.isdir(img.basedir):
                    os.system("rm -fr " + img.basedir + '/*')
            if not os.path.isdir(img.basedir):
                os.makedirs(img.basedir)

            # Now add solname (if any) and time to basedir
            if img.opts.solnname != None:
                img.basedir += img.opts.solnname + '_'
            img.basedir += time.strftime("%d%b%Y_%H.%M.%S")

            # Make the final output directory
            if not os.path.isdir(img.basedir):
                os.makedirs(img.basedir)

        # Check for zeros and blank if img.opts.blank_zeros is True
        if img.opts.blank_zeros:
            zero_pixels = N.where(img.image[0] == 0.0)
            mylog.info('Blanking %i zeros in image' % len(zero_pixels[1]))
            img.image[0][zero_pixels] = N.nan

        img.completed_Ops.append('readimage')
        return img

    def init_wcs(self, img):
        """Initialize wcs pixel <=> sky conversion routines.
        """
        from math import pi
        import warnings

        hdr = img.header

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            from pywcs import WCS
            t = WCS(hdr)
            t.wcs.fix()

        acdelt = [abs(hdr['cdelt1']), abs(hdr['cdelt2'])]

        # Here we define p2s and s2p to allow celestial coordinate
        # transformations. Transformations for other axes (e.g.,
        # spectral) are striped out.
        def p2s(self, xy):
            xy = list(xy)
            for i in range(t.wcs.naxis - 2):
                xy.append(0)
            xy_arr = N.array([xy])
            sky = self.wcs_pix2sky(xy_arr, 0)#, ra_dec_order=True)
            return sky.tolist()[0][0:2]

        def s2p(self, rd):
            rd = list(rd)
            for i in range(t.wcs.naxis - 2):
                rd.append(0)
            rd_arr = N.array([rd])
            pix = self.wcs_sky2pix(rd_arr, 0)#, ra_dec_order=True)
            return pix.tolist()[0][0:2]

        # Here we define functions to transform Gaussian parameters (major axis,
        # minor axis, pos. angle) from the image plane to the celestial sphere.
        # These transforms are valid only at the Gaussian's center and ignore
        # any change across the extent of the Gaussian.
        def gaus2pix(x, location=None, use_wcs=True):
            """ Converts Gaussian parameters in deg to pixels.

            x - (maj [deg], min [deg], pa [deg])
            location - specifies the location in pixels (x, y) for which
                transform is desired
            Input beam angle should be degrees CCW from North.
            The output beam angle is degrees CCW from the +y axis of the image.
            """
            if use_wcs:
                bmaj, bmin, bpa = x
                brot = self.get_rot(img, location) # rotation delta CCW (in degrees) between N and +y axis of image

                s1 = self.angdist2pixdist(img, bmaj, bpa, location=location)
                s2 = self.angdist2pixdist(img, bmin, bpa + 90.0, location=location)
                th = bpa + brot
                if s1 < s2:
                    s1, s2 = s2, s1
                    th += 90.0
                th = divmod(th, 180)[1] ### th lies between 0 and 180
                return (s1, s2, th)
            else:
                return img.beam2pix(x)

        def pix2gaus(x, location=None, use_wcs=True):
            """ Converts Gaussian parameters in pixels to deg.

            x - (maj [pix], min [pix], pa [deg])
            location - specifies the location in pixels (x, y) for which
                transform is desired
            Input beam angle should be degrees CCW from the +y axis of the image.
            The output beam angle is degrees CCW from North.
            """
            if use_wcs:
                s1, s2, th = x
                brot = self.get_rot(img, location) # rotation delta CCW (in degrees) between N and +y axis of image

                th_rad = th / 180.0 * N.pi
                bmaj = self.pixdist2angdist(img, s1, th, location=location)
                bmin = self.pixdist2angdist(img, s2, th + 90.0, location=location)
                bpa = th - brot
                if bmaj < bmin:
                    bmaj, bmin = bmin, bmaj
                    bpa += 90.0
                bpa = divmod(bpa, 180)[1] ### bpa lies between 0 and 180
                return (bmaj, bmin, bpa)
            else:
                return img.pix2beam(x)

        def pix2coord(pix, location=None, use_wcs=True):
            """Converts size along x and y (in pixels) to size in RA and Dec (in degrees)

            Currently, this function is only used to convert errors on x, y position
            to errors in RA and Dec.
            """
            if use_wcs:
                # Account for projection effects
                x, y = pix
                brot = self.get_rot(img, location) # rotation delta CCW (in degrees) between N and +y axis of image
                ra_dist_pix = N.sqrt( (x * N.cos(brot * N.pi / 180.0))**2 + (y * N.sin(brot * N.pi / 180.0))**2 )
                dec_dist_pix = N.sqrt( (x * N.sin(brot * N.pi / 180.0))**2 + (y * N.cos(brot * N.pi / 180.0))**2 )
                s1 = self.pixdist2angdist(img, ra_dist_pix, 90.0 - brot, location=location)
                s2 = self.pixdist2angdist(img, dec_dist_pix, 0.0 - brot, location=location)
            else:
                x, y = pix
                s1 = abs(x * cdelt1)
                s2 = abs(y * cdelt2)
            return (s1, s2)

        instancemethod = type(t.wcs_pix2sky)
        t.p2s = instancemethod(p2s, t, WCS)
        instancemethod = type(t.wcs_sky2pix)
        t.s2p = instancemethod(s2p, t, WCS)

        img.wcs_obj = t
        img.wcs_obj.acdelt = acdelt
        img.pix2sky = t.p2s
        img.sky2pix = t.s2p
        img.gaus2pix = gaus2pix
        img.pix2gaus = pix2gaus
        img.pix2coord = pix2coord


    def init_beam(self, img):
        """Initialize beam parameters, and conversion routines
        to convert beam to/from pixel coordinates"""
        from const import fwsig
        mylog = mylogger.logging.getLogger("PyBDSM.InitBeam")

        hdr = img.header
        cdelt1, cdelt2 = img.wcs_obj.acdelt[0:2]

        ### define beam conversion routines:
        def beam2pix(x):
            """ Converts beam in deg to pixels. Use when no dependence on
            position is appropriate.

            Input beam angle should be degrees CCW from North at image center.
            The output beam angle is degrees CCW from the +y axis of the image.
            """
            bmaj, bmin, bpa = x
            s1 = abs(bmaj / cdelt1)
            s2 = abs(bmin / cdelt2)
            th = bpa
            return (s1, s2, th)

        def pix2beam(x):
            """ Converts beam in pixels to deg. Use when no dependence on
            position is appropriate.

            Input beam angle should be degrees CCW from the +y axis of the image.
            The output beam angle is degrees CCW from North at image center.
            """
            s1, s2, th = x
            bmaj = abs(s1 * cdelt1)
            bmin = abs(s2 * cdelt2)
            bpa = th
            if bmaj < bmin:
                bmaj, bmin = bmin, bmaj
                bpa += 90.0
            bpa = divmod(bpa, 180)[1] ### bpa lies between 0 and 180
            return [bmaj, bmin, bpa]

        def pixel_beam():
            """Returns the beam in sigma units in pixels"""
            pbeam = beam2pix(img.beam)
            return (pbeam[0] / fwsig, pbeam[1] / fwsig, pbeam[2])

        def pixel_beamarea():
            """Returns the beam area in pixels"""
            pbeam = beam2pix(img.beam)
            return 1.1331 * pbeam[0] * pbeam[1]

        ### Get the beam information from the header
        found = False
        if img.opts.beam is not None:
            beam = img.opts.beam
        else:
            try:
                beam = (hdr['BMAJ'], hdr['BMIN'], hdr['BPA'])
                found = True
            except:
                ### try see if AIPS as put the beam in HISTORY as usual
               for h in hdr.get_history():
                  # Check if h is a string or a FITS Card object (long headers are
                  # split into Cards as of PyFITS 3.0.4)
                  if not isinstance(h, str):
                    hstr = h.value
                  else:
                    hstr = h
                  if N.all(['BMAJ' in hstr, 'BMIN' in hstr, 'BPA' in hstr, 'CLEAN' in hstr]):
                    try:
                        dum, dum, dum, bmaj, dum, bmin, dum, bpa = hstr.split()
                    except ValueError:
                        try:
                            dum, dum, bmaj, dum, bmin, dum, bpa, dum, dum = hstr.split()
                        except ValueError:
                            break
                    beam = (float(bmaj), float(bmin), float(bpa))
                    found = True
            if not found: raise RuntimeError("No beam information found in image header.")

        ### convert beam into pixels (at image center)
        pbeam = beam2pix(beam)
        pbeam = (pbeam[0] / fwsig, pbeam[1] / fwsig, pbeam[2])  # IN SIGMA UNITS

        ### and store it
        img.pix2beam = pix2beam
        img.beam2pix = beam2pix
        img.beam = beam   # FWHM size in degrees
        img.pixel_beam = pixel_beam   # IN SIGMA UNITS in pixels
        img.pixel_beamarea = pixel_beamarea
        mylogger.userinfo(mylog, 'Beam shape (major, minor, pos angle)',
                          '(%s, %s, %s) degrees' % (round(beam[0], 5),
                                                    round(beam[1], 5),
                                                    round(beam[2], 1)))

    def init_freq(self, img):
        """Initialize frequency parameters and store them.

        Basically, PyBDSM uses two frequency parameters:

            img.frequency - the reference frequency in Hz of the ch0 image
            img.freq_pars - the crval, crpix, and cdelt values for the
                            frequency axis in Hz

        If the input frequency info (in the WCS) is not in Hz, it is
        converted.
        """
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            from pywcs import WCS

        mylog = mylogger.logging.getLogger("PyBDSM.InitFreq")
        if img.opts.frequency_sp != None and img.image.shape[1] > 1:
            # If user specifies multiple frequencies, then let
            # collapse.py do the initialization
            img.frequency = img.opts.frequency_sp[0]
            img.freq_pars = (0.0, 0.0, 0.0)
            mylog.info('Using user-specified frequencies.')
        elif img.opts.frequency != None and img.image.shape[1] == 1:
            img.frequency = img.opts.frequency
            img.freq_pars = (img.frequency, 0.0, 0.0)
            mylog.info('Using user-specified frequency.')
        else:
            spec_indx = img.wcs_obj.wcs.spec
            if spec_indx == -1:
                raise RuntimeError('No frequency information found in image header.')
            else:
                # Here we define p2f and f2p to allow pixel to frequency
                # transformations. Transformations for other axes (e.g.,
                # celestial) are striped out.
                #
                # First, convert frequency to Hz if needed:
                img.wcs_obj.wcs.sptr('FREQ-???')
                naxis = img.wcs_obj.wcs.naxis
                def p2f(self, spec_pix):
                    spec_list = []
                    for i in range(naxis):
                        spec_list.append(0)
                    spec_list[spec_indx] = spec_pix
                    spec_pix_arr = N.array([spec_list])
                    freq = self.wcs_pix2sky(spec_pix_arr, 0)
                    return freq.tolist()[0][spec_indx]
                def f2p(self, freq):
                    freq_list = []
                    for i in range(naxis):
                        freq_list.append(0)
                    freq_list[spec_indx] = freq
                    freq_arr = N.array([freq_list])
                    pix = self.wcs_sky2pix(freq_arr, 0)
                    return pix.tolist()[0][spec_indx]
                instancemethod = type(img.wcs_obj.wcs_pix2sky)
                img.wcs_obj.p2f = instancemethod(p2f, img.wcs_obj, WCS)
                instancemethod = type(img.wcs_obj.wcs_sky2pix)
                img.wcs_obj.f2p = instancemethod(f2p, img.wcs_obj, WCS)

                if img.opts.frequency != None:
                    img.frequency = img.opts.frequency
                else:
                    img.frequency = img.wcs_obj.p2f(0)

    def get_equinox(self, img):
        """Gets the equinox from the header.

        Returns float year with code, where code is:
        1 - EQUINOX, EPOCH or RADECSYS keyword not found in header
        0 - EQUINOX found as a numeric value
        1 - EPOCH keyword used for equinox (not recommended)
        2 - EQUINOX found as  'B1950'
        3 - EQUINOX found as  'J2000'
        4 - EQUINOX derived from value of RADECSYS keyword
            'ICRS', 'FK5' ==> 2000,  'FK4' ==> 1950
        """
        code = -1
        year = None
        hdr = img.header
        if 'EQUINOX' in hdr:
            year = hdr['EQUINOX']
            if isinstance(year, str):     # Check for 'J2000' or 'B1950' values
                tst = year[:1]
                if (tst == 'J') or (tst == 'B'):
                    year = float(year[1:])
                    if tst == 'J': code = 3
                    if tst == 'B': code = 2
            else:
                code = 0
        else:
            if 'EPOCH' in hdr: # Check EPOCH if EQUINOX not found
                year = float(hdr['EPOCH'])
                code = 1
            else:
                if 'RADECSYS' in hdr:
                    sys = hdr['RADECSYS']
                    code = 4
                    if sys[:3] == 'ICR': year = 2000.0
                    if sys[:3] == 'FK5': year = 2000.0
                    if sys[:3] == 'FK4': year = 1950.0
        return year, code

    def get_rot(self, img, location=None):
        """Returns CCW rotation angle (in degrees) between N and +y axis of image

        location specifies the location in pixels (x, y) for which angle is desired
        """
        if location == None:
            x1 = img.image.shape[2] / 2.0
            y1 = img.image.shape[3] / 2.0
        else:
            x1, y1 = location
        ra, dec = img.pix2sky([x1, y1])
        delta_dec = self.pixdist2angdist(img, 1.0, 0.0, location=[x1, y1])  # approx. size in degrees of 1 pixel
        if dec + delta_dec > 90.0:
            # shift towards south instead
            delta_dec *= -1.0
        x2, y2 = img.sky2pix([ra, dec + delta_dec])
        try:
            rot_ang_rad = N.arctan2(y2-y1, x2-x1) - N.pi / 2.0
            if delta_dec < 0.0:
                rot_ang_rad -= N.pi
        except:
            rot_ang_rad = 0.0
        return rot_ang_rad * 180.0 / N.pi

    def angdist2pixdist(self, img, angdist, pa, location=None):
        """Returns the distance in pixels for a given angular distance in degrees

        pa - position angle in degrees east of north
        location - x and y location of center
        """
        import functions as func

        if location == None:
            x1 = int(img.image.shape[2] / 2.0)
            y1 = int(img.image.shape[3] / 2.0)
        else:
            x1, y1 = location

        pa_pix = self.get_rot(img, location)
        x0 = x1 - 10.0 * N.sin( (pa + pa_pix) * N.pi / 180.0 )
        y0 = y1 - 10.0 * N.cos( (pa + pa_pix) * N.pi / 180.0 )
        ra0, dec0 = img.pix2sky([x0, y0])
        x2 = x1 + 10.0 * N.sin( (pa + pa_pix) * N.pi / 180.0 )
        y2 = y1 + 10.0 * N.cos( (pa + pa_pix) * N.pi / 180.0 )
        ra2, dec2 = img.pix2sky([x2, y2])

        angdist12 = func.angsep(ra0, dec0, ra2, dec2) # degrees
        pixdist12 = N.sqrt( (x0 - x2)**2 + (y0 - y2)**2 ) # pixels
        if angdist12 > 0.0:
            result = angdist * pixdist12 / angdist12
            if N.isnan(result) or result <= 0.0:
                result = N.mean(img.wcs_obj.acdelt[0:2])
        else:
            result = N.mean(img.wcs_obj.acdelt[0:2])
        return result

    def pixdist2angdist(self, img, pixdist, pa, location=None):
        """Returns the angular distance in degrees for a given distance in pixels

        pa - position angle in degrees CCW from +y axis
        location - x and y location of center
        """
        import functions as func

        if location == None:
            x1 = int(img.image.shape[2] / 2.0)
            y1 = int(img.image.shape[3] / 2.0)
        else:
            x1, y1 = location

        x0 = x1 - pixdist / 2.0 * N.sin(pa * N.pi / 180.0)
        y0 = y1 - pixdist / 2.0 * N.cos(pa * N.pi / 180.0)
        ra0, dec0 = img.pix2sky([x0, y0])
        x2 = x1 + pixdist / 2.0 * N.sin(pa * N.pi / 180.0)
        y2 = y1 + pixdist / 2.0 * N.cos(pa * N.pi / 180.0)
        ra2, dec2 = img.pix2sky([x2, y2])

        angdist12 = func.angsep(ra0, dec0, ra2, dec2) # degrees
        return angdist12



