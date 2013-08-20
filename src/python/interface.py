"""Interface module.

The interface module handles all functions typically needed by the user in an
interactive environment such as IPython. Many are also used by the
custom IPython shell defined in pybdsm.py.

"""

def process(img, **kwargs):
    """Find and measure sources in an image.

    This function is used by process_image in __init__.py and by process_image
    in pybdsm.py. It is also used as a method of the Image object in image.py
    to allow reprocessing of existing Image objects with the command
    img.process().

    Any options given as keyword arguments will override existing ones stored
    in img.opts.
    """
    from . import default_chain, _run_op_list
    from image import Image
    import mylogger

    # First, reset img to initial state (in case img is being reprocessed)
    if hasattr(img, 'use_io'): del img.use_io
    if hasattr(img, 'sources'): del img.sources
    if hasattr(img, 'dsources'): del img.dsources
    if hasattr(img, 'gaussians'): del img.gaussians
    if hasattr(img, 'atrous_gaussians'): del img.atrous_gaussians
    if hasattr(img, 'islands'): del img.islands
    if hasattr(img, 'ch0'): del img.ch0
    if hasattr(img, 'image'): del img.image
    if hasattr(img, 'rms'): del img.rms
    if hasattr(img, 'mean'): del img.mean
    if hasattr(img, 'rms_QUV'): del img.rms_QUV
    if hasattr(img, 'mean_QUV'): del img.mean_QUV
    if hasattr(img, 'resid_gaus'): del img.resid_gaus
    if hasattr(img, 'model_gaus'): del img.model_gaus
    if hasattr(img, 'resid_shap'): del img.resid_shap
    if hasattr(img, 'model_shap'): del img.model_shap
    if hasattr(img, 'mask'): del img.mask
    if hasattr(img, 'rms_mask'): del img.rms_mask
    if hasattr(img, 'completed_Ops'): del img.completed_Ops

    try:
        # set options if given
        if len(kwargs) > 0:
            set_pars(img, **kwargs)
    except RuntimeError, err:
        # Catch and log error
        mylog.error(str(err))

        # Re-throw error if the user is not in the interactive shell
        if img._is_interactive_shell:
            return False
        else:
            raise

    # Start up logger. We need to initialize it each time process() is
    # called, in case the quiet or debug options have changed
    log = img.opts.filename + '.pybdsm.log'
    img.log = ''
    mylogger.init_logger(log, quiet=img.opts.quiet,
                         debug=img.opts.debug)
    add_break_to_logfile(log)
    mylog = mylogger.logging.getLogger("PyBDSM.Process")
    mylog.info("Processing "+img.opts.filename)

    # Run all the op's
    try:
        # Run op's in chain
        op_chain = get_op_chain(img)
        _run_op_list(img, op_chain)
        return True
    except RuntimeError, err:
        # Catch and log error
        mylog.error(str(err))

        # Re-throw error if the user is not in the interactive shell
        if img._is_interactive_shell:
            return False
        else:
            raise
    except KeyboardInterrupt:
        mylogger.userinfo(mylog, "\n\033[31;1mAborted\033[0m")
        return False

def get_op_chain(img):
    """Determines the optimal Op chain for an Image object.

    This is useful when reprocessing an Image object. For example,
    if Gaussians were already fit, but the user now wants to use
    shapelets, we do not need to re-run Op_gausfit, etc. At the
    moment, this just returns the default Op chain from __init__.py.
    """
    from . import default_chain

    return default_chain
#     prev_opts = img._prev_opts
#     new_opts = img.opts.to_dict()
#
#     # Find whether new opts differ from previous opts
#     for k, v in prev_opts.iteritems():
#         if v != new_opts[k]:
#             if k == 'rms_box':

    # If filename, beam, trim_box differ, start from readimage
    # Elif shapelet_do, etc. differ, start from there


def load_pars(filename):
    """Load parameters from a save file or dictionary.

    The file must be a pickled opts dictionary.

    filename - name of options file to load.
    Returns None (and original error) if no file can be loaded successfully.
    """
    from image import Image
    import mylogger
    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    # First, check if input is a dictionary
    if isinstance(filename, dict):
        timg = Image(filename)
        return timg, None
    else:
        try:
            pkl_file = open(filename, 'rb')
            pars = pickle.load(pkl_file)
            pkl_file.close()
            timg = Image(pars)
            print "--> Loaded parameters from file '" + filename + "'."
            return timg, None
        except Exception, err:
            return None, err

def save_pars(img, savefile=None, quiet=False):
    """Save parameters to a file.

    The save file is a "pickled" opts dictionary.
    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    import tc
    import sys

    if savefile == None or savefile == '':
        savefile = img.opts.filename + '.pybdsm.sav'

    # convert opts to dictionary
    pars = img.opts.to_dict()
    output = open(savefile, 'wb')
    pickle.dump(pars, output)
    output.close()
    if not quiet:
        print "--> Saved parameters to file '" + savefile + "'."

def list_pars(img, opts_list=None, banner=None, use_groups=True):
    """Lists all parameters for the Image object.

    opts_list - a list of the parameter names to list;
                if None, all parameters are used.
    banner - banner text to place at top of listing.
    use_groups - whether to use the group information for each
                 parameter.
    """
    import tc
    import sys

    # Get all options as a list sorted by name
    opts = img.opts.to_list()

    # Filter list
    if opts_list != None:
        opts_temp = []
        for o in opts:
            if o[0] in opts_list:
                opts_temp.append(o)
        opts = opts_temp

    # Move filename, infile, outfile to front of list
    for o in opts:
        if o[0] == 'filename' or o[0] == 'infile' or o[0] == 'outfile':
            opts.remove(o)
            opts.insert(0, o)

    # Now group options with the same "group" together.
    if use_groups:
        opts = group_opts(opts)

    # Finally, print options, values, and doc strings to screen
    print_opts(opts, img, banner=banner)


def set_pars(img, **kwargs):
    """Set parameters using arguments instead of using a dictionary.

    Allows partial names for parameters as long as they are unique. Parameters
    are set to default values if par = ''.
    """
    import re
    import sys
    from image import Image

    # Enumerate all options
    opts = img.opts.get_names()

    # Check that parameters are valid options and are unique
    full_key = []
    for i, key in enumerate(kwargs):
        chk_key = checkpars(opts, key)
        if chk_key == []:
            raise RuntimeError("Input parameter '" + key + "' not recognized.")
        if len(chk_key) > 1 and key not in opts:
            raise RuntimeError("Input parameter '" + key + "' matches to more than one "\
                         "possible parameter:\n " + "\n ".join(chk_key))
        if key in opts:
            full_key.append(key)
        else:
            full_key.append(chk_key[0])

    # Build options dictionary
    pars = {}
    for i, key in enumerate(kwargs):
        if kwargs[key] == '':
            temp_img = Image({'filename':''})
            opt_names = temp_img.opts.get_names
            for k in opt_names:
                if key == k:
                    kwargs[key] = temp_img.opts.__getattribute__(k)
        pars.update({full_key[i]: kwargs[key]})

    # Finally, set the options
    img.opts.set_opts(pars)


def group_opts(opts):
    """Sorts options by group (as defined in opts.py).

    Returns a list of options, with suboptions arranged in a list inside the
    main list and directly following the main options. Options belonging to the
    "hidden" group are excluded from the returned list (as defined in opts.py).
    """
    groups = []
    gp = []
    for i in range(len(opts)):
        grp = opts[i][1].group()
        if grp != None and grp not in groups:
            groups.append(opts[i][1].group())

    groups.sort()

    # Now, make a list for each group with its options. Don't include
    # "hidden" options, as they should never by seen by the user.
    for g in groups:
        g_list = []
        for i in range(len(opts)):
            if isinstance(opts[i], tuple):
                if g == str(opts[i][1].group()):
                    g_list.append(opts[i])
        for gs in g_list:
            opts.remove(gs)

        for i in range(len(opts)):
            if g == str(opts[i][0]) and g != 'hidden':
                opts.insert(i+1, g_list)
                break
    return opts


def print_opts(grouped_opts_list, img, banner=None):
    """Print options to screen.

    Options can be sorted by group (defined in opts.py) previously defined by
    group_opts. Output of grouped items is suppressed if parent option is
    False. The layout is as follows:

      [20 spaces par name with ...] = [at least 49 spaces for value]
                                      [at least 49 spaces for doc]

    When more than one line is required for the doc, the next line is:

      [25 blank spaces][at least 47 spaces for doc]

    As in casapy, print non-defaults in blue, options with suboptions in
    47m and suboptions in green. Option Values are printed in bold, to help
    to distinguish them from the descriptions. NOTE: in iTerm, one needs
    to set the bold color in the profiles to white, as it defaults to red,
    which is a bit hard on the eyes in this case.
    """
    from image import Image
    import os
    import functions as func

    termy, termx = func.getTerminalSize() # note: returns row, col -> y, x
    minwidth = 28 # minimum width for parameter names and values

    # Define colors for output
    dc = '\033[1;34m' # Blue: non-default option text color
    ec = '\033[0;47m' # expandable option text color
    sc = '\033[0;32m' # Green: suboption text color
    nc = '\033[0m'    # normal text color
    ncb = '\033[1m'    # normal text color bold

    if banner != None:
        print banner
    spcstr = ' ' * minwidth # spaces string for second or later lines
    infix = nc + ': ' + nc # infix character used to separate values from comments
    print '=' * termx # division string for top of parameter listing
    for indx, o in enumerate(grouped_opts_list):
        if isinstance(o, tuple):
            # Print main options, which are always tuples, before printing
            # suboptions (if any).
            k = o[0]
            v = o[1]
            val = img.opts.__getattribute__(k)
            v1 = v2 = ''
            if val == v._default:
                # value is default
                v1 = ncb
                v2 = nc
            else:
                # value is non-default
                v1 = dc
                v2 = nc
            if isinstance(val, str):
                valstr = v1 + repr(val) + v2
                if k == 'filename':
                    # Since we can check whether filename is valid,
                    # do so here and print in red if not.
                    if not os.path.exists(val):
                        valstr = '\033[31;1m' + repr(val) + nc
                width_par_val = max(minwidth, len(k) + len(str(val)) + 5)
            else:
                if isinstance(val, float):
                    val = round_float(val)
                if isinstance(val, tuple):
                    val = round_tuple(val)
                valstr = v1 + str(val) + v2
                width_par_val = max(minwidth, len(k) + len(str(val)) + 4)
            width_desc = max(termx - width_par_val - 3, 44)
            # Get the option description text from the doc string, which
            # is defined in opts.py. By convention, print_opts will only
            # show the short description; help('option_name') will
            # print both the short and long description. The versions
            # are separated in the doc string by '\n', which is split
            # on here:
            desc_text = wrap(str(v.doc()).split('\n')[0], width_desc)
            fmt = '%' + str(minwidth) + 's' + infix + '%44s'

            # Now loop over lines of description
            if indx < len(grouped_opts_list)-1:
                # Here we check if next entry in options list is a tuple or a
                # list.  If it is a list, then the current option has
                # suboptions and should be in the ec color. Since we check the
                # next option, we can't do this if we let indx go to the end.
                if isinstance(grouped_opts_list[indx+1], tuple):
                    parvalstr = nc + k + nc + ' ..'
                else:
                    parvalstr = ec + k + nc + ' ..'
            else:
                # Since this is the last entry in the options list and is a
                # tuple, it cannot be an expandable option, so make it nc color
                 parvalstr = nc + k + nc + ' ..'
            if "'" in valstr:
                len_without_formatting = len(k) + len(str(val)) + 5
            else:
                len_without_formatting = len(k) + len(str(val)) + 4
            for i in range(len_without_formatting, minwidth):
                parvalstr += '.'
            parvalstr += ' ' + valstr
            if "'" not in valstr:
                parvalstr += ' '
            for dt_indx, dt in enumerate(desc_text):
                if dt_indx == 0:
                    print fmt % (parvalstr.ljust(minwidth), dt.ljust(44))
                else:
                    print nc + spcstr + '   %44s' % dt.ljust(44)
        else:
            # Print suboptions, indented 2 spaces from main options in sc color
            parent_opt = grouped_opts_list[indx-1]
            parent_val = img.opts.__getattribute__(parent_opt[0])
            if parent_val == True:
                for og in o:
                    k = og[0]
                    v = og[1]
                    val = img.opts.__getattribute__(k)
                    v1 = v2 = ''
                    if val == v._default:
                        # value is default
                        v1 = ncb
                        v2 = nc
                    else:
                        # value is non-default
                        v1 = dc
                        v2 = nc
                    if isinstance(val, str):
                        valstr = v1 + repr(val) + v2
                        width_par_val = max(minwidth, len(k) + len(str(val)) + 7)
                    else:
                        if isinstance(val, float):
                            val = round_float(val)
                        if k == 'beam_spectrum' and val != None:
                            val = round_list_of_tuples(val)
                        if k == 'frequency_sp' and val != None:
                            val = round_list(val)
                        valstr = v1 + str(val) + v2
                        width_par_val = max(minwidth, len(k) + len(str(val)) + 6)
                    width_desc = max(termx - width_par_val - 3, 44)
                    desc_text = wrap(str(v.doc()).split('\n')[0], width_desc)
                    fmt = '  ' + '%' + str(minwidth) + 's' + infix + '%44s'
                    parvalstr = sc + k + nc + ' ..'
                    if "'" in valstr:
                        len_without_formatting = len(k) + len(str(val)) + 7
                    else:
                        len_without_formatting = len(k) + len(str(val)) + 6
                    for i in range(len_without_formatting, minwidth):
                        parvalstr += '.'
                    parvalstr += ' ' + valstr
                    if "'" not in valstr:
                        parvalstr += ' '
                    for dt_indx, dt in enumerate(desc_text):
                        if dt_indx == 0:
                            print fmt % (parvalstr.ljust(minwidth-2), dt.ljust(44))
                        else:
                            print nc + spcstr + '   %44s' % dt.ljust(44)


def wrap(text, width=80):
    """Wraps text to given width and returns list of lines."""
    lines = []
    for paragraph in text.split('\n'):
        line = []
        len_line = 0
        for word in paragraph.split(' '):
            word.strip()
            len_word = len(word)
            if len_line + len_word <= width:
                line.append(word)
                len_line += len_word + 1
            else:
                lines.append(' '.join(line))
                line = [word]
                len_line = len_word + 1
        lines.append(' '.join(line))
    return lines


def checkpars(lines, regex):
    """Checks that parameters are unique"""
    import re
    result = []
    for l in lines:
        match = re.match(regex,l)
        if match:
            result += [l]
    return result


def in_ipython():
    """Checks if interpreter is IPython."""
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True


def raw_input_no_history(prompt):
    """Removes user input from readline history."""
    import readline
    input = raw_input(prompt)
    if input != '':
        readline.remove_history_item(readline.get_current_history_length()-1)
    return input


# The following functions just make the printing of
# parameters look better
def round_tuple(val):
    valstr_list = []
    for v in val:
        vstr = '%s' % (round(v, 5))
        if len(vstr) > 7:
            vstr = '%.5f' % (v,)
        valstr_list.append(vstr)
    valstr = '(' + ','.join(valstr_list) + ')'
    return valstr

def round_float(val):
    vstr = '%s' % (round(val, 5))
    if len(vstr) > 7 and val < 1e3:
        vstr = '%.5f' % (val,)
    elif len(vstr) > 7 and val >= 1e3:
        vstr = '%.2e' % (val,)
    return vstr

def round_list(val):
    valstr_list = []
    for v in val:
        valstr_list.append('%.2e' % (v,))
    valstr = '[' + ','.join(valstr_list) + ']'
    return valstr

def round_list_of_tuples(val):
    valstr_list = []
    valstr_list_tot = []
    for l in val:
        for v in l:
            vstr = '%s' % (round(v, 5))
            if len(vstr) > 7:
                vstr = '%.5f' % (v,)
            valstr_list.append(vstr)
        valstr = '(' + ','.join(valstr_list) + ')'
        valstr_list_tot.append(valstr)
    valstr = '[' + ','.join(valstr_list_tot) + ']'
    return valstr

# The following functions give convenient access to the output functions in
# output.py
def export_image(img, outfile=None, img_format='fits',
                 img_type='gaus_resid', clobber=False):
    """Write an image to a file. Returns True if successful, False if not.

    outfile - name of resulting file; if None, file is
    named automatically.
    img_type - type of image to export; see below
    img_format - format of resulting file: 'fits' or 'casa'
    incl_wavelet - include wavelet Gaussians in model
                     and residual images?
    clobber - overwrite existing file?

    The following images may be exported:
        'ch0' - image used for source detection
        'rms' - rms map image
        'mean' - mean map image
        'pi' - polarized intensity image
        'gaus_resid' - Gaussian model residual image
        'gaus_model' - Gaussian model image
        'shap_resid' - Shapelet model residual image
        'shap_model' - Shapelet model image
        'psf_major' - PSF major axis FWHM image (FWHM in arcsec)
        'psf_minor' - PSF minor axis FWHM image (FWHM in arcsec)
        'psf_pa' - PSF position angle image (degrees east of north)
        'psf_ratio' - PSF peak-to-total flux ratio (in units of 1/beam)
        'psf_ratio_aper' - PSF peak-to-aperture flux ratio (in units of 1/beam)
    """
    import os
    import functions as func
    from const import fwsig
    import mylogger

    mylog = mylogger.logging.getLogger("PyBDSM."+img.log+"ExportImage")

    # First some checking:
    if not 'gausfit' in img.completed_Ops and 'gaus' in img_type:
        print '\033[91mERROR\033[0m: Gaussians have not been fit. Please run process_image first.'
        return False
    elif not 'shapelets' in img.completed_Ops and 'shap' in img_type:
        print '\033[91mERROR\033[0m: Shapelets have not been fit. Please run process_image first.'
        return False
    elif not 'polarisation' in img.completed_Ops and 'pi' in img_type:
        print '\033[91mERROR\033[0m: Polarization properties have not been calculated. Please run process_image first.'
        return False
    elif not 'psf_vary' in img.completed_Ops and 'psf' in img_type:
        print '\033[91mERROR\033[0m: PSF variations have not been calculated. Please run process_image first.'
        return False
    elif not 'collapse' in img.completed_Ops and 'ch0' in img_type:
        print '\033[91mERROR\033[0m: ch0 image has not been calculated. Please run process_image first.'
        return False
    elif not 'rmsimage' in img.completed_Ops and ('rms' in img_type or 'mean' in img_type):
        print '\033[91mERROR\033[0m: Mean and rms maps have not been calculated. Please run process_image first.'
        return False
    elif not 'make_residimage' in img.completed_Ops and ('resid' in img_type or 'model' in img_type):
        print '\033[91mERROR\033[0m: Residual and model maps have not been calculated. Please run process_image first.'
        return False
    format = img_format.lower()
    if (format in ['fits', 'casa']) == False:
        print '\033[91mERROR\033[0m: img_format must be "fits" or "casa"'
        return False
    if format == 'casa':
        print "\033[91mERROR\033[0m: Only img_format = 'fits' is supported at the moment"
        return False
    filename = outfile
    if filename == None or filename == '':
        filename = img.imagename + '_' + img_type + '.' + format
    if os.path.exists(filename) and clobber == False:
        print '\033[91mERROR\033[0m: File exists and clobber = False.'
        return False
    if format == 'fits':
        use_io = 'fits'
    if format == 'casa':
        use_io = 'rap'
    bdir = ''
    try:
        if img_type == 'ch0':
            func.write_image_to_file(use_io, filename,
                                     img.ch0, img, bdir,
                                     clobber=clobber)
        elif img_type == 'rms':
            func.write_image_to_file(use_io, filename,
                                     img.rms, img, bdir,
                                     clobber=clobber)
        elif img_type == 'mean':
            func.write_image_to_file(use_io, filename,
                                     img.mean, img, bdir,
                                     clobber=clobber)
        elif img_type == 'pi':
            func.write_image_to_file(use_io, filename,
                                     img.ch0_pi, img, bdir,
                                     clobber=clobber)
        elif img_type == 'psf_major':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_maj*fwsig, img, bdir,
                                     clobber=clobber)
        elif img_type == 'psf_minor':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_min*fwsig, img, bdir,
                                     clobber=clobber)
        elif img_type == 'psf_pa':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_pa, img, bdir,
                                     clobber=clobber)
        elif img_type == 'psf_ratio':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_ratio, img, bdir,
                                     clobber=clobber)
        elif img_type == 'psf_ratio_aper':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_ratio_aper, img, bdir,
                                     clobber=clobber)
        elif img_type == 'gaus_resid':
            im = img.resid_gaus
            func.write_image_to_file(use_io, filename,
                                     im, img, bdir,
                                     clobber=clobber)
        elif img_type == 'gaus_model':
            im = img.model_gaus
            func.write_image_to_file(use_io, filename,
                                     im, img, bdir,
                                     clobber=clobber)
        elif img_type == 'shap_resid':
            func.write_image_to_file(use_io, filename,
                                     img.resid_shap, img, bdir,
                                     clobber=clobber)
        elif img_type == 'shap_model':
            func.write_image_to_file(use_io, filename,
                                     img.model_shap, img, bdir,
                                     clobber=clobber)
        else:
            print "\n\033[91mERROR\033[0m: img_type not recognized."
            return False
        if filename == 'SAMP':
            print '--> Image sent to SMAP hub'
        else:
            print '--> Wrote file ' + repr(filename)
        return True
    except RuntimeError, err:
        # Catch and log error
        mylog.error(str(err))

        # Re-throw error if the user is not in the interactive shell
        if img._is_interactive_shell:
            return False
        else:
            raise
    except KeyboardInterrupt:
        mylogger.userinfo(mylog, "\n\033[31;1mAborted\033[0m")
        return False


def write_catalog(img, outfile=None, format='bbs', srcroot=None, catalog_type='gaul',
               bbs_patches=None, incl_chan=False, incl_empty=False, clobber=False,
               force_output=False, correct_proj=True):
    """Write the Gaussian, source, or shapelet list to a file. Returns True if
    successful, False if not.

    filename - name of resulting file; if None, file is
               named automatically.
    catalog_type - type of catalog
        "gaul"  - Gaussian list
        "srl"   - Source list
        "shap"  - Shapelet list ("fits" format only)
    format - format of output list. Supported formats are:
        "fits"  - FITS binary table
        "ascii" - ASCII text file
        "bbs"   - BBS sky model (Gaussian list only)
        "ds9"   - ds9 region file
        "star"  - AIPS STAR file (Gaussian list only)
        "kvis"  - kvis file (Gaussian list only)
        "sagecal" - Sagecal file (Gaussian list only)
    srcroot - root for source and patch names (BBS/ds9 only);
              if None, the srcroot is chosen automatically
    bbs_patches - type of patches to use:
        None - no patches
        "gaussian" - each Gaussian gets its own patch
        "single"   - all Gaussians are put into a single
                     patch
        "source"   - sources are grouped by source into patches
    incl_chan - Include fluxes for each channel?
    incl_empty - Include islands without any valid Gaussians (source list only)?
    sort_by - Property to sort output list by:
        "flux" - sort by total integrated flux, largest first
        "indx" - sort by Gaussian and island or source index, smallest first
    force_output - Force the creation of a catalog, even if it is empty
    correct_proj - Correct source parameters for image projection effects (BBS only)?
    clobber - Overwrite existing file?
    """
    import output

    # First some checking:
    if not 'gausfit' in img.completed_Ops:
        print '\033[91mERROR\033[0m: Image has not been fit. Please run process_image first.'
        return False
    if catalog_type == 'shap' and not 'shapelets' in img.completed_Ops:
            print '\033[91mERROR\033[0m: Image has not been decomposed into shapelets. Please run process_image first.'
            return False
    if catalog_type == 'srl' and not 'gaul2srl' in img.completed_Ops:
            print '\033[91mERROR\033[0m: Gaussians have not been grouped into sources. Please run process_image first.'
            return False
    format = format.lower()
    patch = bbs_patches
    filename = outfile
    if isinstance(patch, str):
        patch = patch.lower()
    if (format in ['fits', 'ascii', 'bbs', 'ds9', 'star',
                   'kvis', 'sagecal']) == False:
        print '\033[91mERROR\033[0m: format must be "fits", '\
            '"ascii", "ds9", "star", "kvis",  or "bbs"'
        return False
    if (patch in [None, 'gaussian', 'single', 'source']) == False:
        print '\033[91mERROR\033[0m: patch must be None, '\
            '"gaussian", "source", or "single"'
        return False
    if (catalog_type in ['gaul', 'srl', 'shap']) == False:
        print '\033[91mERROR\033[0m: catalog_type must be "gaul", '\
              '"srl", or "shap"'
        return False
    if (len(img.sources) == 0 and not incl_empty) or (len(img.sources) == 0 and len(img.dsources) == 0 and incl_empty):
        if not force_output:
            print 'No sources were found in the image. Output file not written.'
            return False
    if filename == '': filename = None

    # Now go format by format and call appropriate function
    if filename == 'SAMP':
        import tempfile
        import functions as func
        import os
        if not hasattr(img,'samp_client'):
            s, private_key = func.start_samp_proxy()
            img.samp_client = s
            img.samp_key = private_key

        # Broadcast fits table to SAMP Hub
        tfile = tempfile.NamedTemporaryFile(delete=False)
        filename = output.write_fits_list(img, filename=tfile.name,
                                             incl_chan=incl_chan, incl_empty=incl_empty,
                                             clobber=True, objtype=catalog_type)
        table_name = 'PyBDSM '+ catalog_type + ' table'
        if catalog_type == 'srl':
            img.samp_srl_table_url = 'file://' + os.path.abspath(tfile.name)
        if catalog_type == 'gaul':
            img.samp_gaul_table_url = 'file://' + os.path.abspath(tfile.name)
        func.send_fits_table(img.samp_client, img.samp_key, table_name, tfile.name)
        print '--> Table sent to SMAP hub'
        return True
    if format == 'fits':
        filename = output.write_fits_list(img, filename=filename,
                                             incl_chan=incl_chan, incl_empty=incl_empty,
                                             clobber=clobber, objtype=catalog_type)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber = False.'
            return False
        else:
            print '--> Wrote FITS file ' + repr(filename)
            return True
    if format == 'ascii':
        filename = output.write_ascii_list(img, filename=filename,
                                              incl_chan=incl_chan, incl_empty=incl_empty,
                                              sort_by='index',
                                              clobber=clobber, objtype=catalog_type)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber = False.'
            return False
        else:
            print '--> Wrote ASCII file ' + repr(filename)
            return True
    if format == 'bbs':
        if catalog_type != 'gaul':
            print "\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with BBS files."
            return False
        filename = output.write_bbs_gaul(img, filename=filename,
                                            srcroot=srcroot, incl_empty=incl_empty,
                                            patch=patch, correct_proj=correct_proj,
                                            sort_by='flux',
                                            clobber=clobber)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber = False.'
            return False
        else:
            print '--> Wrote BBS sky model ' + repr(filename)
            return True
    if format == 'sagecal':
        if catalog_type != 'gaul':
            print "\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with Sagecal files."
            return False
        filename = output.write_lsm_gaul(img, filename=filename,
                                            srcroot=srcroot, incl_empty=incl_empty,
                                            patch=patch,
                                            sort_by='flux',
                                            clobber=clobber)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber = False.'
            return False
        else:
            print '--> Wrote Sagecal lsm file ' + repr(filename)
            return True
    if format == 'ds9':
        filename = output.write_ds9_list(img, filename=filename,
                                            srcroot=srcroot, incl_empty=incl_empty,
                                            clobber=clobber, objtype=catalog_type)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber = False.'
            return False
        else:
            print '--> Wrote ds9 region file ' + repr(filename)
            return True
    if format == 'star':
        if catalog_type != 'gaul':
            print "\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with star files."
            return False
        filename = output.write_star(img, filename=filename,
                                        clobber=clobber)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber = False.'
            return False
        else:
            print '--> Wrote AIPS STAR file ' + repr(filename)
            return True
    if format == 'kvis':
        if catalog_type != 'gaul':
            print "\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with kvis files."
            return False
        filename = output.write_kvis_ann(img, filename=filename,
                                            clobber=clobber)
        if filename == None:
            print '\033[91mERROR\033[0m: File exists and clobber=False.'
            return False
        else:
            print '--> Wrote kvis file ' + repr(filename)
            return True
    # if format == 'casabox':
    #     filename = output.write_casa_gaul(img, filename=filename,
    #                                   incl_wavelet=incl_wavelet,
    #                                   clobber=clobber)
    #     if filename == None:
    #         print '\033[91mERROR\033[0m: File exists and clobber=False.'
    #     else:
    #         print '--> Wrote CASA clean box file ' + filename

def add_break_to_logfile(logfile):
    f = open(logfile, 'a')
    f.write('\n' + '='*72 + '\n')
    f.close()
