"""Interface module.

The interface module handles all functions typically needed by the user in an
interactive environment such as IPython. Many are also used by the
custom IPython shell defined in pybdsf.

"""
from __future__ import print_function
from __future__ import absolute_import

try:
    # For Python 2, use raw_input() for input()
    input = raw_input
except NameError:
    pass


def process(img, **kwargs):
    """Find and measure sources in an image.

    This function is used by process_image in __init__.py and by process_image
    in pybdsf. It is also used as a method of the Image object in image.py
    to allow reprocessing of existing Image objects with the command
    img.process().

    Any options given as keyword arguments will override existing ones stored
    in img.opts.
    """
    from . import default_chain, _run_op_list
    from .image import Image
    from . import mylogger
    from .functions import set_up_output_paths
    import os

    # Start up logger. We need to initialize it each time process() is
    # called, in case the quiet or debug options have changed
    _, basedir = set_up_output_paths(img.opts)
    basename = os.path.basename(img.opts.filename) + '.pybdsf.log'
    logfilename = os.path.join(basedir, basename)
    img.log = ''
    mylogger.init_logger(logfilename, quiet=img.opts.quiet,
                             debug=img.opts.debug)
    add_break_to_logfile(logfilename)
    mylog = mylogger.logging.getLogger("PyBDSF.Process")
    mylog.info("Processing "+img.opts.filename)

    try:
        # set options if given
        if len(kwargs) > 0:
            set_pars(img, **kwargs)
    except RuntimeError as err:
        # Catch and log error
        mylog.error(str(err))

        # Re-throw error if the user is not in the interactive shell
        if img._is_interactive_shell:
            return False
        else:
            raise

    # Run all the op's
    try:
        # Run op's in chain
        img, op_chain = get_op_chain(img)
        if op_chain is not None:
            _run_op_list(img, op_chain)
            img._prev_opts = img.opts.to_dict()
        return True
    except RuntimeError as err:
        # Catch and log error
        mylog.error(str(err))

        # Re-throw error if the user is not in the interactive shell
        if img._is_interactive_shell:
            return False
        else:
            raise
    except KeyboardInterrupt:
        mylogger.userinfo(mylog, "\n\033[31;1mAborted\033[0m")
        if img._is_interactive_shell:
            return False
        else:
            raise

def get_op_chain(img):
    """Determines the optimal Op chain for an Image object.

    This is useful when reprocessing an Image object. For example,
    if Gaussians were already fit, but the user now wants to use
    shapelets, we do not need to re-run Op_gausfit, etc.

    Note that any new options added to opts.py should also be
    added here. If not, a full reprocessing will be done if the
    new option is changed.
    """
    from . import default_chain
    Op_chain = default_chain[:]
    Op_names = ['readimage',
                'collapse',
                'preprocess',
                'rmsimage',
                'threshold',
                'islands',
                'gausfit',
                'wavelet_atrous',
                'shapelets',
                'gaul2srl',
                'spectralindex',
                'polarisation',
                'make_residimage',
                'psf_vary',
                'outlist',
                'cleanup']
    prev_opts = img._prev_opts
    if prev_opts is None:
        return img, default_chain
    new_opts = img.opts.to_dict()

    # Set the hidden options, which should include any option whose change
    # should not trigger a process_image action
    hidden_opts = img.opts.get_names(group='hidden')
    hidden_opts.append('advanced_opts')
    hidden_opts.append('flagging_opts')
    hidden_opts.append('multichan_opts')
    hidden_opts.append('output_opts')

    # Define lists of options for each Op. Some of these can be defined
    # using the "group" parameter of each option.
    #
    # Op_readimage()
    readimage_opts = ['filename', 'beam', 'trim_box', 'frequency',
                      'beam_spectrum', 'frequency_sp']

    # Op_collapse()
    collapse_opts = img.opts.get_names(group='multichan_opts')
    collapse_opts.append('polarisation_do')
    collapse_opts += readimage_opts

    # Op_preprocess()
    preprocess_opts = ['kappa_clip', 'polarisation_do']
    preprocess_opts += collapse_opts

    # Op_rmsimage()
    rmsimage_opts = ['rms_box', 'rms_box_bright', 'adaptive_rms_box',
                     'mean_map', 'rms_map', 'adaptive_thresh', 'rms_box_bright']
    rmsimage_opts += preprocess_opts

    # Op_threshold()
    threshold_opts = ['thresh', 'thresh_pix', 'thresh_isl']
    threshold_opts += rmsimage_opts

    # Op_islands()
    islands_opts = threshold_opts
    islands_opts.append('minpix_isl')

    # Op_gausfit()
    gausfit_opts = ['verbose_fitting']
    gausfit_opts += islands_opts
    gausfit_opts += img.opts.get_names(group='flagging_opts')

    # Op_wavelet_atrous()
    wavelet_atrous_opts = img.opts.get_names(group='atrous_do')
    wavelet_atrous_opts.append('atrous_do')
    wavelet_atrous_opts += gausfit_opts

    # Op_shapelets()
    shapelets_opts = img.opts.get_names(group='shapelet_do')
    shapelets_opts.append('shapelet_do')
    shapelets_opts += islands_opts

    # Op_gaul2srl()
    gaul2srl_opts = ['group_tol', 'group_by_isl', 'group_method']
    gaul2srl_opts += gausfit_opts
    gaul2srl_opts += wavelet_atrous_opts

    # Op_spectralindex()
    spectralindex_opts = img.opts.get_names(group='spectralindex_do')
    spectralindex_opts.append('spectralindex_do')
    spectralindex_opts += gaul2srl_opts

    # Op_polarisation()
    polarisation_opts = img.opts.get_names(group='polarisation_do')
    polarisation_opts.append('polarisation_do')
    polarisation_opts += gaul2srl_opts

    # Op_make_residimage()
    make_residimage_opts = ['fittedimage_clip']
    make_residimage_opts += gausfit_opts
    make_residimage_opts += wavelet_atrous_opts
    make_residimage_opts += shapelets_opts

    # Op_psf_vary()
    psf_vary_opts = img.opts.get_names(group='psf_vary_do')
    psf_vary_opts.append('psf_vary_do')
    psf_vary_opts += gaul2srl_opts

    # Op_outlist() and Op_cleanup() are always done.

    # Find whether new opts differ from previous opts (and are not hidden
    # opts, which should not be checked). If so, found = True and we reset
    # the relevant image parameters and add the relevant Op to the Op_chain.
    re_run = False
    found = False
    for k, v in prev_opts.items():
        if v != new_opts[k] and k not in hidden_opts:
            re_run = True
            if k in readimage_opts:
                if hasattr(img, 'use_io'): del img.use_io
                if hasattr(img, 'image_arr'): del img.image_arr
                while 'readimage' in img.completed_Ops:
                    img.completed_Ops.remove('readimage')
                found = True
            if k in collapse_opts:
                if hasattr(img, 'mask_arr'): del img.mask_arr
                if hasattr(img, 'ch0_arr'): del img.ch0_arr
                while 'collapse' in img.completed_Ops:
                    img.completed_Ops.remove('collapse')
                found = True
            if k in preprocess_opts:
                while 'preprocess' in img.completed_Ops:
                    img.completed_Ops.remove('preprocess')
                found = True
            if k in rmsimage_opts:
                if hasattr(img, 'rms_arr'): del img.rms_arr
                if hasattr(img, 'mean_arr'): del img.mean_arr
                if hasattr(img, 'rms_Q_arr'): del img.rms_Q_arr
                if hasattr(img, 'mean_Q_arr'): del img.mean_Q_arr
                if hasattr(img, 'rms_U_arr'): del img.rms_U_arr
                if hasattr(img, 'mean_U_arr'): del img.mean_U_arr
                if hasattr(img, 'rms_V_arr'): del img.rms_V_arr
                if hasattr(img, 'mean_V_arr'): del img.mean_V_arr
                if hasattr(img, '_adapt_rms_isl_pos'): del img._adapt_rms_isl_pos
                while 'rmsimage' in img.completed_Ops:
                    img.completed_Ops.remove('rmsimage')
                found = True
            if k in threshold_opts:
                while 'threshold' in img.completed_Ops:
                    img.completed_Ops.remove('threshold')
                found = True
            if k in islands_opts:
                if hasattr(img, 'islands'): del img.islands
                while 'islands' in img.completed_Ops:
                    img.completed_Ops.remove('islands')
                found = True
            if k in gausfit_opts:
                if hasattr(img, 'sources'): del img.sources
                if hasattr(img, 'dsources'): del img.dsources
                if hasattr(img, 'gaussians'): del img.gaussians
                while 'gausfit' in img.completed_Ops:
                    img.completed_Ops.remove('gausfit')
                found = True
            if k in wavelet_atrous_opts:
                if hasattr(img, 'atrous_gaussians'): del img.atrous_gaussians
                if hasattr(img, 'islands'): del img.islands
                if hasattr(img, 'sources'): del img.sources
                if hasattr(img, 'dsources'): del img.dsources
                if hasattr(img, 'gaussians'): del img.gaussians
                while 'islands' in img.completed_Ops:
                    img.completed_Ops.remove('islands')
                while 'gausfit' in img.completed_Ops:
                    img.completed_Ops.remove('gausfit')
                while 'wavelet_atrous' in img.completed_Ops:
                    img.completed_Ops.remove('wavelet_atrous')
                found = True
            if k in shapelets_opts:
                while 'shapelets' in img.completed_Ops:
                    img.completed_Ops.remove('shapelets')
                found = True
            if k in gaul2srl_opts:
                while 'gaul2srl' in img.completed_Ops:
                    img.completed_Ops.remove('gaul2srl')
                found = True
            if k in spectralindex_opts:
                while 'spectralindex' in img.completed_Ops:
                    img.completed_Ops.remove('spectralindex')
                found = True
            if k in polarisation_opts:
                while 'polarisation' in img.completed_Ops:
                    img.completed_Ops.remove('polarisation')
                found = True
            if k in make_residimage_opts:
                if hasattr(img, 'resid_gaus_arr'):
                    del img.resid_gaus_arr
                    img.resid_gaus_arr = None  # set to init state
                if hasattr(img, 'model_gaus_arr'): del img.model_gaus_arr
                if hasattr(img, 'resid_shap_arr'): del img.resid_shap_arr
                if hasattr(img, 'model_shap_arr'): del img.model_shap_arr
                while 'make_residimage' in img.completed_Ops:
                    img.completed_Ops.remove('make_residimage')
                found = True
            if k in psf_vary_opts:
                while 'psf_vary' in img.completed_Ops:
                    img.completed_Ops.remove('psf_vary')
                found = True
            if not found:
                break

    # If nothing has changed, ask if user wants to re-run
    if not found and not re_run:
        prompt = "Analysis appears to be up-to-date. Force reprocessing (y/n)? "
        answ = raw_input_no_history(prompt)
        while answ.lower() not in  ['y', 'n', 'yes', 'no']:
            answ = raw_input_no_history(prompt)
        if answ.lower() in ['y', 'yes']:
            re_run = True # Force re-run
        else:
            return img, None

    # If a changed option is not in any of the above lists,
    # force a re-run of all Ops.
    if not found:
        img.completed_Ops = []
        if hasattr(img, 'use_io'): del img.use_io
        if hasattr(img, 'image_arr'): del img.image_arr
        if hasattr(img, 'mask_arr'): del img.mask_arr
        if hasattr(img, 'ch0_arr'): del img.ch0_arr
        if hasattr(img, 'rms_arr'): del img.rms_arr
        if hasattr(img, 'mean_arr'): del img.mean_arr
        if hasattr(img, 'rms_Q_arr'): del img.rms_Q_arr
        if hasattr(img, 'mean_Q_arr'): del img.mean_Q_arr
        if hasattr(img, 'rms_U_arr'): del img.rms_U_arr
        if hasattr(img, 'mean_U_arr'): del img.mean_U_arr
        if hasattr(img, 'rms_V_arr'): del img.rms_V_arr
        if hasattr(img, 'mean_V_arr'): del img.mean_V_arr
        if hasattr(img, 'islands'): del img.islands
        if hasattr(img, 'sources'): del img.sources
        if hasattr(img, 'dsources'): del img.dsources
        if hasattr(img, 'gaussians'): del img.gaussians
        if hasattr(img, 'atrous_gaussians'): del img.atrous_gaussians
        if hasattr(img, 'resid_gaus_arr'): del img.resid_gaus_arr
        if hasattr(img, 'model_gaus_arr'): del img.model_gaus_arr
        if hasattr(img, 'resid_shap_arr'): del img.resid_shap_arr
        if hasattr(img, 'model_shap_arr'): del img.model_shap_arr
        if hasattr(img, '_adapt_rms_isl_pos'): del img._adapt_rms_isl_pos
        return img, Op_chain

    while 'outlist' in img.completed_Ops:
        img.completed_Ops.remove('outlist')
    while 'cleanup' in img.completed_Ops:
        img.completed_Ops.remove('cleanup')
    for completed_Op in img.completed_Ops:
        if completed_Op in Op_names:
            Op_indx = Op_names.index(completed_Op)
            Op_names.pop(Op_indx)
            Op_chain.pop(Op_indx)

    return img, Op_chain

def load_pars(filename):
    """Load parameters from a save file or dictionary.

    If a file is given, it must be a pickled opts dictionary.

    filename - name of options file to load or a dictionary of opts.
    Returns None (and original error) if no file can be loaded successfully.
    """
    from .image import Image
    from . import mylogger
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
            print("--> Loaded parameters from file '" + filename + "'.")
            return timg, None
        except Exception as err:
            return None, err

def save_pars(img, savefile=None, quiet=False):
    """Save parameters to a file.

    The save file is a "pickled" opts dictionary.
    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    from . import tc
    import sys

    if savefile is None or savefile == '':
        basename = os.path.basename(img.opts.filename) + '.pybdsf.sav'
        savefile = os.path.join(img.basedir, basename)

    # convert opts to dictionary
    pars = img.opts.to_dict()
    output = open(savefile, 'wb')
    pickle.dump(pars, output, protocol=0)
    output.close()
    if not quiet:
        print("--> Saved parameters to file '" + savefile + "'.")

def list_pars(img, opts_list=None, banner=None, use_groups=True):
    """Lists all parameters for the Image object.

    opts_list - a list of the parameter names to list;
                if None, all parameters are used.
    banner - banner text to place at top of listing.
    use_groups - whether to use the group information for each
                 parameter.
    """
    from . import tc
    import sys

    # Get all options as a list sorted by name
    opts = img.opts.to_list()

    # Filter list
    if opts_list is not None:
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
    from .image import Image

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
            opt_names = temp_img.opts.get_names()
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
        if grp is not None and grp not in groups:
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
    from .image import Image
    import os
    from . import functions as func

    termy, termx = func.getTerminalSize() # note: returns row, col -> y, x
    minwidth = 28 # minimum width for parameter names and values

    # Define colors for output
    dc = '\033[1;34m' # Blue: non-default option text color
    ec = '\033[0;47m' # expandable option text color
    sc = '\033[0;32m' # Green: suboption text color
    nc = '\033[0m'    # normal text color
    ncb = '\033[1m'    # normal text color bold

    if banner is not None:
        print(banner)
    spcstr = ' ' * minwidth # spaces string for second or later lines
    infix = nc + ': ' + nc # infix character used to separate values from comments
    print('=' * termx) # division string for top of parameter listing
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
                    print(fmt % (parvalstr.ljust(minwidth), dt.ljust(44)))
                else:
                    print(nc + spcstr + '   %44s' % dt.ljust(44))
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
                        if k == 'beam_spectrum' and val is not None:
                            val = round_list_of_tuples(val)
                        if k == 'frequency_sp' and val is not None:
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
                            print(fmt % (parvalstr.ljust(minwidth-2), dt.ljust(44)))
                        else:
                            print(nc + spcstr + '   %44s' % dt.ljust(44))


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
    userinput = input(prompt)
    if userinput != '':
        readline.remove_history_item(readline.get_current_history_length()-1)
    return userinput


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
def export_image(img, outfile=None, img_format='fits', pad_image = False,
                 img_type='gaus_resid', mask_dilation=0, clobber=False):
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
        'island_mask' - Island mask image (0 = outside island, 1 = inside island)
    """
    import os
    from . import functions as func
    from .const import fwsig
    from . import mylogger

    mylog = mylogger.logging.getLogger("PyBDSF."+img.log+"ExportImage")

    # First some checking:
    if not 'gausfit' in img.completed_Ops and 'gaus' in img_type:
        print('\033[91mERROR\033[0m: Gaussians have not been fit. Please run process_image first.')
        return False
    elif not 'shapelets' in img.completed_Ops and 'shap' in img_type:
        print('\033[91mERROR\033[0m: Shapelets have not been fit. Please run process_image first.')
        return False
    elif not 'polarisation' in img.completed_Ops and 'pi' in img_type:
        print('\033[91mERROR\033[0m: Polarization properties have not been calculated. Please run process_image first.')
        return False
    elif not 'psf_vary' in img.completed_Ops and 'psf' in img_type:
        print('\033[91mERROR\033[0m: PSF variations have not been calculated. Please run process_image first.')
        return False
    elif not 'collapse' in img.completed_Ops and 'ch0' in img_type:
        print('\033[91mERROR\033[0m: ch0 image has not been calculated. Please run process_image first.')
        return False
    elif not 'rmsimage' in img.completed_Ops and ('rms' in img_type or 'mean' in img_type):
        print('\033[91mERROR\033[0m: Mean and rms maps have not been calculated. Please run process_image first.')
        return False
    elif not 'make_residimage' in img.completed_Ops and ('resid' in img_type or 'model' in img_type):
        print('\033[91mERROR\033[0m: Residual and model maps have not been calculated. Please run process_image first.')
        return False
    format = img_format.lower()
    if (format in ['fits', 'casa']) == False:
        print('\033[91mERROR\033[0m: img_format must be "fits" or "casa"')
        return False
    filename = outfile
    if filename is None or filename == '':
        filename = img.imagename + '_' + img_type + '.' + format
    if os.path.exists(filename) and clobber == False:
        print('\033[91mERROR\033[0m: File exists and clobber = False.')
        return False
    if format == 'fits':
        use_io = 'fits'
    if format == 'casa':
        use_io = 'rap'
    bdir = ''
    try:
        if img_type == 'ch0':
            func.write_image_to_file(use_io, filename,
                                     img.ch0_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'rms':
            func.write_image_to_file(use_io, filename,
                                     img.rms_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'mean':
            func.write_image_to_file(use_io, filename,
                                     img.mean_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'pi':
            func.write_image_to_file(use_io, filename,
                                     img.ch0_pi_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'psf_major':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_maj_arr*fwsig, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'psf_minor':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_min_arr*fwsig, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'psf_pa':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_pa_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'psf_ratio':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_ratio_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'psf_ratio_aper':
            func.write_image_to_file(use_io, filename,
                                     img.psf_vary_ratio_aper_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'gaus_resid':
            im = img.resid_gaus_arr
            func.write_image_to_file(use_io, filename,
                                     im, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'gaus_model':
            im = img.model_gaus_arr
            func.write_image_to_file(use_io, filename,
                                     im, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'shap_resid':
            func.write_image_to_file(use_io, filename,
                                     img.resid_shap_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'shap_model':
            func.write_image_to_file(use_io, filename,
                                     img.model_shap_arr, img, bdir, pad_image,
                                     clobber=clobber)
        elif img_type == 'island_mask':
            import numpy as N
            import scipy.ndimage as nd
            island_mask_bool = img.pyrank + 1 > 0
            if mask_dilation > 0:
                # Dilate the mask by specified number of iterations
                island_mask_bool = nd.binary_dilation(island_mask_bool,
                    iterations=mask_dilation)
                # Perform a binary closing to remove small holes/gaps. The
                # structure array is chosen to be about the size of the
                # beam (assuming a normally sampled psf), so that holes/gaps
                # smaller than the beam are removed.
                pbeam = int(round(img.beam2pix(img.beam)[0] * 1.5))
                island_mask_bool = nd.binary_closing(island_mask_bool,
                    structure=N.ones((pbeam, pbeam)))

            # Check for telescope, needed for CASA clean masks
            if img._telescope is None:
                print('\033[91mWARNING\033[0m: Telescope is unknown. Mask may not work correctly in CASA.')
            island_mask = N.array(island_mask_bool, dtype=N.float32)
            func.write_image_to_file(use_io, filename,
                                     island_mask, img, bdir, pad_image,
                                     clobber=clobber, is_mask=True)
        else:
            print("\n\033[91mERROR\033[0m: img_type not recognized.")
            return False
        if filename == 'SAMP':
            print('--> Image sent to SMAP hub')
        else:
            print('--> Wrote file ' + repr(filename))
            if use_io == 'rap':
                # remove the temporary fits file used as a casacore template
                import os
                os.remove(filename+'.fits')

        return True
    except RuntimeError as err:
        # Catch and log error
        mylog.error(str(err))

        # Re-throw error if the user is not in the interactive shell
        if img._is_interactive_shell:
            return False
        else:
            raise
    except KeyboardInterrupt:
        mylogger.userinfo(mylog, "\n\033[31;1mAborted\033[0m")
        if img._is_interactive_shell:
            return False
        else:
            raise


def write_catalog(img, outfile=None, format='bbs', srcroot=None, catalog_type='gaul',
               bbs_patches=None, incl_chan=False, incl_empty=False, clobber=False,
               force_output=False, correct_proj=True, bbs_patches_mask=None):
    """Write the Gaussian, source, or shapelet list to a file. Returns True if
    successful, False if not.

    filename - name of resulting file; if None, file is
               named automatically. If 'SAMP', table is sent to a samp hub
               (must be running already).
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
        "sagecal" - SAGECAL file (Gaussian list only)
    srcroot - root for source and patch names (BBS/ds9 only);
              if None, the srcroot is chosen automatically
    bbs_patches - type of patches to use:
        None - no patches
        "gaussian" - each Gaussian gets its own patch
        "single"   - all Gaussians are put into a single
                     patch
        "source"   - sources are grouped by source into patches
        "mask"     - use a Boolean mask to define the patches
    bbs_patches_mask - file name of mask file if bbs_patches="mask"
    incl_chan - Include fluxes for each channel?
    incl_empty - Include islands without any valid Gaussians (source list only)?
    sort_by - Property to sort output list by:
        "flux" - sort by total integrated flux, largest first
        "indx" - sort by Gaussian and island or source index, smallest first
    force_output - Force the creation of a catalog, even if it is empty
    correct_proj - Correct source parameters for image projection effects (BBS only)?
    clobber - Overwrite existing file?
    """
    from . import output

    # First some checking:
    if not 'gausfit' in img.completed_Ops:
        print('\033[91mERROR\033[0m: Image has not been fit. Please run process_image first.')
        return False
    if catalog_type == 'shap' and not 'shapelets' in img.completed_Ops:
        print('\033[91mERROR\033[0m: Image has not been decomposed into shapelets. Please run process_image first.')
        return False
    if catalog_type == 'srl' and not 'gaul2srl' in img.completed_Ops:
        print('\033[91mERROR\033[0m: Gaussians have not been grouped into sources. Please run process_image first.')
        return False
    format = format.lower()
    patch = bbs_patches
    filename = outfile
    if isinstance(patch, str):
        patch = patch.lower()
    if format not in ['fits', 'ascii', 'bbs', 'ds9', 'star',
                   'kvis', 'sagecal', 'csv', 'casabox']:
        print('\033[91mERROR\033[0m: format must be "fits", '\
            '"ascii", "ds9", "star", "kvis", "csv", "casabox", or "bbs"')
        return False
    if patch not in [None, 'gaussian', 'single', 'source', 'mask']:
        print('\033[91mERROR\033[0m: patch must be None, '\
                '"gaussian", "source", "single", or "mask"')
        return False
    if patch == 'mask':
        if bbs_patches_mask is None:
            print('\033[91mERROR\033[0m: if patch is "mask", bbs_patches_mask must be set to the file name of the mask file')
            return False
    if (catalog_type in ['gaul', 'srl', 'shap']) == False:
        print('\033[91mERROR\033[0m: catalog_type must be "gaul", '\
              '"srl", or "shap"')
        return False
    if catalog_type == 'shap' and format != 'fits':
        print("\033[91mERROR\033[0m: Only format = 'fits' is supported with shapelet output.")
        return False
    if (len(img.sources) == 0 and not incl_empty) or (len(img.sources) == 0 and len(img.dsources) == 0 and incl_empty):
        if not force_output:
            print('No sources were found in the image. Output file not written.')
            return False
    if filename == '':
        filename = None

    # Now go format by format and call appropriate function
    if filename == 'samp' or filename == 'SAMP':
        import tempfile
        from . import functions as func
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
        table_name = 'PyBDSF '+ catalog_type + ' table'
        if catalog_type == 'srl':
            img.samp_srl_table_url = 'file://' + os.path.abspath(tfile.name)
        if catalog_type == 'gaul':
            img.samp_gaul_table_url = 'file://' + os.path.abspath(tfile.name)
        func.send_fits_table(img.samp_client, img.samp_key, table_name, tfile.name)
        print('--> Table sent to SMAP hub')
        return True

    if format == 'fits':
        filename = output.write_fits_list(img, filename=filename,
                                             incl_chan=incl_chan, incl_empty=incl_empty,
                                             clobber=clobber, objtype=catalog_type)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber = False.')
            return False
        else:
            print('--> Wrote FITS file ' + repr(filename))
            return True
    if format == 'ascii' or format == 'csv':
        filename = output.write_ascii_list(img, filename=filename,
                                              incl_chan=incl_chan, incl_empty=incl_empty,
                                              sort_by='index', format = format,
                                              clobber=clobber, objtype=catalog_type)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber = False.')
            return False
        else:
            print('--> Wrote ASCII file ' + repr(filename))
            return True
    if format == 'bbs':
        if catalog_type != 'gaul':
            print("\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with BBS files.")
            return False
        filename = output.write_bbs_gaul(img, filename=filename,
                                            srcroot=srcroot, incl_empty=incl_empty,
                                            patch=patch, correct_proj=correct_proj,
                                            sort_by='flux',
                                            clobber=clobber)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber = False.')
            return False
        else:
            print('--> Wrote BBS sky model ' + repr(filename))
            return True
    if format == 'sagecal':
        if catalog_type != 'gaul':
            print("\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with Sagecal files.")
            return False
        filename = output.write_lsm_gaul(img, filename=filename,
                                            srcroot=srcroot, incl_empty=incl_empty,
                                            patch=patch,
                                            sort_by='flux',
                                            clobber=clobber)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber = False.')
            return False
        else:
            print('--> Wrote Sagecal lsm file ' + repr(filename))
            return True
    if format == 'ds9':
        filename = output.write_ds9_list(img, filename=filename,
                                            srcroot=srcroot, incl_empty=incl_empty,
                                            clobber=clobber, objtype=catalog_type)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber = False.')
            return False
        else:
            print('--> Wrote ds9 region file ' + repr(filename))
            return True
    if format == 'star':
        if catalog_type != 'gaul':
            print("\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with star files.")
            return False
        filename = output.write_star(img, filename=filename,
                                        clobber=clobber)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber = False.')
            return False
        else:
            print('--> Wrote AIPS STAR file ' + repr(filename))
            return True
    if format == 'kvis':
        if catalog_type != 'gaul':
            print("\033[91mERROR\033[0m: Only catalog_type = 'gaul' is supported with kvis files.")
            return False
        filename = output.write_kvis_ann(img, filename=filename,
                                            clobber=clobber)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber=False.')
            return False
        else:
            print('--> Wrote kvis file ' + repr(filename))
            return True
    if format == 'casabox':
        filename = output.write_casa_gaul(img, filename=filename,
                                      incl_empty=incl_empty, clobber=clobber)
        if filename is None:
            print('\033[91mERROR\033[0m: File exists and clobber=False.')
        else:
            print('--> Wrote CASA clean box file ' + filename)

def add_break_to_logfile(logfile):
    f = open(logfile, 'a')
    f.write('\n' + '='*72 + '\n')
    f.close()
