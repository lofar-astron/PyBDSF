"""Interactive PyBDSM shell.

This module initializes the interactive PyBDSM shell, which is a customized
IPython enviroment. It should be called from the terminal prompt using the
"pybdsm" shell script or as "python pybdsm.py".
"""
import lofar.bdsm
from lofar.bdsm.image import Image
import pydoc
import sys
import inspect


###############################################################################
# Functions needed only in the custom IPython shell are defined here. Other
# functions used by both the custom shell and normal Python or IPython
# environments are defined in interface.py.
#
# Before starting the IPython shell, we need to define all the functions and
# variables that we want in the namespace. Note that we adopt the convention
# for this UI of using lines of 72 characters max for doc strings and the
# start-up banner. However, the parameter list will fill the entire available
# terminal width to consume as few vertical lines as possible.
global _img
_img = Image({'filename':''})
_img._is_interactive_shell = True
T = True
F = False
true = True
false = False

def inp(cur_cmd=None):
    """List inputs for current task.

    If a task is given as an argument, inp sets the current task
    to the given task. If no task is given, inp lists the parameters
    of the current task.
    """
    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    if cur_cmd != None:
        if not hasattr(cur_cmd, 'arg_list'):
            print '\033[31;1mERROR\033[0m: not a valid task'
            return
        _set_current_cmd(cur_cmd)
    else:
        if not hasattr(_img, '_current_cmd'):
            print '\033[31;1mERROR\033[0m: no task is set'
            return
    lofar.bdsm.interface.list_pars(_img, opts_list=_img._current_cmd_arg_list,
                             banner=_img._current_cmd_desc,
                             use_groups=_img._current_cmd_use_groups)


def go(cur_cmd=None):
    """Executes the current task.

    If a task is given as an argument, go executes the given task,
    even if it is not the current task. The current task is not
    changed in this case.
    """
    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    if cur_cmd == None:
        if not hasattr(_img, '_current_cmd'):
            print '\033[31;1mERROR\033[0m: no task is set'
            return
        cur_cmd = _img._current_cmd
    if not hasattr(cur_cmd, 'arg_list'):
        print '\033[31;1mERROR\033[0m: not a valid task'
        return
    cur_cmd()


def default(cur_cmd=None):
    """Resets all parameters for a task to their default values.

    If a task name is given (e.g., "default show_fit"), the
    parameters for that task are reset. If no task name is
    given, the parameters of the current task are reset.
    """
    global _img
    if cur_cmd == None:
        if not hasattr(_img, '_current_cmd'):
            print '\033[31;1mERROR\033[0m: no task is set'
            return
        cur_cmd = _img._current_cmd

    if hasattr(cur_cmd, 'arg_list'):
        opts_list = cur_cmd.arg_list
    else:
        print '\033[31;1mERROR\033[0m: not a valid task'
        return
    _img.opts.set_default(opts_list)
    _replace_vals_in_namespace(opt_names=opts_list)


def tget(filename=None):
    """Load processing parameters from a parameter save file.

    A file name may be given (e.g., "tget 'savefile.sav'"), in which case the
    parameters are loaded from the file specified. If no file name is given,
    the parameters are loaded from the file 'pybdsm.last' if it exists.

    Normally, the save file is created by the tput command (try "help tput"
    for more info).

    The save file is a "pickled" python dictionary which can be loaded into
    python and edited by hand. See the pickle module for more information.
    Below is an example of how to edit a save file by hand:

      BDSM [1]: import pickle
      BDSM [2]: savefile = open('savefile.sav', 'w')
      BDSM [3]: pars = pickle.load(savefile)
      BDSM [4]: pars['rms_box'] = (80, 20)  --> change rms_box parameter
      BDSM [5]: pickle.dump(pars, savefile) --> save changes

    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    import os

    global _img

    # Check whether user has given a task name as input (as done in casapy).
    # If so, reset filename to None.
    if hasattr(filename, 'arg_list'):
        filename = None

    if filename == None or filename == '':
        if os.path.isfile('pybdsm.last'):
            filename = 'pybdsm.last'
        else:
            print '\033[31;1mERROR\033[0m: No file name given and '\
                  '"pybdsm.last" not found.\nPlease specify a file to load.'
            return

    if os.path.isfile(filename):
        try:
            pkl_file = open(filename, 'rb')
            pars = pickle.load(pkl_file)
            pkl_file.close()
            _img.opts.set_opts(pars)
            _replace_vals_in_namespace()
            print "--> Loaded parameters from file '" + filename + "'."
        except:
            print "\033[31;1mERROR\033[0m: Could not read file '" + \
                  filename + "'."
    else:
        print "\033[31;1mERROR\033[0m: File '" + filename + "' not found."


def tput(filename=None, quiet=False):
    """Save processing parameters to a file.

    A file name may be given (e.g., "tput 'savefile.sav'"), in which case the
    parameters are saved to the file specified. If no file name is given, the
    parameters are saved to the file 'pybdsm.last'. The saved parameters can
    be loaded using the tget command (try "help tget" for more info).

    The save file is a "pickled" python dictionary which can be loaded into
    python and edited by hand. See the pickle module for more information.
    Below is an example of how to edit a save file by hand:

      BDSM [1]: import pickle
      BDSM [2]: savefile = open('savefile.sav', 'w')
      BDSM [3]: pars = pickle.load(savefile)
      BDSM [4]: pars['rms_box'] = (80, 20)  --> change rms_box parameter
      BDSM [5]: pickle.dump(pars, savefile) --> save changes

    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    if filename == None or filename == '':
        filename = 'pybdsm.last'

    # convert opts to dictionary
    pars = _img.opts.to_dict()
    output = open(filename, 'wb')
    pickle.dump(pars, output)
    output.close()
    if not quiet:
        print "--> Saved parameters to file '" + filename + "'."


def _set_pars_from_prompt():
    """Gets parameters and value and stores them in _img.

    To do this, we extract all the valid parameter names
    and values from the f_locals directory. Then, use
    set_pars() to set them all.

    Returns True if successful, False if not.
    """
    global _img
    f = sys._getframe(len(inspect.stack())-1)
    f_dict = f.f_locals

    # Check through all possible options and
    # build options dictionary
    opts = _img.opts.to_dict()
    user_entered_opts = {}
    for k, v in opts.iteritems():
        if k in f_dict:
            if f_dict[k] == '':
                # Set option to default value in _img and namespace
                _img.opts.set_default(k)
                f_dict[k] = _img.opts.__getattribute__(k)
            user_entered_opts.update({k: f_dict[k]})

    # Finally, set the options
    try:
        _img.opts.set_opts(user_entered_opts)
        return True
    except RuntimeError, err:
        # If an opt fails to set, replace its value in the namespace
        # with its current value in _img. Then print error so user knows.
        err_msg = str(err)
        err_msg_trim = err_msg.split('(')[0]
        indx1 = err_msg_trim.find('"') + 1
        indx2 = err_msg_trim.find('"', indx1)
        k = err_msg_trim[indx1:indx2]
        orig_opt_val = opts[k]
        f_dict[k] = orig_opt_val
        print '\033[31;1mERROR\033[0m: ' + err_msg_trim + \
              '\nResetting to previous value.'
        return False


def _replace_vals_in_namespace(opt_names=None):
    """Replaces opt values in the namespace with the ones in _img.

    opt_names - list of option names to replace (can be string if only one)
    """
    global _img
    f = sys._getframe(len(inspect.stack())-1)
    f_dict = f.f_locals
    if opt_names == None:
        opt_names = _img.opts.get_names()
    if isinstance(opt_names, str):
        opt_names = [opt_names]
    for opt_name in opt_names:
        if opt_name in f_dict:
            f_dict[opt_name] = _img.opts.__getattribute__(opt_name)


def _set_current_cmd(cmd):
    """Sets information about current command in img.

    This function is used to emulate a casapy interface.

    """
    global _img
    cmd_name = cmd.__name__
    doc = cmd.__doc__
    _img._current_cmd = cmd
    _img._current_cmd_name = cmd_name
    _img._current_cmd_desc = cmd_name.upper() + ': ' + doc.split('\n')[0]
    _img._current_cmd_arg_list = cmd.arg_list
    _img._current_cmd_use_groups = cmd.use_groups


###############################################################################
# Next, we define the tasks such that they may be called directly by
# the user if so desired. These functions simply pass on the user-
# specified arguments to the appropriate Image method. Here we also
# define the detailed doc strings used by help, and, after each task
# definition, we define its list of arguments and whether it should
# use the opts 'group' attribute, both needed when inp is called. If
# a new parameter is added to a task, it needs to be added to opts.py
# and to the list of arguments for the task below (the "arg_list")
# attribute.
def process_image(**kwargs):
    """Find and measure sources in an image.

    There are many possible parameters and options for process_image. Use
    "inp process_image" to list them. To get more information about a
    parameter, use help. E.g.,

    > help 'rms_box'

    When process_image is executed, PyBDSM performs the following steps in
    order:

    1. Reads in the image.

    2. Calculates basic statistics of the image and stores them in the Image
    object. Calculates sensible values of processing parameters and stores
    them. First calculates mean and rms, with and without (3-sigma)
    clipping, min and max pixel and values, solid angle. Hereafter, rms
    indicates the 3-sigma clipped measure. Next, the number of beams per
    source is calculated (see help on algorithms for details), using a
    sensible estimate of boxsize and stepsize (which can be set using the
    rms_box parameter). Finally, the thresholds are set. They can either be
    hard-thresholded (by the user or set as 5-sigma for pixel threshold and
    3-sigma for island boundaries internally) or can be calculated using the
    False Detection Rate (FDR) method using an user defined value for
    alpha. If the user does not specify whether hard thresholding or FDR
    should be applied, one or the other is chosen internally based on the
    ratio of expected false pixels and true pixels (the choice is written
    out in the log file).

    3. Calculates rms image. 3-sigma clipped rms and mean are calculated
    inside boxes of size boxsize in steps of stepsize. Intermediate values
    are calculated using bilinear interpolation (it was seen that bicubic
    spline did not yield appreciably better results but is also
    available). Depending on the resulting statistics (see help on
    algorithms for details), we either adopt the rms image or a constant rms
    in the following analysis.

    4. Identifies islands of contiguous emission. First all pixels greater
    than the pixel threshold are identified (and sorted by descending flux
    order). Next, starting from each of these pixels, all contiguous pixels
    (defined by 8-connectivity, i.e., the surrounding eight pixels) higher
    than the island boundary threshold are identified as belonging to one
    island, accounting properly for overlaps of islands.

    5. Fit multiple gaussians and/or shapelets to each island. For each
    island, the subimages of emission and rms are cut out. The number of
    multiple gaussians to be fit can be determined by three different
    methods (see help on algorithms for details). With initial guesses
    corresponding to these peaks, gaussians are simultaneously fit to the
    island using the Levenberg-Marqhardt algorithm. Sensible criteria for bad
    solutions are defined. If multiple gaussians are fit and one of them is
    a bad solution then the number of gaussians is decreased by one and fit
    again, till all solutions in the island are good (or zero in number, in
    which case its flagged). After the final fit to the island, the
    deconvolved size is computed assuming the theoretical beam and the
    statistics in the source area and in the island are computed and
    stored. Errors on each of the fitted parameters are computed using the
    formulae in Condon (1997). Finally all good solutions are written into
    the gaussian catalog as an ascii and binary file. If shapelets are
    required, the program calculates optimal nmax, beta and the centre, and
    stores these and the shapelet coefficients in a file.

    """
    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    # Save current command, as it might be overwritten when process
    # is called by the user directly and is not the current command.
    cur_cmd = _img._current_cmd

    # Run process. Note that process automatically picks up options
    # from the Image object, so we don't need to get_task_kwargs as
    # we do for the other tasks.
    success = _img.process(**kwargs)

    # Now restore parameters and save to pybdsm.last
    if success:
        _set_current_cmd(cur_cmd)
        tput(quiet=True)

task_list = _img.opts.get_names()
process_image.arg_list = task_list
process_image.use_groups = True


def show_fit(**kwargs):
    """Show results of fit.

    Selected plots are displayed to give the user a quick overview of the
    results of the fit. The plots may be zoomed, saved to a file, etc. using
    the controls at the bottom of the plot window.

    In addition, the following commands are available:
      Press "i" ........ : Get integrated flux densities and mean rms
                           values for the visible portion of the image
      Press "m" ........ : Change min and max scaling values
      Press "n" ........ : Show / hide island IDs
      Press "0" ........ : Reset scaling to default
      Press "c" ........ : Change source for SED plot
      Click Gaussian ... : Print Gaussian and source IDs (zoom_rect mode,
                           toggled with the "zoom" button and indicated in
                           the lower right corner, must be off)
                           The SED plot will also show the chosen source.

    Parameters: ch0_image, rms_image, mean_image, ch0_islands,
                gresid_image, sresid_image, gmodel_image,
                smodel_image, source_seds, ch0_flagged, pi_image,
                psf_major, psf_minor, psf_pa, broadcast

    For more information about a parameter, use help.  E.g.,
      > help 'ch0_image'

    """
    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    img_kwargs = _get_task_kwargs(show_fit)
    for k in kwargs:
        # If user enters an argument, use it instead of
        # that in _img
        img_kwargs[k] = kwargs[k]
    try:
        success = _img.show_fit(**img_kwargs)
        if success:
            tput(quiet=True)
    except KeyboardInterrupt:
        print "\n\033[31;1mAborted\033[0m"

show_fit.arg_list = ['ch0_image', 'rms_image', 'mean_image', 'ch0_islands',
                     'gresid_image', 'sresid_image', 'gmodel_image',
                     'smodel_image', 'source_seds', 'ch0_flagged', 'pi_image',
                     'psf_major', 'psf_minor', 'psf_pa', 'broadcast']
show_fit.use_groups = False


def write_catalog(**kwargs):
    """Write the Gaussian, source, or shapelet list to a file.

    The lists can be written in a number of formats. The information
    included in the output file varies with the format used. Use
    "help 'format'" for more information.

    Parameters: outfile, format, srcroot, bbs_patches, incl_wavelet, clobber,
                catalog_type, incl_empty

    For more information about a parameter, use help.  E.g.,
      > help 'bbs_patches'

    """
    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    img_kwargs = _get_task_kwargs(write_catalog)
    for k in kwargs:
        # If user enters an argument, use it instead of
        # that in _img
        img_kwargs[k] = kwargs[k]
    try:
        success = _img.write_catalog(**img_kwargs)
        if success:
            tput(quiet=True)
    except KeyboardInterrupt:
        print "\n\033[31;1mAborted\033[0m"

write_catalog.arg_list = ['bbs_patches', 'format', 'outfile', 'srcroot',
                          'incl_chan', 'clobber', 'catalog_type', 'incl_empty']
write_catalog.use_groups = False


def export_image(**kwargs):
    """Write an image to disk.

    Parameters: filename, img_type, img_format, incl_wavelet, clobber

    For more information about a parameter, use help.  E.g.,
      > help 'img_type'

    """
    global _img
    success = _set_pars_from_prompt()
    if not success:
        return
    img_kwargs = _get_task_kwargs(export_image)
    for k in kwargs:
        # If user enters an argument, use it instead of
        # that in _img
        img_kwargs[k] = kwargs[k]
    try:
        success = _img.export_image(**img_kwargs)
        if success:
            tput(quiet=True)
    except KeyboardInterrupt:
        print "\n\033[31;1mAborted\033[0m"

export_image.arg_list = ['outfile', 'img_type', 'img_format',
                         'clobber']
export_image.use_groups = False


def _get_task_kwargs(task):
    """Returns dictionary of keyword arguments from _img for the given task."""
    global _img
    arg_list = task.arg_list
    kwargs = {}
    for a in arg_list:
        kwargs.update({a: _img.opts.__getattribute__(a)})
    return kwargs


###############################################################################
# Customize the help system for PyBDSM. The user can type "help task" to get
# help on a task (it prints the doc string) or "help 'opt'" to get help on
# a option (it prints the doc string defined in opts.py).
class bdsmDocHelper(pydoc.Helper):
    def help(self, request):
        global _img
        topbar = '_' * 72 + '\n' # 72-character divider
        if hasattr(request, '__name__'):
            pydoc.pager(topbar + 'Help on ' + pydoc.text.bold(request.__name__)
                        + ':\n\n' + pydoc.getdoc(request))
        else:
            opts = _img.opts.__class__.__dict__
            try:
                opt = opts[request]
                desc_list = str(opt.doc()).split('\n')
                desc = '\n\n'.join(desc_list)
                default_val = opt._default
                if isinstance(default_val, str):
                    valstr = "'" + default_val + "'"
                else:
                    valstr = str(default_val)
                default_val_text = 'Default value: ' + valstr
                if opt.group() != None and opt.group() != 'hidden':
                    group_text = '\nBelongs to group: ' + opt.group()
                else:
                    group_text = ''
                desc_text = lofar.bdsm.interface.wrap(desc, 72)
                desc_text = '\n'.join(desc_text)
                pydoc.pager(topbar + 'Help on the ' + pydoc.text.bold(request)
                            + ' parameter:\n\n' + default_val_text
                            + group_text
                            + '\n\n' + desc_text)
            except(KeyError):
                print "Parameter '" + request + "' not recognized."
pydoc.help = bdsmDocHelper(sys.stdin, sys.stdout)


###############################################################################
# Now run the IPython shell with this namespace and a customized autocompleter.
# The custom autocompleter is below. It adds task, command, and option names and
# a few common values to ipython's autocompleter. It also adds files in the
# local directory when they might be needed (but only if the user has started
# to enter a string -- this behavior is to help avoid entering filenames as
# non-strings; this is also done for the help autocomplete).
def _opts_completer(self, event):
    """ Returns a list of strings with possible completions."""
    import os
    import glob
    from lofar.bdsm.image import Image
    img = Image({'filename':''})
    opts = img.opts.get_names()

    # Split the command entered by user when TAB was pressed
    # and check for up to three components (from e.g. "par = val",
    # which gives cmd1 = "par", cmd2 = "=", and cmd3 = "val")
    cmd1 = (event.line).rsplit(None)[0]
    if len((event.line).rsplit(None)) > 1:
        cmd2 = (event.line).rsplit(None)[1]
    else:
        cmd2 = ''
    if len((event.line).rsplit(None)) > 2:
        cmd3 = (event.line).rsplit(None)[2]
    else:
        cmd3 = ''

    # First, check to see if user has entered a parameter name
    # and an equals sign. If so, check parameter type. If Enum
    # or Option, match only to the allowable values.
    # Allowable values are available from v._type.values if v is
    # type Enum (v has no attribute _type.values if not).
    if "=" in cmd1 or "=" in cmd2:
        par_vals = []
        if "=" in cmd1:
            cmd3 = cmd1.split('=')[1]
            cmd1 = cmd1.split('=')[0]
        if cmd1 in opts:
            from lofar.bdsm.tc import tcEnum, tcOption
            v = img.opts.__class__.__dict__[cmd1]
            partype = v._type
            if isinstance(partype, tcOption):
                par_vals = ['None']
            elif isinstance(partype, tcEnum):
                if ('"' in cmd2 or "'" in cmd2 or
                    '"' in cmd3 or "'" in cmd3):
                    par_vals = v._type.values
                    if not isinstance(par_vals, list):
                        par_vals = list(par_vals)
                    if None in par_vals:
                        # Remove None from list
                        pindx = par_vals.index(None)
                        par_vals.pop(pindx)
                else:
                    if None in v._type.values:
                        par_vals.append('None')
                    if True in v._type.values:
                        par_vals.append('True')
                    if False in v._type.values:
                        par_vals.append('False')
            elif v._default == True or v._default == False:
                par_vals = ['True', 'False']
        if cmd1 == 'filename' or cmd1 == 'outfile':
            if ('"' in cmd2 or "'" in cmd2 or
                    '"' in cmd3 or "'" in cmd3):
                # Also add files in current directory
                found = [f.replace('\\','/') for f in glob.glob('*')]
                if len(found) > 0:
                    for fnd in found:
                        par_vals.append(fnd)
        return par_vals
    elif cmd1 == 'inp' or cmd1 == 'go':
        # Match task names only
        cmds = ['process_image', 'write_catalog', 'export_image', 'show_fit']
        return cmds
    elif cmd1 == 'cd' or cmd1 == 'tput' or cmd1 == 'tget' or '!' in cmd1:
        # Match to files in current directory (force use of ' or " with
        # tput and tget, as filename must be a string).
        files = []
        found = [f.replace('\\','/') for f in glob.glob('*')]
        if len(found) > 0:
            for fnd in found:
                files.append(fnd)
        if cmd1 == 'tput' or cmd1 == 'tget' and not ('"' in cmd2 or
                                                     "'" in cmd2):
            # User has not (yet) started to enter a string, so don't
            # return filenames
            return []
        return files
    elif cmd1 == 'help':
        if '"' in cmd2 or "'" in cmd2:
            # User has started to enter a string:
            # Match to parameter names, as they must be strings
            par_vals = opts
            return par_vals
        else:
            # User has not started to enter a string:
            # Match to commands + tasks only
            cmds = ['process_image', 'write_catalog', 'export_image',
                    'show_fit', 'go', 'inp', 'tget', 'tput', 'default',
                    'changelog']
            return cmds
    else:
        # Match to parameter, task, and command names only
        # Add command names
        opts.append('inp')
        opts.append('go')
        opts.append('tget')
        opts.append('tput')
        opts.append('default')
        opts.append('help')

        # Add task names
        opts.append('process_image')
        opts.append('show_fit')
        opts.append('write_catalog')
        opts.append('export_image')
        return opts

# Define the welcome banner to print on startup. Also check if there is a newer
# version on the STRW ftp server. If there is, print a message to the user
# asking them to update.
from lofar.bdsm._version import __version__, __revision__, changelog

# Query the STRW FTP server. Tar file must be named "PyBDSM-version#.tar.gz":
#   e.g., "PyBDSM-1.3.1.tar.gz".
# Check whether called from the LOFAR CEPI/II. If so, skip check.
import os
aps_local_val = os.environ.get('APS_LOCAL')
check_for_newer = True
if aps_local_val == None and check_for_newer:
    try:
        import ftplib
        from distutils.version import StrictVersion
        f = ftplib.FTP()
        f.connect("ftp.strw.leidenuniv.nl")
        f.login()
        file_list = []
        file_list = f.nlst('pub/rafferty/PyBDSM')
        f.close()
        ftp_version = ''
        for file in file_list:
            if 'PyBDSM' in file and '.tar.gz' in file:
                ver_start_indx = file.find('-') + 1
                ver_end_indx = file.find('.tar.gz')
                ftp_version = file[ver_start_indx:ver_end_indx]
        if ftp_version == '':
            # No matching files found, continue without message
            pass
        elif StrictVersion(__version__) < StrictVersion(ftp_version):
            print '\n' + '*' * 72
            print "There appears to be a newer version of PyBDSM available at:"
            print "    ftp://ftp.strw.leidenuniv.nl/pub/rafferty/PyBDSM/"
            print "Please consider updating your installation"
            print '*' * 72
    except:
        pass

divider1 = '=' * 72 + '\n'
divider2 = '_' * 72 + '\n'
banner = '\nPyBDSM version ' + __version__ + ' (LOFAR revision ' + \
         __revision__ + ')\n'\
+ divider1 + 'PyBDSM commands\n'\
'  inp task ............ : Set current task and list parameters\n'\
"  par = val ........... : Set a parameter (par = '' sets it to default)\n"\
'                          Autocomplete (with TAB) works for par and val\n'\
'  go .................. : Run the current task\n'\
'  default ............. : Set current task parameters to default values\n'\
"  tput ................ : Save parameter values\n"\
"  tget ................ : Load parameter values\n"\
'PyBDSM tasks\n'\
'  process_image ....... : Process an image: find sources, etc.\n'\
'  show_fit ............ : Show the results of a fit\n'\
'  write_catalog ....... : Write out list of sources to a file\n'\
'  export_image ........ : Write residual/model/rms/mean image to a file\n'\
'PyBDSM help\n'\
'  help command/task ... : Get help on a command or task\n'\
'                          (e.g., help process_image)\n'\
"  help 'par' .......... : Get help on a parameter (e.g., help 'rms_box')\n"\
'  help changelog ...... : See list of recent changes\n'\
+ divider2

# Go ahead and set the current task to process_image, so that the user does not
# need to enter "inp process_image" as the first step (the first task needed
# after startup will almost always be process_image).
_set_current_cmd(process_image)

# Now start the ipython shell. Due to (non-backward-compatible) changes in
# ipython with version 0.11, we must support both versions until 0.11 or
# greater is in common use.
try:
    # IPython >= 0.11
    from IPython.frontend.terminal.embed import InteractiveShellEmbed
    from IPython.config.loader import Config
    from IPython import __version__ as ipython_version
    cfg = Config()
    prompt_config = cfg.PromptManager
    if ipython_version == '0.11':
        cfg.InteractiveShellEmbed.prompt_in1 = "BDSM [\#]: "
    else:
        prompt_config.in_template = "BDSM [\#]: "
    cfg.InteractiveShellEmbed.autocall = 2
    ipshell = InteractiveShellEmbed(config=cfg, banner1=banner,
                                    user_ns=locals())
    ipshell.set_hook('complete_command', _opts_completer, re_key = '.*')
except ImportError:
    # IPython < 0.11
    from IPython.Shell import IPShellEmbed
    argv = ['-prompt_in1','BDSM [\#]: ','-autocall','2']
    ipshell = IPShellEmbed(argv=argv, banner=banner, user_ns=locals())
    ipshell.IP.set_hook('complete_command', _opts_completer, re_key = '.*')
ipshell()


