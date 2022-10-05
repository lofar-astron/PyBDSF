"""Module image.

Instances of class Image are a primary data-holders for all PyBDSF
operations. They store the image itself together with some meta-information
(such as headers), options for processing modules and all data generated during
processing. A few convenience methods are also defined here for interactive
use: to allow viewing and output of the most important data, to allow listing
and setting of options, and to allow re-processing of Images (these methods are
used by the interactive IPython shell made by pybdsf).

This module also defines class Op, which is used as a base class for all PyBDSF
operations.
"""
from __future__ import print_function
from __future__ import absolute_import

import numpy as N
from .opts import *


class Image(object):
    """Image is a primary data container for PyBDSF.

    All the run-time data (such as image data, mask, etc.)
    is stored here. A number of type-checked properties
    are defined for the most basic image attributes, such
    as image data, mask, header, user options.

    To allow transparent caching of large image data to disk,
    the image data must be stored in attributes ending in
    "_arr". Additionally, setting subarrays does not work
    using the attributes directly (e.g., img.ch0_arr[0:100,0:100]
    = 0.0 will not work). Instead, set the subarray values then set
    the attribute (e.g., ch0[0:100,0:100] = 0.0; img.ch0_arr = ch0).

    There is little sense in declaring all possible attributes
    right here as it will introduce unneeded dependencies
    between modules, thus most other attributes (like island lists,
    gaussian lists, etc) are inserted at run-time by the specific
    PyBDSF modules.
    """
    def __init__(self, opts):
        self._prev_opts = None
        self.extraparams = {}
        self.masked = False
        self.completed_Ops = []
        self.waveletimage = False
        self._pi = False
        self.do_cache = False
        self.bbspatchnum = 0
        self.blankpix = 0
        self.use_io = ''
        self.j = 0
        self.freq_pars = [0.0, 0.0, 0.0]
        self.filename = ''
        self.logfilename = ''
        self.resid_gaus_arr = None
        self._is_interactive_shell = False
        self.opts = Opts(opts)

    def __setstate__(self, state):
        """Needed for multiprocessing"""
        self.thresh_pix = state['thresh_pix']
        self.minpix_isl = state['minpix_isl']
        self.clipped_mean = state['clipped_mean']

    def __getstate__(self):
        """Needed for multiprocessing"""
        state = {}
        state['thresh_pix'] = self.thresh_pix
        state['minpix_isl'] = self.minpix_isl
        state['clipped_mean'] = self.clipped_mean
        return state

    def __getattribute__(self, name):
        from . import functions as func
        if name.endswith("_arr"):
            if self.do_cache:
                map_data = func.retrieve_map(self, name)
                if map_data is not None:
                    return map_data
                else:
                    return object.__getattribute__(self, name)
            else:
                return object.__getattribute__(self, name)
        else:
            return object.__getattribute__(self, name)

    def __setattr__(self, name, value):
        from . import functions as func
        if hasattr(self, 'do_cache'):
            if self.do_cache and name.endswith("_arr") and isinstance(value, N.ndarray):
                func.store_map(self, name, value)
            else:
                super(Image, self).__setattr__(name, value)
        else:
            super(Image, self).__setattr__(name, value)

    def __delattr__(self, name):
        from . import functions as func
        if self.do_cache and name.endswith("_arr"):
            func.del_map(self, name)
        else:
            super(Image, self).__delattr__(name)

    def get_map(self, map_name):
        """Returns requested map."""
        from . import functions as func
        if self.do_cache:
            map_data = func.retrieve_map(self, map_name)
        else:
            map_data = getattr(self, map_name)
        return map_data

    def put_map(self, map_name, map_data):
        """Stores requested map."""
        from . import functions as func
        if self.do_cache:
            func.store_map(self, map_name, map_data)
        else:
            setattr(self, map_name, map_data)

    def list_pars(self):
        """List parameter values."""
        from . import interface
        interface.list_pars(self)

    def set_pars(self, **kwargs):
        """Set parameter values."""
        from . import interface
        interface.set_pars(self, **kwargs)

    def process(self, **kwargs):
        """Process Image object"""
        from . import interface
        success = interface.process(self, **kwargs)
        return success

    def save_pars(self, savefile=None):
        """Save parameter values."""
        from . import interface
        interface.save_pars(self, savefile)

    def load_pars(self, loadfile=None):
        """Load parameter values."""
        from . import interface
        import os
        if loadfile is None or loadfile == '':
            loadfile = self.opts.filename + '.pybdsf.sav'
        if os.path.exists(loadfile):
            timg, err = interface.load_pars(loadfile)
            if timg is not None:
                orig_filename = self.opts.filename
                self.opts = timg.opts
                self.opts.filename = orig_filename # reset filename to original
            else:
                if self._is_interactive_shell:
                    print("\n\033[31;1mERROR\033[0m: '"+\
                    loadfile+"' is not a valid parameter save file.")
                else:
                    raise RuntimeError(str(err))
        else:
            if self._is_interactive_shell:
                print("\n\033[31;1mERROR\033[0m: File '"+\
                loadfile+"' not found.")
            else:
                raise RuntimeError('File not found')

    def show_fit(self, **kwargs):
        """Show results of the fit."""
        from . import plotresults
        if not hasattr(self, 'nisl'):
            print('Image has not been processed. Please run process_image first.')
            return False
        plotresults.plotresults(self, **kwargs)
        return True

    def export_image(self, **kwargs):
        """Export an internal image to a file."""
        from . import interface
        try:
            result = interface.export_image(self, **kwargs)
            return result
        except RuntimeError as err:
            if self._is_interactive_shell:
                print("\n\033[31;1mERROR\033[0m: " + str(err))
            else:
                raise RuntimeError(str(err))

    def write_catalog(self, **kwargs):
        """Write the Gaussian, source, or shapelet list to a file"""
        from . import interface
        try:
            result = interface.write_catalog(self, **kwargs)
            return result
        except RuntimeError as err:
            if self._is_interactive_shell:
                print("\n\033[31;1mERROR\033[0m: " + str(err))
            else:
                raise RuntimeError(str(err))


class Op(object):
    """Common base class for all PyBDSF operations.

    At the moment this class is empty and only defines placeholder
    for method __call__, which should be redefined in all derived
    classes.
    """
    def __call__(self, img):
        raise NotImplementedError("This method should be redefined")
