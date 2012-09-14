"""Module image.

Instances of class Image are a primary data-holders for all PyBDSM
operations. They store the image itself together with some meta-information
(such as headers), options for processing modules and all data generated during
processing. A few convenience methods are also defined here for interactive
use: to allow viewing and output of the most important data, to allow listing
and setting of options, and to allow re-processing of Images (these methods are
used by the interactive IPython shell made by pybdsm.py).

This module also defines class Op, which is used as a base class for all PyBDSM
operations.
"""

import numpy as N
from opts import *

class Image(object):
    """Image is a primary data container for PyBDSM.

    All the run-time data (such as image data, mask, etc.)
    is stored here. A number of type-checked properties
    are defined for the most basic image attributes, such
    as image data, mask, header, user options.

    There is little sense in declaring all possible attributes
    right here as it will introduce unneeded dependencies
    between modules, thus most other attributes (like island lists,
    gaussian lists, etc) are inserted at run-time by the specific
    PyBDSM modules.
    """
    opts   = Instance(Opts, doc="User options")
    image  = NArray(doc="Image data, Stokes I")
    ch0    = NArray(doc="Channel-collapsed image data, Stokes I")
    ch0_Q  = NArray(doc="Channel-collapsed image data, Stokes Q")
    ch0_U  = NArray(doc="Channel-collapsed image data, Stokes U")
    ch0_V  = NArray(doc="Channel-collapsed image data, Stokes V")
    header = Any(doc="Image header")
    mask   = NArray(doc="Image mask (if present and attribute masked is set)")
    masked = Bool(False, doc="Flag if mask is present")
    basedir = String('DUMMY', doc="Base directory for output files")
    completed_Ops = List(String(), doc="List of completed operations")
    _is_interactive_shell = Bool(False, doc="PyBDSM is being used in the interactive shell")


    def __init__(self, opts):
        self.opts = Opts(opts)
        self.extraparams = {}

    def __setstate__(self, state):
        """Needed for multiprocessing"""
        self.pixel_beamarea = state['pixel_beamarea']
        self.pixel_beam = state['pixel_beam']
        self.thresh_pix = state['thresh_pix']
        self.minpix_isl = state['minpix_isl']
        self.clipped_mean = state['clipped_mean']

    def __getstate__(self):
        """Needed for multiprocessing"""
        state = {}
        state['pixel_beamarea'] = self.pixel_beamarea
        state['pixel_beam'] = self.pixel_beam
        state['thresh_pix'] = self.thresh_pix
        state['minpix_isl'] = self.minpix_isl
        state['clipped_mean'] = self.clipped_mean
        return state

    def list_pars(self):
        """List parameter values."""
        import interface
        interface.list_pars(self)

    def set_pars(self, **kwargs):
        """Set parameter values."""
        import interface
        interface.set_pars(self, **kwargs)

    def process(self, **kwargs):
        """Process Image object"""
        import interface
        success = interface.process(self, **kwargs)
        return success

    def save_pars(self, savefile=None):
        """Save parameter values."""
        import interface
        interface.save_pars(self, savefile)

    def load_pars(self, loadfile=None):
        """Load parameter values."""
        import interface
        import os
        if loadfile == None or loadfile == '':
            loadfile = self.opts.filename + '.pybdsm.sav'
        if os.path.exists(loadfile):
            timg, err = interface.load_pars(loadfile)
            if timg != None:
                orig_filename = self.opts.filename
                self.opts = timg.opts
                self.opts.filename = orig_filename # reset filename to original
            else:
                if self._is_interactive_shell:
                    print "\n\033[31;1mERROR\033[0m: '"+\
                    loadfile+"' is not a valid parameter save file."
                else:
                    raise RuntimeError(str(err))
        else:
            if self._is_interactive_shell:
                print "\n\033[31;1mERROR\033[0m: File '"+\
                loadfile+"' not found."
            else:
                raise RuntimeError('File not found')

    def show_fit(self, **kwargs):
        """Show results of the fit."""
        import plotresults
        if not hasattr(self, 'nisl'):
            print 'Image has not been processed. Please run process_image first.'
            return False
        plotresults.plotresults(self, **kwargs)
        return True

    def export_image(self, **kwargs):
          """Export an internal image to a file."""
          import interface
          interface.export_image(self, **kwargs)

    def write_catalog(self, **kwargs):
        """Write the Gaussian, source, or shapelet list to a file"""
        import interface
        interface.write_catalog(self, **kwargs)
    write_gaul = write_catalog # for legacy scripts


class Op(object):
    """Common base class for all PyBDSM operations.

    At the moment this class is empty and only defines placeholder
    for method __call__, which should be redefined in all derived
    classes.
    """
    def __call__(self, img):
        raise NotImplementedError("This method should be redefined")
