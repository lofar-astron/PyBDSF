import bdsf
import sys
import filecmp
from compare_outputs import compare_results

# Process the image
img = bdsf.process_image('tbdsf_process_image.in', ncores=2, output_all=True)

# List of operations that must have been done on `img`.
operations = [
    'readimage', 'collapse', 'preprocess', 'rmsimage', 'threshold',
    'islands', 'gausfit', 'gaul2srl', 'make_residimage', 'wavelet_atrous',
    'shapelets', 'spectralindex', 'polarisation', 'psf_vary', 'cleanup'
]

# Check if the outputs agree with the reference ones. The arguments used here
# are suitable for CI jobs, where libraries, architecture, etc. may vary from
# those used to generate the reference outputs. If the check is done against a
# reference generated on the same system (e.g., for testing by hand), all
# options can be enabled
outputs_agree = compare_results(
    filecmp.dircmp("tbdsf_process_image.in_fits_pybdsf", "reference_outputs"),
    1e-2,  # rtol: a value of 1e-2 works well for CI jobs; 1e-3 can be used otherwise
    2,  # verbosity
    check_images=True,  # compare the rms values (within rtol) of images
    check_catalogs=False,  # compare entries (within rtol) in FITS catalogs; not suitable for CI jobs
    check_other=False  # compare entries in text files; not suitable for CI jobs
)

# Return exit status 0 if everything went fine, otherwise return 1.
if img and all(oper in img.completed_Ops for oper in operations) and outputs_agree:
    sys.exit(0)
else:
    sys.exit(1)
