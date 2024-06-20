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

# Check if the outputs agree with the reference ones
outputs_agree = compare_results(
    filecmp.dircmp("tbdsf_process_image.in_fits_pybdsf", "reference_outputs"),
    1e-3,
    2
)

# Return exit status 0 if everything went fine, otherwise return 1.
if img and all(oper in img.completed_Ops for oper in operations) and outputs_agree:
    sys.exit(0)
else:
    sys.exit(1)
