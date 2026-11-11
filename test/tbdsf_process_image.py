import argparse
import bdsf
import sys
import filecmp
from compare_outputs import compare_results

# The outputs are compared with the reference outputs, within a relative
# tolerance. The default is suitable when the test runs on the same kind of
# system (architecture, libraries) as the one that generated the reference
# outputs (x86-64 Linux). Some fits in this test are ill-conditioned, so that
# tiny differences in rounding can change their results; therefore, CI jobs,
# which may run on other systems, use a larger tolerance (--rtol 0.03).
parser = argparse.ArgumentParser(description="Test bdsf.process_image()")
parser.add_argument(
    "--rtol",
    type=float,
    default=0.001,
    help="Relative tolerance for comparing outputs with the reference "
    "outputs (default=%(default)s)",
)
args = parser.parse_args()

# Process the image
img = bdsf.process_image('tbdsf_process_image.in', ncores=2, output_all=True)

# List of operations that must have been done on `img`.
operations = [
    'readimage', 'collapse', 'preprocess', 'rmsimage', 'threshold',
    'islands', 'gausfit', 'gaul2srl', 'make_residimage', 'wavelet_atrous',
    'shapelets', 'spectralindex', 'polarisation', 'psf_vary', 'cleanup'
]

# Check if the outputs agree with the reference ones. Catalogs and other
# text files are not compared, because small differences in a few fits can
# change them considerably.
outputs_agree = compare_results(
    filecmp.dircmp("tbdsf_process_image.in_fits_pybdsf", "reference_outputs"),
    rtol=args.rtol,
    verbosity=0,  # only print warnings and errors
    check_images=True,  # compare the rms values (within rtol) of images
    check_catalogs=False,  # compare entries (within rtol) in FITS catalogs
    check_other=False  # compare entries in text files
)

# Return exit status 0 if everything went fine, otherwise return 1.
if img and all(oper in img.completed_Ops for oper in operations) and outputs_agree:
    sys.exit(0)
else:
    sys.exit(1)
