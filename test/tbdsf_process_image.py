import bdsf
import sys

# Process the image
img = bdsf.process_image('tbdsf_process_image.in', ncores=2)

# List of operations that must have been done on `img`.
operations = [
    'readimage', 'collapse', 'preprocess', 'rmsimage', 'threshold',
    'islands', 'gausfit', 'gaul2srl', 'make_residimage', 'wavelet_atrous',
    'shapelets', 'spectralindex', 'polarisation', 'psf_vary', 'cleanup'
]

# Return exit status 0 if everything went fine, otherwise return 1.
if img and all(oper in img.completed_Ops for oper in operations):
    sys.exit(0)
else:
    sys.exit(1)
