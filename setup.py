from skbuild import setup  # This line replaces 'from setuptools import setup'
# from setuptools_scm import get_version

setup(
    # name='bdsf',
    # version=get_version(),
    author='David Rafferty',
    author_email='drafferty@hs.uni-hamburg.de',
    url='https://github.com/lofar-astron/PyBDSF',
    description='Blob Detection and Source Finder',
    long_description=open('README.rst', 'rt').read(),
    long_description_content_type='text/x-rst',
    platforms='Linux, Mac OS X',
    packages=['bdsf', 'bdsf.nat'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: C++',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Astronomy'
    ],
    extras_require={
        'ishell': [
            'ipython!=8.11.*,!=8.12.*,!=8.13.*,!=8.14.*,!=8.15.*',
            'matplotlib',
        ],
    },
    install_requires=['backports.shutil_get_terminal_size',
                        'astropy', 'numpy', 'scipy'],
    python_requires=">=3.7",
    entry_points = {
        'console_scripts': [
            'pybdsf = bdsf.pybdsf:main [ishell]',
            'pybdsm = bdsf.pybdsf:main [ishell]'
        ]
    },
    zip_safe=False,
    include_package_data=False,
)
