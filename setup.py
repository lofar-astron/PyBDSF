from skbuild import setup

# Scikit-build requires some options to be passed in the call to `setup()`.
# Hence, these cannot be set in `pyproject.toml`
setup(
    packages=["bdsf", "bdsf.nat"],
    include_package_data=False,
)
