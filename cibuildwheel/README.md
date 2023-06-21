# Using `cibuildwheel` 

The`cibuildwheel` tool is intended to be run on a CI server. However, for testing purposes, it is also possible to run `cibuildwheel` locally. In order to do so, you need to create a Python virtual environment and install `cibuildwheel` in it. Next go to the root directory of the project (i.e. the directory containing your `setup.py`, `setup.cfg`, and/or `pyproject.toml` file). Assuming you're on Linux, enter the following command:

```
cibuildwheel --platform linux
```

For more information check the `cibuildwheel --help` output, or goto https://cibuildwheel.readthedocs.io/en/stable/