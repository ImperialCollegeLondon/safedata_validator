# The safedata_validator package

This package provides methods to validate XLSX files containing SAFE data formatted data and metadata.

See the main documentation:
  [https://safedata_validator.readthedocs.io](https://safedata_validator.readthedocs.io)




## Build notes

- local validation

Testing the command line interface requires the package to be installed. You can build the local distribution:

```{sh}
# Create distribution
python setup.py sdist bdist_wheel
```

That can then be installed using `pip`. The `--no-index` and `--find-links`
options stop `pip` from using the web package index and point to the local distribution.
The `--no-cache-dir` is vital as `pip` will cache the installed package
(probably in `/tmp/`) and *unless the version number changes* will just keep
using that cache. So, in development, this flag is needed to keep `pip` using the
actual local distribution and **not** the outdated cached version.

 `pip install safedata_validator --no-cache-dir --no-index --find-links dist/safedata_validator-2.0.0a0.tar.gz`

 - Using `twine`
