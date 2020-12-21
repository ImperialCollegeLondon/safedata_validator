# The safedata_validator package

This package provides methods to validate XLSX files containing formatted data and metadata using the SAFE Project data format.

See the main documentation for a detailed description and usage:

> [https://safedata_validator.readthedocs.io](https://safedata_validator.readthedocs.io)

The rest of this document describes the project development and building structure.

## Development notes


### Installing the development version

Testing the command line interface requires the package to be installed, but this is conveniently done in 'editable' mode, where it is always looking at the current state of the repo directory.

```bash
pip install -e .
```


### Testing

At present, the package is tested by running `safedata_validate` to check against a set of valid and bad input XLSX files. The code will move towards using the `pytest` framework, but at the moment is not structured to do this well. 

Although the package is registered with Travis CI, at the moment, there is no testing to meaningfully fail, 

### Publishing versions

The repository uses the `git flow` framework for releasing versions. When the `develop` branch contains a version of the code to be released as a new version, `git flow release start x.y.z` is used to start a release branch. When that has passed testing, the command `git flow release finish x.y.z` creates a new tagged version on the `master` branch. 

Note that the `master` branch should _only_ be used for new releases.

To publish a new version, once released, use the following steps.

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

 - Using `twine`:
 
 ```bash
 # Install twine if you don't have it
 pip install twine
 # Build source dist and binary
 python setup.py sdist bdist_wheel
 # Upload just the new versions to the Test Pypi site and check it out
 twine upload -r testpypi dist/safedata_validator-2.0.0*
 # Upload to the real PyPi site.
 twine upload -r pypi dist/safedata_validator-2.0.0*
 
 