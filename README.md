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

Although the package is registered with Travis CI, at the moment, there is no testing to meaningfully fail.

### Releasing a new version

The repository uses the `git flow` framework for releasing versions. When the `develop` branch contains a version of the code to be released as a new version, `git flow release start x.y.z` is used to start a release branch. When that has passed testing, the command `git flow release finish x.y.z` creates a new tagged version on the `master` branch. 

Note that the `master` branch should _only_ be used for new releases.

To publish a new version, first create the source distribution and a binary.

```{sh}
# Create distribution
python setup.py sdist bdist_wheel
```

If the distribution is to be tested locally, then the following command can be used to install the source distribution:

```bash
pip install safedata_validator --no-cache-dir --no-index --find-links dist/safedata_validator-x.y.z.tar.gz
```

The `--no-index` and `--find-links` options stop `pip` from using the web package index and point to the local distribution. The `--no-cache-dir` is vital as `pip` will cache the installed package (probably in `/tmp/`) and *unless the version number changes* will just keep using that cache. So, in development, this flag is needed to keep `pip` using the actual local distribution and **not** the outdated cached version.


### Publishing a new version

Versions are published to PyPi using `twine`.  As a first step, the new package is published to the PyPi **test** site, giving a chance to check the package upload before committing it to the live PyPi archive.
 
```bash
# Install twine if you don't have it
pip install twine
# Build source dist and binary
python setup.py sdist bdist_wheel
# Upload just the new versions to the Test Pypi site and check it out
twine upload -r testpypi dist/safedata_validator-x.y.z*
# Upload to the real PyPi site
twine upload -r pypi dist/safedata_validator-x.y.z*
```

