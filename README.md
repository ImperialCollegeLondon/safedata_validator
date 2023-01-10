# The safedata_validator package

This package provides methods to validate XLSX files containing formatted data and
metadata using the SAFE Project data format.

See the main documentation for a detailed description and usage:

> [https://safedata_validator.readthedocs.io](https://safedata_validator.readthedocs.io)

The rest of this document describes the project development and building structure.

## Development notes

### Installing the development version

This package makes use of the `python` dependency manager
[`poetry`](https://python-poetry.org). This means that the package and all required
dependencies can be installed with a single command.

```bash
poetry install
```

This installed package is 'editable', i.e. it changes with the current state of the repo
directory. The package is installed as part of a virtual environment, so can only be
used when the relevant environment is active. This environment can be activated with a
single `poetry` command.

```bash
poetry shell
```

### Testing

Testing for this package makes use of the `pytest` framework. When new functions are
added to this package unit tests for them must also be added. These operate as a check
that functions still operate in the manner they were originally designed to. Either all
unit tests can be run locally.

```bash
pytest
```

Or a specific testing file can be run.

```bash
pytest test/test_specific_module.py
```

These unit tests are also run as part of our continuous integration workflow, which runs
whenever commits are made to this repository and for all pull requests.

### Releasing a new version

All new package releases should be from the `master` branch, so the changes to `develop`
have to be moved here. This is achieved using a `release` branch.

```bash
git branch release/x.y.z
```

This branch should be pushed to the remote repo to ensure that other developers have
access to it.

```bash
git push --set-upstream origin release/x.y.z
```

This gives all developers an opportunity to fix problems they have noticed, prior to
release. This can be done by making pull requests against the `release` branch. Once
these changes have been made, documentation has been checked, and all tests pass, the
`release` branch is ready to be merged with the `master` branch.

```bash
git switch master
git merge release/x.y.z
```

A tag should be added marking the package version, and then the updated `master` branch
should be pushed to the remote repository.

```bash
git tag x.y.z
git push origin x.y.z
```

Finally, any changes added to the `release` branch should be merged back into `develop`.

```bash
git switch develop
git merge release/x.y.z
git push
```

Note that the `master` branch should _only_ be used for new releases.

To publish a new version, first create the source distribution and a binary.

```{sh}
# Create distribution
python setup.py sdist bdist_wheel
```

If the distribution is to be tested locally, then the following command can be used to
install the source distribution:

```bash
pip install safedata_validator --no-cache-dir --no-index --find-links dist/safedata_validator-x.y.z.tar.gz
```

The `--no-index` and `--find-links` options stop `pip` from using the web package index
and point to the local distribution. The `--no-cache-dir` is vital as `pip` will cache
the installed package (probably in `/tmp/`) and _unless the version number changes_ will
just keep using that cache. So, in development, this flag is needed to keep `pip` using
the actual local distribution and **not** the outdated cached version.

### Publishing a new version

Versions are published to PyPi using `twine`.  As a first step, the new package is
published to the PyPi **test** site, giving a chance to check the package upload before
committing it to the live PyPi archive.

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
