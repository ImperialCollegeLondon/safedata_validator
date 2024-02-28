# Package development

## Installing the development version

First, clone the GitHub repository to download the latest development code:

```bash
git clone https://github.com/ImperialCollegeLondon/safedata_validator.git
```

The package makes use of the `python` dependency manager
[`poetry`](https://python-poetry.org), so you will also need to install poetry:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

Once those two steps are complete, you can run a single command from the package root to
install the package and all dependencies:

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

## Testing

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

## Releasing a new version

All new package releases should be from the `main` branch, so the changes to `develop`
have to be moved here. This is achieved using a `release` branch.

```bash
git branch release/x.y.z
git switch release/x.y.z
```

The version of this branch should be updated to a pre-release version using `poetry`.
There are multiple options here: `premajor` should be used for major versions (e.g.
`2.0.0`), `preminor` for minor (e.g. `2.1.0`), and prepatch for patch versions (e.g.
`2.1.1`).

```bash
poetry version [premajor/preminor/prepatch]
```

The change this causes to `pyproject.toml` should be committed, and the branch is now
ready to be pushed to the remote repo to ensure that other developers have access to it.

```bash
git push --set-upstream origin release/x.y.z
```

This branch ultimately needs to be merged into `main` (and back into `develop`), so once
the release branch is on GitHub, a **pull request should be made against `main`** from
the release branch. This will cause the GitHub continuous integration tests and
documentation building to run, validating the release branch. Other developers can also
look at the PR to review the changes prior to release, and to make commits to the branch
to update the PR or fix issues.

Once everything seems to be running smoothly, the next steps are to make sure that the
publication process and website build work correctly. The package can be published to
the Test PyPi site (see below for details) using:

```bash
poetry build
poetry publish -r test-pypi
```

ReadTheDocs then needs to be updated to build the `release/x.y.z` branch. The branch
needs to be 'Activated' from the Versions tab on the RTD project admin site - it should
be 'Active' but also 'Hidden'.

Once the relevant changes have been made and checks have been performed the final commit
should bump the version using `poetry` so that the correct version is recorded in
`pyproject.toml`.

```bash
poetry version [major/minor/patch]
```

Alternatively, edit `pyproject.toml` by hand if the release tag to be used is unusual in
some way (e.g. `x.y.zrc1`). Commit this last change!

The `release` branch now is ready to be merged with the `main` branch and the pull
request should be accepted and merged online. A tag should be added marking the package
version, and then the updated `main` branch should be pushed to the remote repository.

```bash
git tag x.y.z
git push origin x.y.z
```

A PR should also be created online to merge any changes added to the `release` branch
back into `develop`.

Note that the `main` branch should _only_ be used for new releases.

## Uploading package releases to test PyPi

It can often be the case that a package that appears to build fine locally has errors
that prevent it from uploading properly to PyPi. By first uploading to [test PyPi
site](https://test.pypi.org/) these kind of errors can be caught without clogging up the
real PyPi site with broken packages. Upload of new package versions occurs via `poetry`.
This means that `poetry` must be configured to have access to the test PyPi. To do this,
test PyPi must be added as a repository, and a valid API access token must be associated
with the repository.

```bash
poetry config repositories.testpypi https://test.pypi.org/legacy/
poetry config pypi-token.testpypi my_test_api_token
```

This token should be a **personal** API token for PyPi, these can be generated through
[your test PyPi account](https://test.pypi.org/account/login/). We are now setup to
publish to test PyPi, but before the package is published it must first be built. This
also done using `poetry`.

```bash
poetry build
```

This produces both a `sdist` source distribution, and a `wheel` compiled package. These
can then be published to test PyPi.

```bash
poetry publish -r testpypi
```

It is important to note that a test upload should **always** be done before the package
is uploaded to PyPi. You should perform a test upload when the `release` branch
otherwise ready to merge to `main`. If the `release` branch changes after this point,
use `poetry version prerelease` to increment the version and run another test upload to
confirm that the final version can be uploaded cleanly.

## Documentation

The package documentation is maintained using [MkDocs](https://www.mkdocs.org/)
and deployed to
[https://safedata-validator.readthedocs.io/](https://safedata-validator.readthedocs.io/).
MkDocs is installed automatically to the `poetry` virtual environment.

In order to build and deploy the documentation.

* Edit the source files in the `docs` folder.
* Some of the documentation presents the command line help for the script tools
  in the package. To keep these synchronized with the codebase, the
  `docs/command_line_usage` directory contains a shell script that saves these
  outputs to file, so they can be included in the documentation. If the script
  commands are updated, these inputs need to be recreated.
* From the package root, run `mkdocs build`. This will create the docs site in
  the `site` folder - note that this folder is not included in the git repo.
* Changes made to the `main` branch of the repository will automatically trigger a
  rebuild of the package documentation.
* To build the documentation for specific branches you need to login to [Read the
  Docs](https://readthedocs.org). You can then build whichever branch you require.

## Final publication of new package version

Once a `release` branch has passed all the tests and been merged into `main`, it
should be published to [PyPi](https://pypi.org/). This allows users to install the new
package version via `pip`. As with test PyPi, publication is handled by the `poetry`
package manager. PyPi is automatically configured as the default upload repository, so
in this case you only need to add an API token to the `poetry` configuration for `PyPi`.

```bash
poetry config pypi-token.pypi my_api_token
```

As with the token for test PyPi, this token should be a **personal** API token, these
can be generated through [your PyPi account](https://pypi.org/account/login/). We use
**personal** tokens rather than a project specific token as the standard setup method
with `poetry` only allows one PyPi token to be saved. If necessary this issue can be
circumvented by adding each project as a new repository
[see](https://python-poetry.org/docs/repositories/#publishable-repositories) (with PyPi
remaining the repository published to) and then configuring this duplicate repository to
use the project specific token. We are not using this approach at present as we feel it
introduces unnecessary complexity.

It should be noted that your **personal** token will only allow you to publish new
package versions if you are a maintainer. If you wish to upload a new package version
you should therefore contact the current maintainers to request maintainer status.

Once `poetry` has been setup to allow publication of `safedata_validator` to PyPi, the
new package version can then be published from the `main` branch:

```bash
git switch main
poetry build
poetry publish
git switch develop
```
