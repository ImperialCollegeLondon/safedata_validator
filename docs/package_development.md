# Package development

## Installing the development version

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

All new package releases should be from the `master` branch, so the changes to `develop`
have to be moved here. This is achieved using a `release` branch.

```bash
git branch release/x.y.z
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

This gives all developers an opportunity to fix problems they have noticed, prior to
release. This can be done by making pull requests against the `release` branch. The
following checks should also be carried out:

* All GitHub continuous integration tests pass
* Upload to [Test PyPi](https://test.pypi.org) works correctly
* Online documentation builds properly for the `release` branch

Once the relevant changes have been made and checks have been performed the final commit
should bump the version using `poetry` so that the correct version is recorded in
`pyproject.toml`.

```bash
poetry version [major/minor/patch]
```

The `release` branch now is ready to be merged with the `master` branch. To do so you
should switch to the `master` branch and merge the `release` branch into it.

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
is uploaded to PyPi. You should perform a test upload when the `release` branch is
created. You should also perform another test upload whenever significant changes are
made to the `release` branch.

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
* Changes made to the `master` branch of the repository will automatically trigger a
  rebuild of the package documentation.
* To build the documentation for specific branches you need to login to [Read the
  Docs](https://readthedocs.org). You can then build whichever branch you require.

## TODO - Think of better name (again)

TODO - Write a section here that details final publication procedure

TODO -  Copy and paste this stuff into the test version where relevant

TODO - SOME OF THIS HAS BEEN MENTIONED IN AN ABOVE SECTION, CHANGE TO ADDRESS THAT
Before the package can be uploaded an API token for PyPi must be configured.

```bash
poetry config pypi-token.pypi my_api_token
```

This token should be a **personal** API token for PyPi, these can be generated through
[your PyPi account](https://pypi.org/account/login/). We use **personal** tokens rather
than a project specific token as the standard setup method with `poetry` only allows one
PyPi token to be saved. If necessary this issue can be circumvented by adding each
project as a new repository
[see](https://python-poetry.org/docs/repositories/#publishable-repositories) (with PyPi
remaining the repository published to) and then configuring this duplicate repository to
use the project specific token. We are not using this approach at present as we feel it
introduces unnecessary complexity. Your **personal** token will only allow you to
publish new package versions if you are a maintainer. If you wish to upload a new
package version you should therefore contact the current maintainers to request
maintainer status.

```bash
# Build source dist and binary
poetry build
# Upload just the new versions to the Test PyPi site and check it out
poetry publish -r testpypi
# Upload to the real PyPi site
poetry publish
```
