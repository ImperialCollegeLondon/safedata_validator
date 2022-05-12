# Build notes 

TODO Incomplete

### Configure `git`

It is easier if `git` is configured to push new tags along with commits. This
essentially just means that new releases can be sent with a single commit, which
is simpler and saves Travis from building both the the code commit and then the
tagged version. This only needs to be set once.

```sh
set git config --global push.followTags true
```


Using git-flow and travis (although no tests at present so travis does very
little except for testing the thing builds)

Use git flow to create a release and then bump the version number in
`version.py`.

Check the package builds and installs locally:

```
python setup.py sdist bdist_wheel
```

Once all seems well,  finish the release, go to the master branch and push it to
create the tagged version on github.

Once that is done, switch back to `develop` and bump the version number to add
`.post9000` to show the code is in development again.

## PyPi

To upload the new version to testpypi, checkout master and run

```
python setup.py sdist bdist_wheel
```

Remembering to change the version number, you can then create an account at pypi
and testpypi and use `twine` to test:

```
twine upload -r testpypi dist/*1.2.8*
```

and then - once that seems to have gone ok - release the distribution for use
via `pip`

```
twine upload dist/*1.2.7*
```


## Documentation

The package documentation is maintained using [MkDocs](https://www.mkdocs.org/)
and deployed to
[https://safedata-validator.readthedocs.io/](https://safedata-validator.readthedocs.io/).
MkDocs can be installed using:

```
pip install mkdocs
```


In order to build and deploy the documentation.

* Edit the source files in the `docs` folder.
* Some of the documentation presents the command line help for the script tools
  in the package. To keep these synchronized with the codebase, the
  `docs/command_line_usage` directory contains a shell script that saves these
  outputs to file, so they can be included in the documentation. If the script
  commands are updated, these inputs need to be recreated.
* From the package root, run `mkdocs build`. This will create the docs site in
  the `site` folder - note that this folder is not included in the git repo.
* Once you have checked the local copy of the documentation, then commiting the
  changes to the repository will automatically trigger a rebuild of the package
  documentation.


