# Build notes 

TODO Incomplete

Using git-flow and travis (although no tests at present so travis does very little except for testing the thing builds)

Use git flow to create a release and then finish the release

Go to the master branch

```
python setup.py sdist bdist_wheel
```

Upload the new version to testpypi

```
twine upload -r testpypi dist/*1.2.7*
```

Once that seems to have gone ok,

```
twine upload dist/*1.2.7*
```


## Documentation

The package documentation is maintained using [MkDocs](https://www.mkdocs.org/) and deployed to [https://safedata-validator.readthedocs.io/](https://safedata-validator.readthedocs.io/). MkDocs can be installed using:

```
pip install mkdocs
```


In order to build and deploy the documentation.

* Edit the source files in the `docs` folder.
* From the package root, run `mkdocs build`. This will create the docs site in the `site` folder - note that this folder is not included in the git repo.
* Once you have checked the local copy of the documentation, then commiting the changes to the repository will automatically trigger a rebuild of the package documentation.



