# Configuring the `safedata_validator` package

The `safedata_validator` package needs to be configured to use specific resources for
data validation using an [INI format](https://en.wikipedia.org/wiki/INI_file)
configuration file.

## Configuration file format

The file structure for the configuration file is a simple text file containing
the details below:

```ini
gbif_database = /path/to/local/backbone.sqlite3
ncbi_database = /path/to/local/ncbi_database.sqlite3
gazetteer = /path/to/gazeteer.geojson
location_aliases = /path/to/location_aliases.csv
use_project_ids = True

[extents]
temporal_soft_extent = 2002-02-02, 2030-01-31
temporal_hard_extent = 2002-02-01, 2030-02-01
latitudinal_hard_extent = -90, 90
latitudinal_soft_extent = -4, 6
longitudinal_hard_extent = -180, 180
longitudinal_soft_extent = 110, 120

[zenodo]
community_name = safe
contact_name = The SAFE Project
contact_affiliation = Imperial College London
contact_orcid = 0000-0003-3378-2814
use_sandbox = true
zenodo_api = https://api.zenodo.org
zenodo_token = abc
zenodo_sandbox_api = https://sandbox.zenodo.org
zenodo_sandbox_token = xyz

[metadata]
api = https://safeproject.net
token = xyz
ssl_verify = true
```

## Configuration file locations

You can put configuration files in any location and can even have multiple
configurations that point to different resources, such as different versions of
taxonomic databases. You can always use a specific configuration by providing a
`safedata_validator` tool with the path to that configuration file, using the
[`--resources` option](../command_line_tools/overview.md).

However, if you only use a single configuration, the `safedata_validator` tools can
automatically load that configuration if it is saved to a specific location. Conventions
for configuration file locations differ across operating systems and we use the
conventions used by the `appdirs` package.

The `safedata_validator` package will look for both user and site configuration files.
User configuration are only available for a particular user account, but site
configurations allow a data manager to set up a specific machine with data resources for
all users. If both are present, the user configuration is preferred.

For site and user configurations, the file **must** be called `safedata_validator.cfg`
and the user and site locations are:

<!-- markdownlint-disable MD046 -->

=== "Mac OS X"

    User
    : `/Users/username/Library/Application Support/safedata_validator/safedata_validator.cfg`

    Site
    : `/Library/Application Support/safedata_validator/safedata_validator.cfg`

=== "Windows"

    The repeated directory names are not an error!
    
    User
    : `C:\Users\username\AppData\Local\safedata_validator\safedata_validator\safedata_validator.cfg`

    Site
    : `C:\\ProgramData\\safedata_validator\\safedata_validator\safedata_validator.cfg`

=== "Linux"

    User
    : `/home/username/.config/safedata_validator/safedata_validator.cfg`

    Site
    : `/etc/xdg/safedata_validator/safedata_validator.cfg`

## Configuration components

The configuration file content breaks down into three distinct parts of the
`safedata_validator` workflow:

* **Validation**: the core resource files, such as taxonomic databases and gazetteer,
  and the extents settings for a project, that are used to check a dataset is valid.

* **Publication**: the account settings and access tokens used to publish validated
  datasets to Zenodo.

* **Metadata**: the URL and access tokens to send metadata about a dataset to a
  metadata server, allowing datasets to be accessed using the
  [`safedata`](https://imperialcollegelondon.github.io/safedata/) R package.

!!! Info

    You do not need to configure the publication and metadata sections if you are only
    using `safedata_validator` to validate datasets. You can also publish datasets to 
    Zenodo without needing to set up and configure a metadata server.

### Validation configuration

The validation configuration includes the following components, all of which must be
provided to validate datasets.

**The `gbif_database` element**

: This element provides the path to the [local GBIF backbone
  database](build_local_gbif.md) to be used in this configuration.

**The `ncbi_database` element**

: This element provides the path to the [local NCBI backbone
  database](build_local_ncbi.md) to be used in this configuration.

**The `gazetteer` and `location_aliases` elements**

: These elements provide the paths to the [location database](gazetteer.md) for the
  project.

**The `use_project_ids` element**

: The `safedata` system was initially developed as part of the SAFE project, which
  organises datasets into projects. Each project is identified using a unique
  identification number, and each dataset is associated with one or more projects using
  these project ID numbers. We anticipate that other deployments of the `safedata`
  system may not wish to organise datasets in this manner. If so, the checking of
  project IDs can be turned off by setting the `use_project_ids` option  to `False`.
  With this option off, datasets will **fail validation** if they contain project IDs.

    !!! Warning
        Each deployment of the `safedata` system will have to make a choice of whether
        or not to organise datasets into project. Once this choice has been made it will
        be binding. As such, this option should only be changed by the relevant data 
        manager during the initial setup of the data system.

**The `extents` element**

: The `safedata_validator` package tracks the geographic and temporal extents of
    datasets, which are needed to generate Gemini 2 metadata for a dataset. A project
    can provide both **soft extents**, which cause validation to raise a warning, and
    **hard extents**, which case validation to fail.

    By default, only hard extents on geographic coordinates are applied:

    * `latitudinal_hard_extent`: (-90, 90)
    * `longitudinal_hard_extent`: (-180, 180)

### Publication configuration

The `safedata_zenodo` command line tool provides functionality to for publish validated
datasets to the Zenodo data repository. You will need to:

* Create an [account with Zenodo](https://zenodo.org/signup/).
* Use this account to [create a new Zenodo community](https://zenodo.org/communities)
  that will be used to group all published datasets.
* From your user account, [generate an access
  token](https://zenodo.org/account/settings/applications/) that will allow
  `safedata_zenodo` to authenticate access to the Zenodo API for uploading datasets. The
  token will need to have the `deposit:actions` scope.

It is **extremely highly recommended** that you repeat the sign up steps above using the
[sandbox version of Zenodo](https://sandbox.zenodo.org),
using the same community name. The sandbox site:

* provides an identical environment to the real Zenodo site, so that upload
  function outputs can be checked, without adding unnecessary (or invalid) datasets to
  the official repository.

* allows you to test the `safedata_validator` workflow all the way through to dataset
  publication without generating live DOIs.
  
The configuration element `use_sandbox` can be used to switch between testing and actual
publication. For more information on the sandbox site, see [this
page](https://developers.zenodo.org/#testing).

!!! Important

    The access tokens and credentials generated above **provide administrator** access
    to your Zenodo community and datasets and you should store them securely. As they 
    are included in the `safedata_validator` configuration file, you must be careful 
    about who has access to a computer setup to provide validation.

Once you have been through this process, you can then fill out the following
configuration elements:

**The `community_name` element**

: This sets the Zenodo community to be used for publishing datasets on both the main and
  sandbox Zenodo sites.

**The `contact_name`, `contact_affiliation`, and `contact_orcid` elements**

: All datasets will be published using these contact details and should provide a
  permanent set of contact information for the project datasets. Note that this is
  different from the **dataset authors**.

**The `zenodo_token` and `zenodo_sandbox_token` elements**

: These are the personal access tokens generated for the user accounts on the main and
  sandbox Zenodo sites.

**The `use_sandbox` element**

: When this element is `true`, all datasets will be published to the testing sandbox
  site. Set this to `false` when you are ready to actually start publishing datasets.

### Metadata configuration

Zenodo only allows a fairly limited amount of metadata to be stored for each dataset.
While this is completely adequate to describe the contents of a dataset, more extensive
metadata must be stored elsewhere if detailed searches within datasets are desired.

The [`safedata_server`](https://github.com/ImperialCollegeLondon/safedata_server) web
application allows more detailed metadata to be made available to end users of the data
and provides an API to aid data discovery and downloading. This API is used extensively
by the [`safedata`](https://imperialcollegelondon.github.io/safedata/) R package.

To use this system, you will need to deploy the web application to a publically
accessible URL and then configure the following elements:

**The `api` and `token` elements**

: These provide the URL of the metadata server API and an access token required to
  authenticate data upload to that server.

**The `ssl_verify` element**

: In production use, the metadata server should be set up with a properly validated SSL
  certificate to allow HTTPS, and this is relatively easy using LetsEncrypt. However,
  when setting up and testing a system, requiring a valid certificate can be a road
  block and this element allows SSL verification to be turned off.
