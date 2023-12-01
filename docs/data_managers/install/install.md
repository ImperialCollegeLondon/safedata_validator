# Installing `safedata_validator`

The following steps should allow you to install `safedata_validator`:

* Download and install [Python](https://www.python.org/downloads/). On Windows, make
  sure to check the option to "Add `python.exe` to PATH" when running the installer, as
  this allows you to run Python from the command line.

* For the following commands, you will need to open a command line terminal. On Windows,
  that is the 'Command Prompt' application, which you can find from the search bar. On
  macOS and Linux machines, this will be the 'Terminal' application.

* Once you have a command line terminal, check that you can access the Python program
  using the command `python --version`. You should see something like this:

    ```term
    $ python --version
    Python 3.9.12
    ```

    If that did not work, you may still need to add `python` to your path - look for
    help online for your operating system.

* Install the `safedata_validator` package using the Python `pip` package installer.
  This will also install all of the other Python packages that `safedata_validator`
  needs to run.

    ```term
    pip install safedata_validator
    ```

* Once that has completed, you should now be able to run the following command, to
  verify that the command line tools have been installed correctly. See the
  [usage
  instructions](../command_line_tools/safedata_validate.md#using-safedatavalidate) for
  more information on `safedata_validate`.

    ```term
    $ safedata_validate --version
    safedata_validate 3.0.0a2
    ```

* You will now need to provide a [configuration](configuration.md) for the
  `safedata_validator` tools and install some required resources.
