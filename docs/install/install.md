# Installing `safedata_validator`

The following steps should allow you to install `safedata_validator`:

1. If you're using a Mac or Linux, then it is almost certain that you already
   have Python installed. You may also do if you are using Windows. To find out,
   open a command line window (run Terminal on a Mac, or `cmd` on Windows) and
   then type `python` (on Mac or Unix) or `python.exe` on Windows).

2. If some text and a prompt (`>>>`) appears then you have Python. Check the
   version number in the first line: if it doesn't start Python 2.7 then you
   currently need to install Python 2.7 to run the checker. You can have
   multiple versions of Python installed, but it is going to be more complicated
   than is covered here.

    If you've got Python 2.7 then type `quit()` and skip to step 4.
    
    If you get a line that says the command is not found then you need to
    install Python. Download a copy from here:

    [https://www.python.org/downloads/](https://www.python.org/downloads/)

    The code is currently written to use Python 2.7, so make sure you download
    an installer for the most recent version of Python 2.7 and not the more
    recent Python 3 versions. For Windows, choose one of the MSI installer
    options.

3. Repeat the command line check from the first step: if this still doesn't work
   then you probably just need to tell the computer where to find Python: search
   online for instructions to add python to the `PATH` environment variable. On
   Windows, you will want to add (using the typical install location)
   `C:\Python27` and `C:\Python27\scripts`.

4. The `safedata_validator` package requires some extra Python packages that may
   need to be installed. You can install of the required packages and
   `safedata_validator` itself using the `pip` package installer. At the command
   line, type:

        pip install --user openpyxl requests simplejson shapely appdirs python-dateutil safedata_validator

    The additional packages allow Python to: read Excel files, get validation
    data over the internet, handle JSON formatted data, validate WKT formatted
    GIS vector data, handle configuration file locations and parse dates more
    easily. If you want to install `safedata_validator` for all users on a
    computer, you will need to remove `--user`.

5. The package will have created  a system command `safedata_validate` that you
   can use to validate a file from the command line. Open a command line
   terminal,  and run the following:

        safedata_validate -h

    You should now see the [usage instructions](../usage/usage.md).

 6. In addition to the command line option, you should now be able to `import
    safedata_validator` within Python, which will allow you to use the  Dataset
    class and methods defined  in the package within your own code. 
