# Installing `safedata_validator`

The following steps should allow you to install `safedata_validator`:

1. If you're using a Mac or Linux, then it is almost certain that you already
   have Python installed. You may also do if you are using Windows. To find out,
   open a command line window (run Terminal on a Mac, or `cmd` on Windows) and
   then type `python` (on Mac or Unix) or `python.exe` on Windows).

2. If some text and a prompt (`>>>`) appears then you have Python. Check the
   version number in the first line: the package requires at least Python 3.9.

   If you get a line that says the command is not found, then it could be one of two
   things:

   * You need to install Python. Download a copy from here:

     [https://www.python.org/downloads/](https://www.python.org/downloads/)

     If you do see a Python prompt, but the version is less than 3.9, then you will
     need to update it.

   * Your computer has not been set up to find Python: search
     online for instructions to add python to the `PATH` environment variable. On
     Windows, you will want to add (using the typical install location)
     `C:\Python39` and `C:\Python39\scripts`.

3. You can now install `safedata_validator` and the other packages it depends on using
   the `pip` package installer. At the command
   line, type:

        pip install --user safedata_validator

   If you want to install `safedata_validator` for all users on a computer, you will
   need to remove `--user`.

4. The package will have created  a system command `safedata_validate` that you
   can use to validate a file from the command line. Open a command line
   terminal,  and run the following:

        safedata_validate -h

    You should now see the [usage instructions](../command_line_tools/safedata_validate.md#using-safedatavalidate).

5. In addition to the command line option, you should now be able to `import
   safedata_validator` within Python, which will allow you to use the
   class and methods defined in the package within your own code.
