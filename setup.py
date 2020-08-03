import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = {}
with open("safedata_validator/version.py") as fp:
    exec(fp.read(), version)

setuptools.setup(
    name="safedata_validator",
    version=version['__version__'],
    author="David Orme",
    author_email="d.orme@imperial.ac.uk",
    description="Validator for data files in the SAFE data submission format.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/imperial_college_london/safedata_validator",
    packages=['safedata_validator'],
    entry_points = {
            'console_scripts':
             ['safedata_validate=safedata_validator:_safedata_validator_cli']
    },
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
