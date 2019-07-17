import setuptools

setuptools.setup(
    name="safe_dataset_checker",
    version="1.1.0",
    author="David Orme",
    author_email="d.orme@imperial.ac.uk",
    description="Code used to validate data files in the SAFE data submission format.",
    long_description="# TODO",
    long_description_content_type="text/markdown",
    url="https://github.com/imperial_college_london/safe_dataset_checker",
    #packages=setuptools.find_packages(),
    py_modules=['safe_dataset_checker'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
