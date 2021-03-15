---
layout: post
title: "How to start a new python project that you wish to support for long."
tag: 
  - misc
typora-root-url: ../../website
---

Most of the times I have to modify my initial code base so as to use it in future work. This kind of approach is not good as it requires significant amount change in the initial code base. Also, any kind of improvements added to the modified code are not reflected in the base code. All in all this approach of always changing the source code is highly inefficient.

Instead, I should have, from day one, tried to make a single repository and should have modularized the code. Now I am trying to do that.

This requires following a directory structure that I have seen many times on GitHub. For the sake of this project I will base my directory structure on popular libraries that are based on FEniCS. 

```bash
.
├── _data
├── _demo
├── _doc
├── _docker
├── _flats
├── _phase
├── _scripts
├── _test
└── setup.py
```



An overview of what each of these does:

| File/Directory | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `_data`        | This directory contains the data that will be used for test scripts. These will be mesh files created in gmsh |
| `_demo`        | This directory will contain the demos supplied with the package |
| `_doc`         | This directory will contain the documentation for the package. This will be automatically generated from python files |
| `_docker`      | This directory will contain the DockerFile to build the docker image for this project |
| `_flats`       | This directory will contain the flat programs based on which modules will be developed |
| `_phase`       | This is the main directory for the source code of the package |
| `_scripts`     | This directory will contain scripts to generate documentation |
| `_test`        | This directory will contain the unit tests                   |
| `setup.py`     | This file will help to install the package with `pip`        |

Once we have setup the repository in this fashion, I have to update the `setup.py` with the following code

```python
import codecs
import os

from setuptools import find_packages, setup

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))
about = {}
with open(os.path.join(base_dir, "pacakage", "__about__.py"), "rb") as f:
    exec(f.read(), about)


def read(fname):
    return codecs.open(os.path.join(base_dir, fname), encoding="utf-8").read()


setup(
    name="phase",
    version=about["__version__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    packages=find_packages(),
    description="...",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/iitrabhi/pacakage",
    license=about["__license__"],
    platforms="any",
    install_requires=["numpy", "matplotlib"],
    python_requires=">=3",
    # extras_require={
    #     "all": ["netCDF4", "h5py", "lxml"],
    #     "exodus": ["netCDF4"],
    #     "hdf5": ["h5py"],  # MED, MOAB, XDMF formats
    #     "xml": ["lxml"],  # Dolfin, VTU, XDMF, SVG
    # },
    classifiers=[
        about["__status__"],
        about["__license__"],
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    entry_points={"console_scripts": ["phase = phase._cli:main"]},
)
```

I don't know weather all of this is required but this works for now. This setup file will help us install the package with pip.

Now to develop this package we just need to go to the directory where we have `setup.py` and run the following command.

```bash
pip install -e .
```

This will install the package in editable mode so that we can debug it while developing. Once the package is install we can do the following 

``` python
from package import *	
```



