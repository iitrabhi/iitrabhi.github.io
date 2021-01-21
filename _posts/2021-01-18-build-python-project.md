---
layout: post
title: "How to start a new python project that you wish to support for long."
categories: misc
typora-root-url: ../../website
---

Most of the times I have to modify my initial code base so as to use it in future work. This kind of approach is not good as it requires significant amount change in the initial code base. Also, any kind of improvements added to the modified code are not reflected in the base code. All in all this approach of always changing the source code is highly inefficient.

Instead, I should have, from day one, tried to make a single repository and should have modularized the code. Now I am trying to do that.

This requires following a directory structure that I have seen many times on GitHub. For the sake of this project I will base my directory structure on popular libraries that are based on FEniCS. 

```
 .
 +-- _data
 +-- _demo
 +-- _doc
 +-- _docker
 +-- _flats
 +-- _phase
 +-- _scripts
 +-- _test
 +-- setup.py
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



