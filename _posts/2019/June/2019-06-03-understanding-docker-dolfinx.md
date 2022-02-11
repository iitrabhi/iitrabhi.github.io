---
layout: post
title: "Understanding the dockerfile of DOLFIN-X"
description: "Handling a complex project like FEniCS required Docker."
tag: 
  - gsoc
typora-root-url: ../../../../website
---

There are excellent articles available on the web explaining the need and use of Docker. The interested reader is requested to start with the [official documentation of Docker](https://docs.docker.com/engine/docker-overview/). In this article I will try to present my **limited** personal understanding of Docker and why is it necessary for the project FEniCS. As per the official wiki

> Docker is a set of coupled software-as-a-service and platform-as-a-service products that use operating-system-level virtualization to develop and deliver software in packages called containers. The software that hosts the containers is called Docker Engine.

This means that docker engine is kind of a virtual PC onto which we can load an OS of our choice (preferably UBUNTU) and then install further application packages into it. Immediately two questions come to mind:

1. What is the difference between docker engine and something more conventional like a virtual machine?
2. Why do we even need something like docker?

For the first question, the answer is quite simple, docker engine and docker images are lightweight as compared to their counterpart. Just for an example, the official ubuntu “.iso” file is well over 1GB in size whereas the official [docker image for Ubuntu is less than 100MB in size](https://hub.docker.com/_/ubuntu?tab=tags).

The answer to why do we need docker is a bit more complicated. If you have ever developed any scientific package in any programming language you would have used a few of the many libraries available on the internet. Most of the scientific programs that are build using python use at least numpy, scipy or pandas. When the end user of your software package installs your software, it is required that s/he also installs all the dependencies. This could be tricky and as an individual developer or a small organization, you have no guarantee that your library and its dependencies would work on all the available varity of operating systems and hardware available in the market. This is especially true with compiled languages such as C++ which require a platform specific compiler to work.

Docker helps in this situation by providing operating-system-level virtualization to develop and deliver software in packages called containers. Development with docker could be broken down into the following steps:

- Install docker engine in your system.
- Load the required OS (UBUNTU) into docker engine. This is done by using a docker image, which you can think of as an alternative to the “.iso” image file.
- Develop your application inside docker over the virtual OS.

This process, in essence, helps you to forget about OS and hardware specific programming and helps you to focus on the core development of your library.

We can build our own docker image with the help of Dockerfile. More on it [here ](https://docs.docker.com/engine/reference/builder/)and [here](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/). The official repo of DOLFIN-X provides us with a Dockerfile. Dockerfile has the instructions to install all the required dependencies inside an OS and then create a docker image from it. I will try to explain all the commands in the [official Dockerfile of DOLFIN-X.](https://github.com/FEniCS/dolfinx/blob/master/Dockerfile)

------

**ARG** variables are used to specify variables that are available to the RUN instruction. Environment variables defined using the ENV instruction always override an ARG instruction of the same name. These are used to specify the version number of packages that were used at the time of development. The docker image is built with the following version of packages.

```
ARG GMSH_VERSION=4.2.2
ARG PYBIND11_VERSION=2.2.4
ARG SPDLOG_VERSION=1.3.1
ARG PETSC_VERSION=3.10.4
ARG SLEPC_VERSION=3.10.2
ARG PETSC4PY_VERSION=3.10.1
ARG SLEPC4PY_VERSION=3.10.0
ARG TINI_VERSION=v0.18.0

ARG MAKEFLAGS
ARG PETSC_SLEPC_OPTFLAGS="-02 -g"
ARG PETSC_SLEPC_DEBUGGING="yes"
```

The **FROM** instruction initializes a new build stage and sets the **Base Image** for subsequent instructions. You can think of this as specifying the operating system in which you wish to install further applications. FROM command pulls the image specified from docker-hub and then sets it as the base image i.e. the OS for our installs. You have to look for the exact name of the image on docker-hub that you wish to use.

You can add labels to your image to help organize images by project, record licensing information, to aid in automation, or for other reasons. For each label, add a line beginning with LABEL and with one or more key-value pairs.

```
FROM ubuntu:18.04 as base
LABEL maintainer="fenics-project <fenics-support@googlegroups.org>"
LABEL description="Base image for real and complex FEniCS test environments"
```

An **ARG** declared before a **FROM** is outside of a build stage, so it can’t be used in any instruction after a FROM. To use the default value of an ARG declared before the first FROM use an ARG instruction without a value inside of a build stage:

```
ARG GMSH_VERSION
ARG PYBIND11_VERSION
ARG SPDLOG_VERSION
```

The **WORKDIR** instruction sets the working directory for any RUN, CMD, ENTRYPOINT, COPY and ADD instructions that follow it in the Dockerfile. If the WORKDIR doesn’t exist, it will be created even if it’s not used in any subsequent Dockerfile instruction.

The WORKDIR instruction can be used multiple times in a Dockerfile. If a relative path is provided, it will be relative to the path of the previous WORKDIR instruction.

This directory is created in the root of the Ubuntu image. To visualize the whole process think that you are logged into ubuntu on your desktop pc and you wish to create a directory named tmp in the root directory.

```
WORKDIR /tmp
```

The **ENV** instruction sets the environment variable to the value. This value will be in the environment for all subsequent instructions in the build stage and can be replaced inline in many as well.

```
ENV OPENBLAS_NUM_THREADS=1 \
    OPENBLAS_VERBOSE=0
```

The RUN instruction will execute any commands in a new layer on top of the current image and commit the results. The resulting committed image will be used for the next step in the Dockerfile.

This command is used to run the commands that we usually run in a shell (i.e. terminal(MAC) or cmd(WINDOWS)). We just need to append the command RUN in front of the instructions that we wish to execute in the shell.

```
# Install dependencies available via apt-get.
# First set of packages are required to build and run FEniCS.
# Second set of packages are recommended and/or required to build documentation or tests.
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get -qq update && \
    apt-get -yq --with-new-pkgs -o Dpkg::Options::="--force-confold" upgrade && \
    apt-get -y install \
        cmake \
        g++ \
        gfortran \
        libboost-dev \
        libboost-filesystem-dev \
        libboost-iostreams-dev \
        libboost-math-dev \
        libboost-program-options-dev \
        libboost-system-dev \
        libboost-thread-dev \
        libboost-timer-dev \
        libeigen3-dev \
        libhdf5-mpich-dev \
        liblapack-dev \
        libmpich-dev \
        libopenblas-dev \
        mpich \
        ninja-build \
        pkg-config \
        python3-dev \
        python3-matplotlib \
        python3-numpy \
        python3-pip \
        python3-scipy \
        python3-setuptools && \
    apt-get -y install \
        doxygen \
        git \
        graphviz \
        valgrind \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Download Install gmsh############################################################
# Here we are entering the directory /usr/local in ubuntu and then downloading and unziping gmsh into it.
RUN cd /usr/local && \
    wget -nc --quiet http://gmsh.info/bin/Linux/gmsh-${GMSH_VERSION}-Linux64.tgz && \
    tar -xf gmsh-${GMSH_VERSION}-Linux64.tgz
# Adding the folder to PATH variable ensures that the system would recoganize a call to gmsh
ENV PATH=/usr/local/gmsh-${GMSH_VERSION}-Linux64/bin:$PATH
####################################################################################
# Install Python packages (via pip)
# First set of packages are required to build and run FEniCS.
# Second set of packages are recommended and/or required to build documentation or run tests.
RUN pip3 install --no-cache-dir mpi4py numba && \
    pip3 install --no-cache-dir cffi decorator flake8 pygmsh pytest pytest-xdist sphinx sphinx_rtd_theme

# Install pybind11
RUN wget -nc --quiet https://github.com/pybind/pybind11/archive/v${PYBIND11_VERSION}.tar.gz && \
    tar -xf v${PYBIND11_VERSION}.tar.gz && \
    cd pybind11-${PYBIND11_VERSION} && \
    mkdir build && \
    cd build && \
    cmake -DPYBIND11_TEST=False ../ && \
    make install && \
    rm -rf /tmp/*
```

Here the process of installation of the OS and all the required dependencies for the base image ends and we return back to the root directory.

```
WORKDIR /root
```

Next, there are commands to create different flavors of the installation of FEniCS based on the “**base**” image. All the commands follow a similar pattern.
