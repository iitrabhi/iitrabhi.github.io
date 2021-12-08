---
layout: post
title: "How to build petsc with complex support."
tag: 
  - fenics
typora-root-url: ../../website
---

## The problem 

Need complex number support in PETSc. The standard version of FEniCS is compiled with real version of petsc. For the eigenvalue problem we need complex numbers. I have gone to the petsc4py repo and found some useful information. 

## The resolution

First of all we need to check the installation of petsc. 

```bash
from petsc4py import PETSc
print(PETSc.ScalarType)

output >> <class 'numpy.float64'>
<class 'numpy.complex128'>
```

Now as per [Drew Parsons](https://bitbucket.org/petsc/petsc4py/issues/115/is-concurrent-installation-of-real-and), *the structure of the PETSc installation makes it possible to install different versions and different configurations side by side. You can choose which one you want to work with at run time (or build time) using PETSC_DIR, choosing between real and complex number support, for instance.*

Also as per the official repo of [petsc4py](https://gitlab.com/petsc/petsc4py) - *"The petsc4py project now lives in the [main PETSc repository](https://gitlab.com/petsc/petsc) and can be installed by configuring `--with-petsc4py` or via the [Python package](https://pypi.org/project/petsc4py/)."* 

For this we first need to install PETSc with different `scalar-types` and mark those installation with `PETSC_ARCH` environment variable. Then we can switch between the different installations.

```bash
# Real, 32-bit int
    python3 ./configure \
    PETSC_ARCH=linux-gnu-real-32 \
    --with-scalar-type=real && \
    make PETSC_DIR=/usr/local/petsc PETSC_ARCH=linux-gnu-real-32 ${MAKEFLAGS} all && \
# Complex, 32-bit int
    python3 ./configure \
    PETSC_ARCH=linux-gnu-complex-32 \
    --with-scalar-type=complex && \
    make PETSC_DIR=/usr/local/petsc PETSC_ARCH=linux-gnu-complex-32 ${MAKEFLAGS} all && \
# Real, 64-bit int
    python3 ./configure \
    PETSC_ARCH=linux-gnu-real-64 \
    --with-scalar-type=real && \
    make PETSC_DIR=/usr/local/petsc PETSC_ARCH=linux-gnu-real-64 ${MAKEFLAGS} all && \
# Complex, 64-bit int
    python3 ./configure \
    PETSC_ARCH=linux-gnu-complex-64 \
    --with-scalar-type=complex && \
    make PETSC_DIR=/usr/local/petsc PETSC_ARCH=linux-gnu-complex-64 ${MAKEFLAGS} all && \
 
```

Once PETSc is installed with all of the above `ARCH` we can install `petsc4py` and tell it about the different `ARCH`. The source for `petsc4py` is now in the `PETSc` directory

```bash
# Install petsc4py
    cd src/binding/petsc4py && \
    PETSC_ARCH=linux-gnu-real-32:linux-gnu-complex-32:linux-gnu-real-64:linux-gnu-complex-64 pip3 install --no-cache-dir .
```

Now by changing the environment variable we can activate the specific version of petsc and use it. For example to use the complex build we just need to set the environment variable as such

```bash
export PETSC_ARCH=linux-gnu-complex-32
```

After this by running `Scalartype` we will get complex support

```bash
from petsc4py import PETSc
print(PETSc.ScalarType)

output >> <class 'numpy.complex128'>
```





## References

- https://bitbucket.org/petsc/petsc4py/issues/115/is-concurrent-installation-of-real-and
- https://bitbucket.org/petsc/petsc4py/issues/83/how-to-create-complex-matrix-by-using
- https://gitlab.com/petsc/petsc
- [Major understanding of installation](https://github.com/FEniCS/dolfinx/blob/main/docker)

