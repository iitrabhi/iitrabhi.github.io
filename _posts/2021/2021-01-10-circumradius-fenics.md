---
layout: post
title: "How to get hmin for each cell in FEniCS"
categories: [coding]
tags: [fenics]
typora-root-url: ../../../website
---
FEniCS is full of useful commands that gives the user full control over the underlying mesh. One of the things that I have to constantly do in my research is to check the minimum element size in the mesh. This is very easy to check with FEniCS with `hmin()` command that **"computes the  maximum cell size in mesh, measured greatest distance between any two vertices of a cell"**. Thus, for a right angled triangle it will give us the length of the hypotenuses as the minimum element size in the mesh.

```python
from dolfin import *
%matplotlib inline
mesh = UnitSquareMesh(2,2)
plot(mesh)
mesh.hmin()
```

Output: 0.7071067811865476

![image-20210110161133917](/assets/images/image-20210110161133917.png)

But, sometimes I do require a vector containing the element size of all the elements inside the mesh. This is also easily possible with the `Circumradius(mesh)` command that **returns symbolic cell circumradius for given mesh**. Now using circumradius along with a discontinuous galerkin function space I can calculate the element size on each element of the cells with the following commands.

```python
c = Circumradius(mesh)
c_vector=project(c,FunctionSpace(mesh, "DG", 0)).vector()[:] * 2
```

This will give us a vector containig the  circumdiameter ( circumradius * 2) of each cell.

```
array([ 0.70710678,  0.70710678,  0.70710678,  0.70710678,  0.70710678, 0.70710678,  0.70710678,  0.70710678])
```

