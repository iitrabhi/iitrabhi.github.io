---
layout: post
title: "FEniCS: Completion of phase one"
description: "Successful completion of phase one of GSoC2019"
categories: [gsoc,coding]
typora-root-url: ../../../../website
---

The coding period officially started on May 27th, 2019.  I have already discussed my project in detail in the previous posts and would just summarize it here. The goal of the project is to preserve the named metadata in `.geo` files through to the tagged regions (MeshValueCollection) in dolfinx. This would help the end user to easily model with a large number of tagged regions.

A typical mesh workflow with FEniCS is as follows:

- **Creation of geometric model:**
  This task is carried out in a popular meshing software [gmsh](http://gmsh.info/).
- **Conversion of geometric model:**
  Then we convert the .msh file created by gmsh into XDMF file that could be read by DOLFIN using [meshio](https://github.com/nschloe/meshio).

- **Assignment and Application:**
  Next, we need to assign the material properties and apply loads and boundary condition. This is done with the help of “[Mesh Function](https://fenicsproject.org/docs/dolfinx/dev/cpp/d0/d94/classdolfin_1_1MeshFunction.html)” and “[Mesh Value Collection](https://fenicsproject.org/docs/dolfinx/dev/cpp/d0/db6/classdolfin_1_1MeshValueCollection.html)” in FEniCS.

The problem with the current implementation of the mesh workflow is that if a user has many different tagged regions in the mesh that s/he wishes to utilize in the simulation then there is no direct method of doing so. What was required by my GSoC project was to create a new mesh workflow that would seamlessly allow users to integrate different boundaries into their simulation. The agreed-upon method was to use the `<Information>` tag of XDMF to carry the information about various tagged regions from `gmsh` to `dolfin` via meshio.  So, in essence, we could keep the current process of working with a mesh and add methods to meshio and dolfinx.

The changes to be made to different parts of the mesh workflow is as follows:

- **Creation of geometric model:**
  There is no need for change in this part.
- **Conversion of geometric model:**
  We need to update [meshio](https://github.com/nschloe/meshio) to write the map of string tag to number tag in `<Information>` element of the XDMF file. The required changes to meshio are done and[ here is the CircleCI build corresponding to the change](https://circleci.com/gh/iitrabhi/GSoC2019/45).

- **Assignment and Application:**
  The last part of the project was to add the capabilities to dolfin to read `<Information>` tag data. A new function named read_tags() was added to the XDMFFile.cpp and the functionality of it was tested both in C++ as well as python version of the code. [Here is the CircleCI build for the same.](https://circleci.com/gh/iitrabhi/GSoC2019/45)

A basic skeleton of the new mesh workflow has been prepared in the first phase. In the next phase of programming, I would try to make these methods much more robust and add further functionality.
