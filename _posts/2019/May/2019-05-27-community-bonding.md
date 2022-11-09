---
layout: post
title: "May 6th – May 27th, Community Bonding Period"
description: "Meeting the commutiny of GSoC2019"
categories:   [gsoc,coding]
typora-root-url: ../../../../website
---

The program starts with the community bonding period that lasts for almost a month. This is the phase of the program in which you have to know the working of your organization and get introduced to the community. You are also expected to become familiar with community practices and processes. For this phase of my GSoC journey I had planned the following:

- Get completely familiar with the code architecture of FEniCS-X, especially DOLFIN-X.
- Get familiar with the code architecture of gmsh.
- Study about the data structures of gmsh’s MSH, and, XDMF.
- Discuss with my mentor and the FEniCS community on the topic and decide on the approach to be followed for the project.

Initially, I started with trying to understand the current implementation of the mesh workflow pipeline. For a typical FEM analysis we need to follow the following steps:

- **Creation of geometric model:
  This task is carried out in a popular meshing software [gmsh](http://gmsh.info/).**
- **Conversion of geometric model:
  **Then we convert the .msh file created by gmsh into XDMF file that could be read by DOLFIN using [meshio](https://github.com/nschloe/meshio).

- **Assignment and Application:**
  Next, we need to assign the material properties and apply loads and boundary condition. This is done with the help of “[Mesh Function](https://fenicsproject.org/docs/dolfinx/dev/cpp/d0/d94/classdolfin_1_1MeshFunction.html)” and “[Mesh Value Collection](https://fenicsproject.org/docs/dolfinx/dev/cpp/d0/db6/classdolfin_1_1MeshValueCollection.html)” in FEniCS.

In gmsh different mesh entities (vertex, cells, face, volume) are grouped into “[Physical Groups](http://gmsh.info/doc/texinfo/gmsh.html#Geometry-module-1)”. Their only purpose is to assemble elementary entities into larger groups so that they can be referred to later as single entities. We create these groups so that we can refer to them later in FEniCS as an individual entity. This is somewhat similar to the creation of “sets” in ABAQUS. We use these groups to assign material properties or apply loads and boundary condition to the model.

The problem with the current implementation of the mesh processing pipeline is that it is very difficult to handle models with a large number of groups. When we create “physical groups” in gmsh we can assign them some string name (top, bottom, load_point, etc). When we create a `mesh` object from a `.msh` file in meshio, it stores these string tags in ‘field_data’ and then assigns them some number tags. When we write to an XDMF file the information regarding string tags is lost and we are left with only number tags. So there is no way for a FEniCS user to refer to those groups with the tags that s/he had assigned. There is also no map of information available in FEniCS regarding which “numerical tag” refers to which “string tag”. The main motivation for this project is to get those string tags from gmsh into DOLFIN-X so that the end user could easily refer to them in their code.

I planned to approach the problem in stages.

- **Understand the workings of GMSH
  There are excellent [screencasts](http://gmsh.info/screencasts/) provided on the official site of gmsh that got me up and running with the package in no time. Then the next task was the understand the underlying data structure of “.msh” file type. This was also easy as all the information was provided in their [documentation](http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-1).**

- **Understand the workings of meshio**
  Meshio is fairly easy to work with and the implementation is quite straight forward as explained on the offical github page. Even though I got it working on my first try, I was not able to visualize the XDMF file created by it in PARAVIEW. As soon as I opened up the file, PARAVIEW crashed with a segfault error. After digging up a bit on the internet and on the official repo on Github I found out that this was a known issue and the author of meshio [had filed it at the official repo](https://gitlab.kitware.com/paraview/paraview/issues/17945) of PARAVIEW. I tried to look up into the issue and maybe found out a bug in meshio.

  Initially, I started with the official documentation of [XDMF](http://www.xdmf.org/index.php/XDMF_Model_and_Format). The initial thing that I found out was that the definition of elements in XDMF and MSH is different and meshio is trying to map the data from MSH to XDMF. Next,  I tried to export `vtk` form meshio, and, for both 2D and well as 3D geometries it worked fine. I then exported the geometries to XDMF from Paraview. Upon inspection, I have found out that there is a difference in files produced by paraview and meshio. In the elements table, paraview has repeating element type (integer) values for elements *vertex* as well as *line*, whereas meshio is repeating the element type only for *line* and not *vertex*. I believed this mismatch resulted in a segfault.

  Upon digging a lot more into the [documentation of XDMF](http://www.xdmf.org/index.php/XDMF_Model_and_Format#XML_Element_.28Xdmf_ClassName.29_and_Default_XML_Attributes) I found out that , the element type is not repeated and the first number gives us the `TopologyType` and second gives us the `NodesPerElement` which in case of `Polyvertex` is both 1 and in case of a `Polyline` is 2, thus the repetition. This is required by the reader and is missing in meshio output.

  I got the file to work by using physical groups in my gmsh model. Performing this exercise gave me a better understanding of the workings of meshio.

- **Understanding the use of physical groups in DOLFIN-X**
  The last stage was to understand how to use the data extracted from MSH to assign material properties and boundary conditions in DOLFIN-X. For this stage, I started to work on the conversion of the code provided by Michal. I had previously worked with the Docker installation of FEniCS but that worked out of the box. With FEniCS-X, the problem started with the docker installation itself as some of the libraries necessary for the working of the code were not installed in the docker image.

  Thus, I had to understand the process of creation of a Docker image and how to build one myself.[ Getting started guide of docker](https://docs.docker.com/get-started/) is an excellent resource for understanding the principles of docker and how to build your own docker image. I would try to breakdown and explain the commands in the official docker image of DOLFIN-X in a separate post.  After understanding the concept behind docker I was able to install the required libraries in my own docker container. This was a temporary fix and eventually, my mentor Dr. Jack S. Hale made the required changes to the official Dockerfile.

  Once the container was up and running properly the next step was to read and understand the functions “[Mesh Function](https://fenicsproject.org/docs/dolfinx/dev/cpp/d0/d94/classdolfin_1_1MeshFunction.html)”(MF) and “[Mesh Value Collection](https://fenicsproject.org/docs/dolfinx/dev/cpp/d0/db6/classdolfin_1_1MeshValueCollection.html)”(MVC) that are used mark cells in the model using the information provided by the XDMF file. MVC associates the data with  entities through the corresponding *cell index* and local entity number, that we could retrieve via the `mvc.values()` function whereas MF actually creates an array of the size equal to the number of entities of a particular topological dimension and then marks the locations corresponding to the entity that we wish to mark. More details could be understood by going through [this jupyter notebook](https://github.com/iitrabhi/GSoC2019/blob/master/Notebooks/Understanding Mesh Workflow.ipynb).

  I was successful in properly converting the code snippet provided by Michal to DOLFIN-X. [Here is the final file](https://github.com/iitrabhi/GSoC2019/tree/master/Scripts/python).

The whole process gave me a fair bit of understanding of the code base. I had learned so much about opensource development, working as a team on Github, creation of Docker images, parsing XML files, about XDMF and HDF5 file format and most important I got introduced to a passionate group of people with whom I share a common interest.
