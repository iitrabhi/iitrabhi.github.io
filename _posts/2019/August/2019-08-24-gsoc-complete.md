---
layout: post
title: "GSoC complete: FEniCS-The mesh workflow"
description: "This is the final post for my google summer of code 2019"
categories: 
  - [gsoc,coding]
typora-root-url: ../../../../website
---

This is my final post for the GSoC2019 program. The primary goal of the project was to ensure that the meshing package of choice gmsh, DOLFIN, and the preferred visualization package, [Paraview](http://paraview.org/) work seamlessly together. The intention was to make improvements to the process of preserving the information about tagged regions of the mesh when importing it to DOLFIN. Two approaches were finalized for the project.

- Preserve the string tags when converting from `.msh` to `.xdmf`.
- Add a constructor to the class MeshValueCollection to support its creation from primitive arrays.

The targets were achieved with the following pull requests.

### MESHIO

- **PR #425**

  -Add methods to read and write field_data to XDMF – 

  merged on Aug 9

  - This PR adds methods to preserve the mapping between the string tags and int tags when converting from `.msh` to `.xdmf`. The idea is to store this mapping inside the `<Information>` element of XDMF.

### DOLFINX

- **PR #439**

   – Add function to read Information tag data in XDMF – 

  open

  - This PR extends the functionality of XDMFFile interface with methods `read_information()` and `write_information()`. The methods are designed to read and write `<Information>` tag data of XDMF file format. The syntax is inline with XDMF standard and works well with PARAVIEW.

- **PR #467**

   – Method to construct MeshValueCollection from arrays – 

  open

  - This PR extends the functionality of MeshValueCollection class with constructor that supports its creation from primitive arrays.

## Demonstration

We can solve many different forms of PDE on simple as well as complex domains in FEniCS. We have many different prebuilt geometries in FEniCS that helps new users to get up and running with simple FEniCS classes and methods. Even though these built-in meshes provide various methods for their construction and refinement, they are limited to simple shapes. To work with complex geometrical structures it is recommended that the user follows the following mesh workflow.

![image-20220211110737972](/assets/images/image-20220211110737972.png)

In this post, I am going to walk you through the process of solving a PDE in FEniCS. I will demonstrate the whole mesh workflow in detail. The domain under consideration in a unit square as presented below. The domain is made up of two different materials and we need to apply different boundary conditions on all of its edges.

![poisson_subdomain](/assets/images/poisson_subdomain.png)

## Creating mesh in Gmsh

There are excellent tutorials and documentation available on the official website of Gmsh. In this small screencast, I have described the process of creation and tagging of the above mesh using gmsh.

https://youtu.be/C4W_28NjF58

## Convert mesh to XDMF using meshio

The .msh file created by gmsh could be converted to .xdmf by using the package meshio. The package could be easily installed by the following command:

```
pip install meshio
```

Once you have the package installed you can use the following command to convert the mesh to xdmf. Right now FEniCS does not support mixed topologies so you have to individually export mesh entities of different dimension to different XDMF files. Thus for the current mesh, we need to export one mesh of 2D `triangle` elements and the other of 1D ‘line’ elements.

```python
mesh_of_triangles = meshio.Mesh(points=points[:, :2],
                                cells={'triangle': cells['triangle']},
                                cell_data={'triangle': 
                                    {'subdomain': 
                                        cell_data['triangle']
                                        ['gmsh:physical']}},
                                field_data=field_data) 

meshio.write("poisson_subdomain_triangle.xdmf", mesh_of_triangles )

mesh_of_lines = meshio.Mesh(points=points[:, :2],
                                cells={'line': cells['line']},
                                cell_data={'line': {'subdomain':
                                    cell_data['triangle']
                                    ['gmsh:physical']}},
                                field_data=field_data) 

meshio.write("poisson_subdomain_line.xdmf", mesh_of_lines )
```

## Import XDMF to FEniCS

Once we have the XDMF files we can import them to FEniCS as follows:

```python
with XDMFFile(MPI.comm_world,
              "poisson_subdomain_triangle.xdmf") as xdmf_infile:
    mesh = xdmf_infile.read_mesh(cpp.mesh.GhostMode.none)
    tag_info = xdmf_infile.read_information_int()
    print(tag_info)
```

FEniCS provide us with two classes `MeshFunction` and `MeshValueCollection` that help us to mark different mesh entities (vertex, line, facet, cell). The way in which the `MeshFunction` class works is that it assigns a marker to every element of a particular dimension (vertex = 0, line = 1, etc.) and then use that marker to differentiate between different regions. `MeshValueCollection` class differs from the `MeshFunction` class in two ways. First, data does not need to be associated with all entities (only a subset). Second, data is associated with entities through the corresponding cell index and local entity number (relative to the cell), not by global entity index, which means that data may be stored robustly to file. You can get a much better understanding of them by exploring them individually as I have done [here](https://github.com/iitrabhi/GSoC2019/blob/master/Notebooks/Understanding MeshFunction and MeshValueCollection.ipynb).

## Visualize the solution

We could visualize the solution using opensource package paraview. Once again the filetype of choice is XDMF and we can easily write the solution to the file using the **write** method of class **XDMFFile**.

```python
with XDMFFile(dolfin.MPI.comm_world, "output.xdmf") as xdmf_outfile:
    xdmf_outfile.write(u)
```

You can find the complete code [here](https://github.com/iitrabhi/dolfinx/blob/iitrabhi/mvc-xdmf/python/demo/poisson-subdomain/demo_poisson_subdomain.py).

## Summary

The mesh workflow of FEniCS is not as straight forward as that of some commercial packages like Abaqus or Ansys. But, FEniCS provides you with greater flexibility in creating a computational model and would also help you understand the mathematics behind FEM. I believe that if understood and used properly, FEniCS could help researchers code and test out their mathematical models with a speed and ease that is just simply not possible with other packages.
