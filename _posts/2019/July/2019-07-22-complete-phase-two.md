---
layout: post
title: "FEniCS: Completion of phase two"
description: "The second phase of google summer of code is complete."
categories: [gsoc,coding]
typora-root-url: ../../../../website
---

Goals for the second phase of the programming were two-fold:

- Create a method to construct MVC from arrays
- Develope code to test the integration of gmsh API with new mesh workflow.

Creating MVC directly from arrays would allow bypassing saving of the mesh file onto the hard disk. The idea is to directly create MVC from pygmsh data. The functionality was achieved with and [the following PR was made](https://github.com/FEniCS/dolfinx/pull/467). Here is a snippet to create MVC from pygmsh data.

```cpp
#include <cfloat>
#include <dolfin.h>
#include "main.h"

using namespace dolfin;

int main(int argc, char *argv[]) {

    auto dummy_mesh = DummyMesh();
    auto vect = dummy_mesh.points;
    std::map<std::vector<double>, size_t> co_ord_point_map;


    dolfin::io::XDMFFile xdmf_file(MPI_COMM_WORLD, "../input/mesh_2d.xdmf");
    auto mesh_2d = std::make_shared<mesh::Mesh>(xdmf_file.read_mesh(mesh::GhostMode::none));

    auto mvc = dolfin::mesh::MeshValueCollection<int>(mesh_2d, 1,
                                                          dummy_mesh.cells, dummy_mesh.cell_data);

    return 0;
}
```

You can find the complete code [here](https://github.com/iitrabhi/GSoC2019/tree/master/Scripts/cpp/mesh-workflow-no-xdmf).

## Python + Pygmsh

```python
from pygmsh.built_in.geometry import Geometry
from pygmsh import generate_mesh
import dolfin
import dolfin.io

geom = Geometry()

mesh_ele_size = 0.5
p0 = geom.add_point([0, 0, 0], lcar=mesh_ele_size)
p1 = geom.add_point([1, 0, 0], lcar=mesh_ele_size)
p2 = geom.add_point([1, 1, 0], lcar=mesh_ele_size)
p3 = geom.add_point([0, 1, 0], lcar=mesh_ele_size)

l0 = geom.add_line(p0, p1)
l1 = geom.add_line(p1, p2)
l2 = geom.add_line(p2, p3)
l3 = geom.add_line(p3, p0)

ll = geom.add_line_loop(lines=[l0, l1, l2, l3])
ps = geom.add_plane_surface(ll)

# Tag line and surface
geom.add_physical(l0, label="BOTTOM")
geom.add_physical(l1, label="RIGHT")
geom.add_physical(l2, label="TOP")
geom.add_physical(l3, label="LEFT")
geom.add_physical(ps, label="DOMAIN")


mesh = generate_mesh(geom)
points, cells, cell_data, boundary = mesh.points, mesh.cells, mesh.cell_data, mesh.field_data

comm = dolfin.MPI.comm_world
rank = comm.Get_rank()
print('My rank is ', rank)

mesh_from_array = dolfin.cpp.mesh.Mesh(
    dolfin.MPI.comm_world,
    dolfin.cpp.mesh.CellType.triangle,
    points,
    cells['triangle'],
    [],
    dolfin.cpp.mesh.GhostMode.none)


mvc_from_array = dolfin.MeshValueCollection("size_t",
                                            mesh_from_array,
                                            1,
                                            cells["line"],
                                            cell_data["line"]['gmsh:physical'])
```

## Python + Gmsh API

```python
import gmsh
import sys
import numpy
import meshio

from dolfin import (MPI, MeshValueCollection, cpp)

# Map from gmsh to type that is also used by meshio
mapping = {1: 'line',
           2: 'triangle',
           4: 'tetra',
           8: 'line3',
           9: 'triangle6',
           11: 'tetra10',
           15: 'vertex'
           }
nodes = {'line': 2,
         'triangle': 3,
         'tetra': 4,
         'line3': 3,
         'triangle6': 6,
         'tetra10': 10,
         'vertex': 1
         }

# Generate mesh
mesh_ele_size = 0.5

gmsh.initialize(sys.argv)
gmsh.model.add('unit_square')
gmsh.model.geo.addPoint(0.0, 0.0, 0.0, mesh_ele_size, 1)
gmsh.model.geo.addPoint(1.0, 0.0, 0.0, mesh_ele_size, 2)
gmsh.model.geo.addPoint(1.0, 1.0, 0.0, mesh_ele_size, 3)
gmsh.model.geo.addPoint(0.0, 1.0, 0.0, mesh_ele_size, 4)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 6)

gmsh.model.addPhysicalGroup(1, [1], 1)
gmsh.model.setPhysicalName(1, 1, "BOTTOM")

gmsh.model.addPhysicalGroup(1, [3], 2)
gmsh.model.setPhysicalName(1, 2, "TOP")

gmsh.model.addPhysicalGroup(1, [4], 3)
gmsh.model.setPhysicalName(1, 3, "LEFT")

gmsh.model.addPhysicalGroup(1, [2], 4)
gmsh.model.setPhysicalName(1, 4, "RIGHT")

gmsh.model.addPhysicalGroup(2, [6], 6)
gmsh.model.setPhysicalName(2, 6, "DOMAIN")

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate()
nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()
dim = 3

# Reshape to get right format
points = numpy.reshape(coord, (int(coord.size / dim), dim))
element_types, element_tags, node_tags = gmsh.model.mesh.getElements()
cells = {}
cell_data = {}

# Generate cells and cell_data dict
for num, element in enumerate(element_types):
    name = mapping[element]
    num_nodes = nodes[name]
    # since nodes are numbered starting from 0
    cells[name] = node_tags[num] - 1
    cells[name] = numpy.reshape(cells[name],
                                (numpy.int(cells[name].size / num_nodes),
                                    num_nodes))
    # If mesh contains the physical group of dimension dim
    # -1 to match with gmsh api dimensions
    dim_tags = gmsh.model.getPhysicalGroups(num_nodes - 1)
    if dim_tags:
        for dim_tag in dim_tags:
            dim = dim_tag[0]
            tag = dim_tag[1]
            element_type, element_tag, node_tag = gmsh.model.mesh.getElements(dim, tag)
            for ele_num in element_tag[0]:
                element_tags[dim - 1][element_tags[dim - 1] == ele_num] = tag
        cell_data[name] = element_tags[dim - 1]

gmsh.write("unit_square.geo_unrolled")

# End Generate mesh
gmsh.finalize()

meshio.write("unit_square_gmsh.xdmf", meshio.Mesh(
    points=points,
    cells={"triangle": cells["triangle"]}))

mesh = cpp.mesh.Mesh(MPI.comm_world,
                     cpp.mesh.CellType.triangle, points,
                     cells['triangle'],
                     [], cpp.mesh.GhostMode.none)

mvc_from_array = MeshValueCollection("size_t",
                                     mesh,
                                     1,
                                     cells["line"],
                                     cell_data["line"])
```

