---
layout: post
title: "How to visualize degrees of freedom in FEniCS."
tag: 
	-  fenics
typora-root-url: ../../website
---

When we create or import a mesh in FEniCS, it automatically creates the degree of freedom table for us. Now that table is based on multiple considerations, which I do not understand, but the thing to know is that, degrees of freedom are not directly derived from the vertex numbers. Thus, there is no logical relation between vertex numbers and dof numbers. But, fortunately, FEniCS provides us with the vertex to degree of freedom table, i.e. for a particular vertex number we can find out the degree of freedom numbers.

```python
from dolfin import *

mesh = BoxMesh(Point(0.,0.,0.),Point(2,1,1), 2, 2, 2)
V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)
v2d=vertex_to_dof_map(V)
v2d = v2d.reshape((-1, mesh.geometry().dim()))
```

The array `v2d` gives us the degrees of freedom corresponding to a particular vertex. 

![image-20210125142158802](/assets/images/image-20210125142158802.png)

To see the degrees of freedom in paraview we can use the following code snippet

```python
from dolfin import *

mesh = BoxMesh(Point(0.,0.,0.),Point(2,1,1), 2, 2, 2)
V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)

num_dof = mesh.num_vertices()*V.dofmap().num_entity_dofs(0)
dof_map = Function(V,name="dof")
dof_map.vector()[:] = [int(i) + 1 for i in np.linspace(0,num_dof-1,num_dof)]

with XDMFFile("dof.xdmf") as xdmf:
    xdmf.write(dof_map)
```

The thing to remember here is that, when we call the `vector()` method on any function, the return vector is ordered as per the `degrees of freedom`. Thus, we can just put in the row number of the vector as degree of freedom number and visualize it in paraview. The `Id` in paraview is the vertex number assigned by FEniCS.![image-20210125141950089](/assets/images/image-20210125141950089.png)

> Note that in this post I have added +1 to the degrees of freedom. This is done just to compare the output with that of MATLAB. If you wish to use this with FEniCS than delete the +1.

```python
from dolfin import *

mesh = BoxMesh(Point(0.,0.,0.),Point(2,1,1), 2, 2, 2)
V = VectorFunctionSpace(mesh, 'Lagrange', degree=1)

num_dof = mesh.num_vertices()*V.dofmap().num_entity_dofs(0)
dof_map = Function(V,name="dof")
dof_map.vector()[:] = [int(i) for i in np.linspace(0,num_dof-1,num_dof)]

with XDMFFile("dof.xdmf") as xdmf:
    xdmf.write(dof_map)
```

