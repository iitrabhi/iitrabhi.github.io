---
layout: post
title: "How to modularize fenics code."
categories: fenics
typora-root-url: ../../website
---
To date, most of the codes that I have designed using FEniCS were based on flat style. Those were all contained in a single file. FEniCS documentation includes an excellent post on how to [improve the poison's solver](https://fenicsproject.org/pub/tutorial/html/._ftut1016.html#ch:poisson0:impl2). I would recommend anyone reading this post first to go and understand that post.

```python
from fenics import *

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
usol = Function(V)

problem = LinearVariationalProblem(a, L, usol, bc)
solver = LinearVariationalSolver(problem)
solver.solve()
```

A typical FEniCS code starts with definition of a `FunctionSpace` and all the other entities that are defined on this space, contains the original definition of the `FunctionSpace`

![untitled@2x (2)](/assets/images/untitled@2x%20(2).png)

When we are modularizing the code we have to keep in consideration the fact that any kind of further manipulation of the variables will require the original definition of the function space. **Thus, when we put the above code into some python function we have to return all the variables that we wish to use outside the function**. 