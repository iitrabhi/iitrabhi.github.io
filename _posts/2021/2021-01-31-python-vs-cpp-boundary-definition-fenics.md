---
layout: post
title: "Difference between python and C++ definition of boundaries in FEniCS."
categories: [coding]
tags: [fenics]
typora-root-url: ../../../website
---

In FEniCS there are two ways of defining the sub-domains or different zones in a mesh. One is the python way and other is the C++ way. In this post I will find the speed difference between the two. I will follow along [this article](https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1005.html) from the official FEniCS website. 

```python
from fenics import *

lang = "p"

# Create mesh and define function space
mesh = UnitSquareMesh(1000, 1000)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

if lang == "p":
	print("using python")
	def boundary(x, on_boundary):
	    return on_boundary
else:
	print("using c++")
	boundary = CompiledSubDomain("on_boundary")

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)
```

For expressions, I went throught the discourse post [here](https://fenicsproject.discourse.group/t/cpp-based-expression/929).