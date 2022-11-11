---
layout: post
title: "How to parallelize FEniCS in jupyter notebook."
description: "I have been trying to debug a portion of my FEniCS code that required parallelization for faster runtime for debugging. Just found out that we can parallelize in jupyter itself 🤯."
categories: [fenics,coding]
typora-root-url: ../../../../website

---

## The problem

Designing and debugging a FEniCS code in a jupyter notebook is effortless. I love the interactive environment and the markdown support that allows writing equations inside the notebook; this helps debug the application for the correctness of the code. But, recently, I had to debug a portion of the code that was taking far too much time to run. Somehow, I had to parallelize that portion of the code itself in the notebook.

## The solution

- Install the package [ipyparallel](https://ipyparallel.readthedocs.io/en/latest/)

- The tutorial presented in the documentation is quite straighforward. You need to only make a function out of the portion of the code that you need to parallelize and  then call it using `ipyparallel`

  ```python
  import ipyparallel as ipp
  
  def mpi_example():
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      return f"Hello World from rank {comm.Get_rank()}. total ranks={comm.Get_size()}"
  
  # request an MPI cluster with 4 engines
  with ipp.Cluster(engines='mpi', n=4) as rc:
      view = rc.broadcast_view()
      r = view.apply_sync(mpi_example)
      print("\n".join(r))
  ```

## The execution

From the official documentation, it seems pretty straightforward. Parallelization is as simple as setting the number of processors to parallelize. But, from experience, I know that it is never this simple for custom applications. Anyhow, as a starting point, if I am not trying to get anything back from the parallel run, it should just work. Here is my first try with the Poisson equation.

```python
import ipyparallel as ipp
import time
```
```python
def poisson():
    
    import dolfin as df    
    comm = df.MPI.comm_world
    
    n = 600
    mesh = df.RectangleMesh(comm, df.Point(0,0), df.Point(1,1), n, n)
    
    V = df.FunctionSpace(mesh, "Lagrange", 1)
    # Define Dirichlet boundary (x = 0 or x = 1)
    boundary = df.CompiledSubDomain("near(x[0], 0) || near(x[0],1)")
    # Define boundary condition
    u0 = df.Constant(0.0)
    bc = df.DirichletBC(V, u0, boundary)
    # Define variational problem
    u = df.TrialFunction(V)
    v = df.TestFunction(V)
    f = df.Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)",degree = 2)
    g = df.Expression("sin(5*x[0])",degree=2)
    a = df.inner(df.grad(u), df.grad(v))*df.dx
    L = f*v*df.dx + g*v*df.ds
    # Compute solution
    u = df.Function(V)
    df.solve(a == L, u, bc)
    with df.XDMFFile("output.xdmf") as xdmf:
        xdmf.write(u)
    return mesh.num_cells()
```
```python
# request an MPI cluster with 4 engines
start = time.time()
with ipp.Cluster(engines='mpi', n=6) as rc:
    view = rc.broadcast_view()
    r = view.apply_sync(poisson)
end = time.time()
print("Run time in seconds ----------------",int(end-start))
```
```python
print("Elements with each processor: ",r)
```
``` 
Elements with each processor:  [119007, 119491, 122701, 117966, 120324, 120511]
```
As we can see from the output, each processor has a different portion of the mesh and thus a different number of elements. 
Right now, I am unable to see significant computational gains with parallelization. But, anyhow, this experiment is half successful.

## Reference

[Using IPython for parallel computing](https://ipyparallel.readthedocs.io/en/latest/)