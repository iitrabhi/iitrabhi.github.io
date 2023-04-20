---
layout: post
title: "Demystifying High-Performance Computing for Finite Element Analysis Simulations"
description: "This article aims to explain the fundamentals of High-Performance Computing (HPC) and its role in Finite Element Method (FEM) simulations"
categories:  [coding, simulation]
tags: []
typora-root-url: ../../../../website
---

High-performance computing (HPC) has become an increasingly important tool in modern scientific research, engineering, and other fields that require the processing of large amounts of data and complex calculations. Finite Element Method (FEM) simulations, which are used to model and analyze physical systems and structures, are one area where HPC has had a significant impact.

As someone who's interested in FEM simulations, I've found that understanding HPC can be pretty tough, even for experts. So, in this article, I'll explain the basics of HPC and how it's used in FEM simulations. By the end, you'll have a better idea of how HPC can help you run more accurate and efficient simulations.

## How do we simulate reality?

In engineering we simulate reality by creating mathematical models that describe the behavior of real-world phenomena. These models can be based on fundamental physical laws and principles, as in the case of physics-based simulations used in scientific research and engineering. Alternatively, they can be based on statistical models or machine learning algorithms that learn patterns and relationships from data.

Once a mathematical model is created, it can be implemented on a computer and solved using matrix algebra. Matrix algebra is a powerful tool for solving systems of linear equations, which are commonly used in many types of simulations. In matrix algebra, equations are represented using matrices, which are arrays of numbers that can be manipulated using mathematical operations such as addition, subtraction, multiplication, and inversion.

By using matrix algebra to solve mathematical models, we can simulate the behavior of real-world phenomena with a high degree of accuracy and efficiency. However, creating accurate and reliable models requires a deep understanding of the underlying physical principles and a careful consideration of the data used to inform the model.

## Why simulation is hard?

Achieving high-fidelity simulations of reality is a difficult task due to the trade-off between accuracy and computational cost. In order to accurately model complex physical phenomena, we need to use Finite Element Analysis (FEA) simulations that have a large number of degrees of freedom. These degrees of freedom represent the number of parameters used to describe the system being analyzed, such as the position, velocity, and acceleration of each point in a structure. They also represent the size of the matrices used in analysis.

However, increasing the number of degrees of freedom also increases the computational cost of the simulation. This is because more degrees of freedom require more memory and processing power to accurately model the system. As a result, running high-fidelity simulations can take a long time and require significant computing resources.

To overcome this challenge, researchers have developed various techniques to optimize the computational efficiency of FEA simulations. These include the use of parallel computing, which allows multiple processors to work together to solve a single simulation, and adaptive meshing, which dynamically adjusts the resolution of the simulation based on the complexity of the system being analyzed. By using these techniques, it's possible to achieve high-fidelity simulations of reality while minimizing computational cost and run time.

## What is high performance computing?

Normal computing refers to the use of general-purpose computing systems, such as desktops, laptops, or mobile devices, to perform everyday tasks such as browsing the internet, sending emails, creating documents, and running basic applications. The computers are designed to handle light weight tasks and sometimes basic parallel computations.

Whereas high performance computing is a type of computing that is designed to perform complex calculations and handle large amounts of data at very high speeds. It typically involves the use of specialized hardware and software optimized for parallel processing, which allows multiple computations to be performed simultaneously. HPC systems are used in a wide range of fields, including scientific research, engineering, and finance, among others.

## How is HPC different from normal computing?

HPC can be thought of as a large dam that is used to control and manage the flow of water in a river. It is designed to manage and control large amounts of data and computing power, much like a dam manages the flow of water. Just as a dam is made up of multiple components that work together to control the flow of water, an HPC system is made up of multiple components such as processors, memory, storage, and networking infrastructure that work together to manage and control large-scale computational workloads.

Workstations, on the other hand, can be compared to a pipeline system that delivers a large volume of water to multiple users at once. Workstations are designed for heavy-duty computing tasks that require high-performance computing power, memory, and storage resources. They are typically desktop computers with multiple cores, high-speed storage, and powerful graphics cards that can handle demanding tasks such as video editing, 3D rendering, and scientific simulations.

Laptops can be thought of as portable water bottles that can be carried around and used to deliver a small amount of water to an individual user. Like a water bottle, a laptop is designed to be compact, lightweight, and easy to carry around. It has limited computational power and memory resources, but it's suitable for performing everyday tasks such as browsing the web, sending emails, and light productivity work.

## What are some HPC linear algebra libraries?

HPC linear algebra libraries are software packages that provide optimized implementations of linear algebra algorithms for HPC architectures. These libraries are designed to take advantage of the parallelism inherent in HPC systems to perform matrix computations efficiently on large-scale data sets.

Some popular HPC linear algebra libraries include:

1.  **BLAS (Basic Linear Algebra Subprograms):** a standard library that provides a set of low-level linear algebra operations such as matrix multiplication, matrix-vector operations, and vector-vector operations.
2.  **LAPACK (Linear Algebra Package):** a higher-level library that builds on top of BLAS and provides routines for solving systems of linear equations, eigenvalue problems, and singular value decomposition.   
3.  **ScaLAPACK (Scalable Linear Algebra Package):** a parallel implementation of LAPACK designed to run efficiently on distributed-memory systems.
4.  **PETSc (Portable, Extensible Toolkit for Scientific Computation):** a library that provides a set of scalable linear and nonlinear solvers, including direct solvers based on LU and Cholesky factorizations and iterative solvers such as Conjugate Gradient.
5.  **Trilinos:** a collection of libraries that provide scalable algorithms for solving linear and nonlinear problems, including linear algebra, optimization, and partial differential equations.

These libraries are widely used in scientific and engineering applications, including computational fluid dynamics, weather modeling, and molecular dynamics simulations. By leveraging these libraries, researchers and practitioners can achieve high performance and scalability in their numerical simulations and data analytics.

## Can I use HPC libraries on my laptop?

HPC libraries such as BLAS, LAPACK, and PETSc can be installed and run on laptops, but the performance of these libraries will be limited by the processing power and memory capacity of the laptop.

HPC libraries are like high-performance engines in cars. Just like how high-performance engines are designed to produce more power and speed than standard engines, HPC libraries are designed to perform complex computations more efficiently and quickly than standard software libraries.

Installing HPC libraries on a laptop is like installing a high-performance engine in a small car. While the high-performance engine may be capable of producing a lot of power, the size and weight of the car will limit its performance. Similarly, while HPC libraries can be installed on laptops, the limited processing power and memory capacity of the laptop will limit the performance of the libraries.

In contrast, using HPC libraries on a dedicated HPC cluster is like using a high-performance engine in a sports car. The sports car is designed to handle the power and speed of the high-performance engine, allowing it to achieve peak performance. Similarly, dedicated HPC clusters are designed to handle the complex computations performed by HPC libraries, allowing them to achieve high performance and scalability.

Nonetheless, even though the performance of HPC libraries may be limited on a laptop, they can still be useful for smaller-scale computations or for prototyping and testing purposes. For example, if you are developing a new algorithm that requires linear algebra computations, you could start by testing it on your laptop using HPC libraries before moving to a more powerful HPC system for larger-scale simulations.

HPC systems can be very expensive to build, operate, and maintain. The cost of purchasing and maintaining the necessary hardware, such as high-end CPUs, GPUs, and networking equipment, can be significant. Additionally, the cost of energy and cooling required to run HPC systems can be high, particularly for larger systems. By starting with a smaller-scale system, you can test and refine your application before investing in a larger, more expensive HPC system.

## Message passing interface

MPI stands for Message Passing Interface which is meant to be a tool to pass messages among various processes to carry out certain task. Processes corresponds to the physical cores available in your system. You can look for the number of cores available in your own computer by going into the task manager (On a windows computer) and looking for it in the Performance tab.

![](assets/images/Pasted%20image%2020230420123236.png)

On my system, I have 6 physical cores. This means that I can parallelize my code to a maximum of 6 processes. Now the question is – how to achieve that with python?

Before that we need to understand a bit about MPI and its terminology.

**MPI_COMM_WORLD:** MPI uses objects called communicators, which are a collection of processes. The default communicator is called MPI_COMM_WORLD and it encompasses all the processes available. In my case the MPI_COMM_WORLD will look something like this:

![](assets/images/Pasted%20image%2020230420123247.png)

**World Size:** This would tell the program about the number of processors available in the world.

**Processor Rank:** This is a unique number assigned to each processor inside the world. The numbering starts from 0 as shown in the above figure.

**Barrier:** As the name suggests this acts as a barrier in the parallel execution. This forces MPI to execute all the commands before the barrier by all the processes. This is required in the situation where you need a certain variable generated by the program in all processes. This will become clear in the example presented below.

In MPI, processes are organized into groups called communicators. Communicators define a group of processes that can communicate with each other using MPI messages. Think of communicators as a group of people working together to solve a puzzle.

Each process in MPI is identified by a unique rank within the communicator. The rank of a process is like a person's seat at the table while working on a puzzle. Just as each person has their own seat at the table, each process has its own rank within the communicator.

For example, imagine that you are working on a puzzle with a group of friends. You are all part of the same communicator, which is like the puzzle-solving group. Within the group, each person has a unique rank that identifies their seat at the table. If you are sitting in seat 1, your rank within the communicator is 1.

Processes use their ranks to communicate with each other. For example, one process might send a message to another process with a specific rank. This message might contain data that the receiving process needs to complete its part of the computation.

## How do you run a script in parallel?

First we need to install a few packages for MPI to work. I have compiled all of them into an image that you can use by installing docker on your system.

- Go to https://www.docker.com/ and download and install docker desktop
- After installation open terminal/command prompt in you computer and then type the following command
```
docker run -v host_system_path:/root/ -w /root/ -it iitrabhi/fenics
```

*Note: you should replace the variable `host_system_path` with the path of the folder that contains your code. e.g. If `D:\Codes` contains your code then to start the command line interface you have to run:*

```
docker run -v D:\Codes:/root/ -w /root/ -it iitrabhi/fenics
```

You can download the example code at https://github.com/iitrabhi/parallel-fenics.

You can run any python script with MPI by typing the following command in terminal:

```
mpirun -np 6 python3 main.py
```

This will tell MPI to run test-parallel.py on 6 processors. The thing to understand here is that even though you are running the program on 6 processors, you are not actually doing parallel computations. You are just doing the same computation 6 times. To actually do parallel computations, you need to manually split the code to work parallelly. When you type the above command the system creates 6 different copies of the program file and sends it to 6 different processes.

![](assets/images/Pasted%20image%2020230420123306.png)

Thus we need to identify the processor number inside the program and execute the commands accordingly. We can identify the processor number by first getting a handle to the world communicator by using command

```python
comm = MPI.COMM_WORLD
```

and then get the size of the world and the rank of the processor by using the commands

```python
rank = comm.Get_rank()
size = comm.Get_size()
```

Based on this information we can modify our logic to run on multiple processors. This simple program sums the numbers from a to b and gives us the result. This logic is parallelizable. We can split the whole domain of calculation into the number of processors available and then add the numbers in that domain. Finally, we can add the results of all the processors. The program below is self-explanatory and you can run the same on your machine with the help of `mpirun` command.

```python
import numpy
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

a = 1
b = 10000000

num_per_rank = b // size # the floor division // rounds the result down to the nearest whole number.
summ = numpy.zeros(1)

temp = 0
lower_bound = a + rank * num_per_rank
upper_bound = a + (rank + 1) * num_per_rank
print("This is processor ", rank, "and I am summing numbers from", lower_bound," to ", upper_bound - 1, flush=True)

comm.Barrier()
start_time = time.time()

for i in range(lower_bound, upper_bound):
    temp = temp + i

summ[0] = temp

if rank == 0:
    total = numpy.zeros(1)
else:
    total = None

comm.Barrier()
# collect the partial results and add to the total sum
comm.Reduce(summ, total, op=MPI.SUM, root=0)

stop_time = time.time()

if rank == 0:
    # add the rest numbers to 1 000 000
    for i in range(a + (size) * num_per_rank, b + 1):
        total[0] = total[0] + i
    print("The sum of numbers from 1 to 1 000 000: ", int(total[0]))
    print("time spent with ", size, " threads in milliseconds")
    print("-----", int((time.time() - start_time) * 1000), "-----")
```

The only thing to notice is that the input to the loop changes according to the the processor number (rank). Thus instead of looping `b` times, each processor has to loop only `num_per_rank` times.

Running the above script on a single core result in the following output:

```
➜  pre-processing (master) time mpirun -np 1 python3 test-parallel.py
This is processor  0 and I am summing numbers from 1  to  60000000
The sum of numbers from 1 to 1 000 000:  1800000030000000
time spent with  1  threads in milliseconds
----- 5851 -----

real    0m6.824s
user    0m6.800s
sys     0m0.010s
```

and running the same code on 6 cores results in the following output

```
➜  pre-processing (master) time mpirun -np 6 python3 test-parallel.py
This is processor  0 and I am summing numbers from 1  to  10000000
This is processor  2 and I am summing numbers from 20000001  to  30000000
This is processor  3 and I am summing numbers from 30000001  to  40000000
This is processor  4 and I am summing numbers from 40000001  to  50000000
This is processor  1 and I am summing numbers from 10000001  to  20000000
This is processor  5 and I am summing numbers from 50000001  to  60000000
The sum of numbers from 1 to 1 000 000:  1800000030000000
time spent with  6  threads in milliseconds
----- 1668 -----

real    0m1.945s
user    0m11.250s
sys     0m0.200s
```

Thus we have achieved almost a $3.5 \times$ boost in speed by running the code in parallel. Now, this was a trivial example but in real calculations, we can expect greater speed boosts.

## Running finite element simulation in parallel using open source HPC package FEniCS

```python
from dolfin import *
set_log_level(50)

comm = MPI.comm_world
rank = MPI.rank(comm)


def mprint(*argv):
    if rank == 0:
        out = ""
        for arg in argv:
            out = out + str(arg)
        # this forces program to output when run in parallel
        print(out, flush=True)

def epsilon(v):
    return sym(grad(v))

def sigma(v):
    return lmbda*tr(epsilon(v))*Identity(2) + 2.0*mu*epsilon(v)


lx, ly, mul = 10, 1, 3
nx, ny = int(lx * mul), int(ly * mul)
mesh = RectangleMesh(comm, Point(0., 0.), Point(lx, ly), nx, ny, "crossed")

E, nu, rho = Constant(2e11), Constant(0.3), Constant(7850) # E=N/m2 and rho=kg/m3
mu, lmbda = E/2/(1+nu), E*nu/(1+nu)/(1-2*nu)
lmbda = 2*mu*lmbda/(lmbda+2*mu)

rho_g = rho * 9.8 * 100
f = Constant((0, -rho_g))
U = VectorFunctionSpace(mesh, 'Lagrange', degree=1)
u, v = TrialFunction(U), TestFunction(U)
a = inner(sigma(u), epsilon(v))*dx
l = inner(f, v)*dx

mprint("DoF: ", (U.dim()))

left = CompiledSubDomain("near(x[0], 0.0, tol) && on_boundary", tol=1e-14)

bc = DirichletBC(U, Constant((0., 0.)), left)
u_sol = Function(U, name="displacement")
problem = LinearVariationalProblem(a, l, u_sol, bc)
solver = LinearVariationalSolver(problem)

prm = solver.parameters
prm["linear_solver"] = 'cg'
prm["preconditioner"] = 'hypre_euclid'
solver.solve()

mprint("Maximal deflection    : ", -
       round(MPI.min(comm, u_sol.vector().min()), 3))
mprint("Beam theory deflection: ", round(3*rho_g*lx**4/2/E/ly**3, 3))

file_results = XDMFFile(comm, "output/elasticity_results.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True
file_results.parameters["rewrite_function_mesh"] = False
file_results.write(u_sol)

# solvers = (
#     "bicgstab",
#     "cg",
#     "default",
#     "gmres",
#     "minres",
#     "mumps",
#     "petsc",
#     "richardson",
#     "superlu",
#     "tfqmr",
#     "umfpack",
# )

# preconditioners = (
#     "amg",
#     "default",
#     "hypre_amg",
#     "hypre_euclid",
#     "hypre_parasails",
#     "icc",
#     "ilu",
#     "jacobi",
#     "none",
#     "petsc_amg",
#     "sor",
# )

# linesearch = ("basic", "bt", "cp", "l2", "nleqerr")

# list_timings(TimingClear.clear, [TimingType.wall])


```