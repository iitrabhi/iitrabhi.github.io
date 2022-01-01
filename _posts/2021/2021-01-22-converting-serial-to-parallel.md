---
layout: post
title: "How to convert serial FEniCS code to parallel."
tag: 
  - fenics
typora-root-url: ../../../website
---

FEniCS supports parallel computation out of the box, but that requires certain considerations from the side of the programmer. To convert a FEniCS code into parallel, we have to add the following to the code

## 1. Get the MPI communicator

FEniCS utilizes MPI for parallelization. The command `MPI`  is built into dolfin. When we import  dolfin we import MPI with it. To understand MPI in detail [read this post](https://computationalmechanics.in/parallelizing-for-loop-in-python-with-mpi/). We start by getting the MPI communicator and the rank of the current processor on which the code is placed.

```python
comm = MPI.comm_world
rank = MPI.rank(comm)
```

## 2. Define MPI IO functions

The input output operations could be a bit tricky in MPI. This is because if you will send the code to 5 processors with a print statement, the print will happen 5 times on the console. To prevent this whenever we wish to print something we should use the following `mprint` command.

```python
def mprint(*argv):
    if rank == 0:
        out = ""
        for arg in argv:
            out = out + str(arg)
        # this forces program to output when run in parallel
        print(out, flush=True)
```

writing to XDMF file is handled by default in FEniCS, but if we wish to write to a text file then we use the following method.

```python
def mwrite(filename, my_list):
    MPI.barrier(comm)
    if rank == 0:
        with open(filename, "w") as f:
            for item in my_list:
                f.write("%s" % item)
```

The input for the above command should be a list of strings in the following format.

```python
my_list = [
    "force \t displacement \n",
    "1 \t 2 \n",
    "2 \t 4 \n"
]
```

##  3. Keep in mind the mesh partition

When we use MPI in FEniCS, the mesh is divided into `n` number of parts, where `n` is the number of processors. Then, each processors gets its own part of the mesh. The continuity of the field variables between different parts of the mesh is maintained by FEniCS. But, while preparing a parallel code, we need to keep in mind that the operations that we perform will happen only on a part of mesh. 

One of the most common mistake that I have made repeatedly was to ask the value of a function at a particular point from FEniCS in parallel. This will always throw an error, because, that point will exist only with  a single processor and all other processors will not have any knowledge of that point. Thus we have to either `allow_extrapolation` in our code or catch the exception and handle it in other processors.

![image-20210122140442143](/assets/images/image-20210122140442143.png)

In the above image the white part is with first processor  and the blue part is with second processor. This mesh is generated from a file that was executed on two processors.

Also, whenever we are retrieving any kind of `max`, `min`, `sum` values from the code, we need to ask the communicator also to do the same operation, i.e. we need to type

```python
MPI.max(comm, max(array))
```

`max(array)` will give the maximum value of the array from a single processor, then `MPI.max()` will give us the maximum of all the processors.