---
layout: post
title: "FEniCS is incredible"
categories: misc
typora-root-url: ../../website
---
> FEniCS is incredible, but only if you have some basic understanding of the finite element method and linear algebra.

I have been using FEniCS in my research for the past two years, and I have to say that it is just pure gold in scientific computing with differential equations. It helps you think way beyond the constructions of domain-specific implementations of finite element method. You will not understand its powers until you have faced difficulties with commercial packages like Abaqus or Ansys or have tried to make your custom implementation in some programming language.

Now, don't get me wrong. Commercial packages are fantastic to carry out FE simulations using established FEM frameworks rapidly. I have used Abaqus and similar commercial packages for seven years in my research work at the Department of Atomic Energy. But, once you start to implement your custom formulations, it becomes too difficult and overly restrained to code them with the tools provided in commercial packages.

On the other hand, coding your custom implementation gives you a greater understanding of the whole structure of FEM and the coding principles. But, when you start to add complex functionalities into your code, like optimizing for speed, storage, ram usage, and many more things that are automatically handled in commercial packages, you start to feel the burden. These shortcomings, may also significantly hamper the progress of your research. Also, there is always a chance of mistake in custom implementation that might go unnoticed until you arrive at the problem.

In here comes well maintained open-source libraries that are thoroughly tested by the open-source community. They bring in the best of both worlds. You have significant control over your implementation and, at the same time, have access to inbuilt functionalities that are at par with commercial packages.

If you have ever worked with any programming languages, you might have used specific libraries to carry out a particular set of tasks for you. Like in Python, we use NumPy, which is a community maintained library to perform matrix manipulations. This library is highly efficient and provides access to various methods to work with matrices efficiently.  Numpy mimics the functionality of the popular software package MATLAB but gives it in open source so that anyone can work with it without purchasing any licenses.

> From an accuracy and reliability point of view, using well maintained open-source libraries to perform scientific computation is perfectly fine. 

From an accuracy and reliability point of view, using well maintained open-source libraries to perform scientific computation is perfectly fine. These packages follow rigorous procedures to maintain the code and are designed based on sound software design principles. Many such libraries are frequently used in research.

FEniCS is one such library, that is maintained by a passionate group of researchers. I love working with FEniCS so much because my research deals with differential equations that are still not built into the commercial packages, and coding them entirely on my own would have taken a significant amount of time. 

With FEniCS, I get the following benefits.

### 1. Easy conversion of differential equations to discrete form

Suppose I have to code my custom implementation for a particular finite element form (say the bar problem's finite element form). In that case, I have to first start from the variational-form to the weak-form and then derive the discrete-form. Once I have the discrete system, I have to code the whole system into my computer, which usually is of the form

$$ K= \int B^TDBdx$$

$$F=\int N^T F dx$$

And then, I have to code the whole FEM framework, which includes calculation of Jacobians and performing integration with Gauss Quadrature.

Whereas in FEniCS, I can give in the GDE or the weak form, and FEniCS will assemble the whole system for me. I can not stress enough how huge this is. Under the hood, this thing is not as simple as performing symbolic calculations with, say, sympy. Based on the differential equation, the mesh, and the element type you have chosen for the mesh, FEniCs will give you the assembled K matrix and F vector in just under ten lines.

This kind of coding could easily take 100's if not thousands of lines worth of code in a custom implementation. 

### 2. FEniCS is parallel

Handling a large structural system necessitates parallel computing, and anyone who has written any parallel code understands the complexity that comes with parallelization. FEniCS by design is parallel and thus allows the user to run the same code on parallel systems. Once you get a hold of parallel computing, it is tough to go serial. FEniCS has forced me to better understand parallel computation, and I have made blog posts about it that you can read [here](https://computationalmechanics.in/parallelizing-for-loop-in-python-with-mpi/) or [here](https://computationalmechanics.in/parallel-computing-with-fenics/). Going parallel on workstation could give you a significant reduction in run times with FEniCS. I have achieved upto 10X reduction, and the simulation that initially took ten days to complete in the serial was completed in one day on parallel. Achieving such a kind of scalability is very hard with custom codes.

### 3. FEniCS is designed for  HPC

FEniCS utilizes PeTSc as its linear algebra backend and supports output in XDMF format, which can handle TB's of data without corruption. Thus, as mentioned on their site, we can code our implementation on our laptops and that same code would work on the workstation or computational clusters that supports HPC. Moreover, the input and output is also supported in parallel with XDMF.

### 4. FEniCS is based on mixed programming.

This one has got to be my favorite. The computational engine of FEniCS is written in C++, and the front end supports coding in Python. Personally, I'm not too fond of compiled languages. I believe that they are not useful for non-CS researchers in the initial stage of research. Any language comes with two costs.

- Programming cost:  The amount of time that it takes to code with the language
- Runtime cost: The amount of time it takes to run the program written with the language

Compiled languages such as FORTRAN, C++ win in the runtime speed and run significantly faster than interpreted languages. But, programming with compiled languages is too tedious and takes away all the fun from coding. 

> In the initial stages of designing or developing an algorithm, you should use interpreted languages, and once the design is perfected, shift to compiled.

Since most of the heavy lifting for designing a sound parallel implementation is carried out by the FEniCS developers, researchers using it could benefit from the speed of C++ while coding their implementation in Python.

### 5. FEniCS has one of the best documentation in opensource

FEniCS is under active development and has a strong support community, which is quite active on [discourse](https://fenicsproject.discourse.group). If you ever get any kind of doubts regarding FEniCS, you can ask them in the discourse. 

 In the end, I would like to say that if anyone is at an intermediate level in their understanding of the FEM and have gone through the linear algebra courses by Dr. Gilbert Strang, they should seriously consider FEniCS for their research. You will surely get an even better understanding of the method and coding finite element.  