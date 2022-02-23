---
layout: post
title: "55 Line topology optimization code"
description: "Most compact and readable code for topology optimization."
langs: [python, c++]
year: "July 2020"
date: 2020-07-01
typora-root-url: ../../../../website
---

Even before starting my PhD, I had a strong inclination towards python and had a desire to do my work using only opensource environment. Fortunately, I found FEniCS, while searching for some documentation on PETSc. It has best of both worlds, the speed of C++ and the ease of python. I had a strong belief that coding the topology optimization problem in FEniCS could result in a highly readable code as opposed to other existing implementations. This resulted in the development of the smallest code for 2D+3D parallel implementation of topology optimization-[The 55-line code of topology optimization](https://www.researchgate.net/publication/347300347_A_55-line_code_for_large-scale_parallel_topology_optimization_in_2D_and_3D). 

![55line](/assets/images/55line.png)

The current implementation’s primary purpose is to help researchers expand their applications of TO with complex design domains on parallel HPC systems. We believe this implementation will bridge the gap between open-source codes that are limited to regular shapes and commercial packages that might not be accessible by many. Moreover, the use of HPC library PETSc - which has a very easy-to-understand interface in FEniCS - would also allow researchers to study the effect of different solvers and optimization algorithms on the TO solution. These solvers can be configured with less than ten code lines but could lead to substantial gains in terms of computational speed and accuracy. The implementation has been thoroughly tested with respect to existing literature, and the results are provided in detail in the text of the manuscript. The benefits of the current implementation are:

- This implementation is capable of handling both 2D as well as 3D domains in 55 lines.

- This implementation is not limited to unit cells since we have used the Euclidean distance matrices and actual cell centers to evaluate the filter’s weight matrix in 55 lines.

- Vectorization of most of the operations has also led to a reduction in the length of the code while increasing the code’s efficiency and maintaining its readability.

- Moreover, the code is based on variational forms instead of conventional matrix forms; it allows for a generalized approach to solving the problem and could be extended easily for dif- ferent materials models and multi-physics topology without significant change in the primary codebase.

- The code is also highly optimized and benefits from python’s readability while maintaining the speed of C++.

- The code can also handle large-scale complex structural configurations in 2D and 3D with support for parallelization and could help researchers develop new insights using parallel computing and complex design domain configurations.

Coding a custom parallel implementation in any low-level programming language is a daunting task for beginners in mechanics and is considerably time-consuming. Moreover, the parallel execution requires multiple tests to check the errors that could arise out of parallelization. It also requires a deeper understanding of the CPU/GPU architecture. In this implementation, we advocate using the high-level software tool FEniCS that handles the tedious task of parallelization and helps researchers focus on the physics of the problem. It also allows us to cut the development times by a substantial amount, as can be inferred by the length of the code presented. 
