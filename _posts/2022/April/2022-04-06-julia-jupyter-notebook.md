---
layout: post
title: "How to run julia code in jupyter notebook."
description: "It is very easy to install Julia in Jupyter."
tag: 
  - code
typora-root-url: ../../../../website
---

## The problem

I am learning the adjoint method and the open source code available is written in Julia. I wish to run it in a familiar Ipython notebook.

## The solution

- Install Julia: [Download here](https://julialang.org/downloads/)

- Open terminal and type the following commands

  ```bash
  julia
  using Pkg
  Pkg.add("IJulia")
  ```

- Close and reopen terminal after the installation is complete. Now open the notebook with the command

  ```bash
  jupyter notebook
  ```

- Inside the notebook from the kernels select Julia
  ![image-20220406081339376](/assets/images/image-20220406081339376.png) 

## References

- [Here is the complete guide](https://datatofish.com/add-julia-to-jupyter/)

