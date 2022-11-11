---
layout: post
title: "Got selected to GSoC 2019 to work on FEniCS"
description: "I am so much happy to have been officially accepted to contribute to the FEniCS project."
categories: [coding]
tags: [gsoc]
typora-root-url: ../../../../website
---

I am so much happy to have been officially accepted to contribute to the FEniCS project. For those who don’t know, [FEniCS](https://fenicsproject.org/) (Finite Element Computational Software) is a popular open-source (LGPLv3) computing platform for solving partial differential equations (PDEs). According to [wikiversity](https://en.wikiversity.org/wiki/Introduction_to_Python_for_FEniCS),

> FEniCS is an acronym with FE representing Finite Element, CS representing Computational Software, and according to Anders Logg, a Senior Research Scientist with the FEniCS Project, “ni sits nicely in the middle”. Andy Terrel also notes that the FEniCS software package was originally compiled at the [University of Chicago](http://www.uchicago.edu/index.shtml), whose mascot is a phoenix, which likely inspired the name.

# Why FEniCS?

I am a second year Ph.D. student at the Indian Institute of Technology Roorkee. My primary areas of interest are the finite element method (FEM) and isogeometric analysis (IGA). I have worked with many different FEM software packages (such as ABAQUS, STAAD, SAP, SAFE, etc) in the past.

For my Ph.D. work, I had decided to work with some open-source package as that would help me better understand the method and how do we actually go from differential equations to numerical implementation. Initially, I started with some small open source codes written by individuals to understand different parts of FEM.

Then, I came across FEniCS, and really liked the simplicity of the way in which we could code a variational form. Also, it had a fairly vast set of documentation and tutorials available that helped me get started. I was hooked to the package and started implementing some simple FEM problems with.

I have to admit that the initial learning curve is quite steep especially if you are used to some GUI based application, but as you get confident you can explore many different variational forms without the need of explicitly coding different operators.

I came to know about GSoC through the official website of FEniCS where they had announced that they will be participating in Google Summer of Code 2019 under the umbrella of NumFOCUS. I was really excited about this opportunity as this would help me to better understand the working of this software and moreover give me an opportunity to work with its creators.

# The application process

The organizers had listed three different project ideas on their [GitHub webpage](https://github.com/FEniCS/gsoc/blob/fenics/ideas-2019/2019/ideas-list-fenics.md). Out of those the project idea “*Gmsh/XDMF/DOLFIN mesh pipeline*” intrigued me the most as this was something that was directly related to my Ph.D. work. As per the [official guidelines](https://github.com/numfocus/gsoc/blob/master/CONTRIBUTING-students.md) I had to primarily do the following:

1. Write a proposal
2. Try to fix some existing issue on the official repo and make a pull request.

NumFocus has an [official template](https://github.com/numfocus/gsoc/blob/master/templates/proposal.md) which the student has to follow. To write a good proposal, first of all, you have to understand the project and try to execute a part of it so that you have an understanding of what it is that you are required to do. Working on the project beforehand also helps you to come up with realistic deadlines in your proposal. Thus, before writing the proposal I tried to work with Gmsh, Meshio, Dolfin, and Paraview on a few problems.

I went through the official documentation of XDMF and Gmsh and tried to understand the existing issues with the mesh processing pipeline. This gave me a fairly good idea as to what is required to be done in the project. The most difficult part was trying to understand the code base of DOLFIN. I made the first draft and then moved towards the second step and that was to make a pull request in the official repo.

On a daily basis, I work mostly with Python and have a good understanding of the language. I have not used C++ in any of my previous projects. Thus, the second part of the application was quite tough as DOLFIN (the core component of FEniCS) is written in C++. Initially, I tried to read some books on C++ and made some small codes using C++.

Parallely, I was also working with meshio. I found a small bug in meshio (luckily written in python), corrected it and made the[ pull request](https://github.com/nschloe/meshio/pull/374). I found a similar bug in the testing module of dolfinx, corrected it and made the [pull request](https://github.com/FEniCS/dolfinx/pull/375/files). Both of these pull request got merged into the main code. I was so much happy as this was the first time I had contributed to an open source project.

After sharing the proposal with the organizing committee, I got approval and then uploaded the proposal to the GSoC website. And finally, on, 7th May I got a mail from GSoC stating that my proposal has been accepted.
