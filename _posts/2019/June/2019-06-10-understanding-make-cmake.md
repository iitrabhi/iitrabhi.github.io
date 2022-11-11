---
layout: post
title: "Understanding make and cmake"
description: "Understanding the compilation process of C++"
categories: [coding]
tags: [gsoc]
typora-root-url: ../../../../website
---

Up until this point in my life, I have never worked on a project in C++. Starting work on the cpp version of FEniCS was quite difficult as  I had no idea how even a simple program works in cpp. Thus, I planned to first properly learn the syntax of cpp and then try and understand the workings of FEniCS cpp code. As is with most of the programming languages the syntax for conditional and loop statements was similar but there were a few specialties in cpp. The first one being the compilation and the availability of different compilers.

C++ is a compiled language whereas python is interpreted. This means that we need special programs known as compilers to build the source code into an executable format and there are many different such programs (compilers) available in the market. To make things much worse, C++ compilers are not compatible with each other. The safest bet was to use ‘gcc’.

There are three ways in which we can compile a program:

1. ### **Direct command line** 

   > g++ a.cpp b.cpp and so on

   This way we have to specify the name of all source files one by one.

2. ### **Makefiles** 

   With make files, we write the whole command in blocks. We specify which compiler we would like to use, which libraries we would like to use and what are the source files. In a large project, we put the make files in individual folders and add one to the main folder. An example makefile could be found [here](https://github.com/iitrabhi/GSoC2019/blob/master/Scripts/cpp/learning-cpp/make/makefile). To understand more about makefiles you can go through the references.

3. ### **Cmake** 

   Cmake is a generator of makefiles. If we wish to execute our C++ program on various machines we would have to create makefiles specific to each one of them as they might have different compilers installed. This would be a cumbersome process without cmake. Cmake generates platform specific make files which makes it easy to compile our source code on different machines. An example makefile could be found [here](https://github.com/iitrabhi/GSoC2019/blob/master/Scripts/cpp/learning-cpp/cmake/CMakeLists.txt).

## References

- [CPP programming](https://www.programiz.com/cpp-programming)
- [Understanding compilation](https://www.toptal.com/c-plus-plus/c-plus-plus-understanding-compilation)
- [Understanding class definition](https://gist.github.com/darkstalker/eeb7e48a45f1b78db4a2c6ebfd01e926)
- [Understanding make](https://www.tutorialspoint.com/makefile/makefile_quick_guide.htm)
- [Understanding cmake 1](https://arne-mertz.de/2018/05/hello-cmake/)
- [Understanding cmake 2](http://derekmolloy.ie/hello-world-introductions-to-cmake/)
