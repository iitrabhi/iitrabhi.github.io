---
layout: post
title: "Some surprising findings on direct vs iterative in MATLAB."
categories: 
  - [coding]
typora-root-url: ../../../website
---

 I am trying to solve a million DoF system with PETSc as well as MATLAB. In both cases the system has been assembled and I have the K matrix and the F vector.

- MATLAB - takes 4 seconds to solve the system using sparsity and decomposition.
- PETSc is unable to handle the system with direct solvers and iterative solvers take 60 seconds to solve the system with 17 iterations. 

This has baffled me, as till now I was under the impression that I will always get better speeds with PETSc. From matlab documentation:

- The speed of solving a linear system with a direct method strongly depends on the density and fill pattern of the coefficient matrix.
- The speed of solving a linear system with an indirect method does not depend as strongly on the fill pattern of the coefficient matrix as a direct method. **However, using an iterative method typically requires tuning parameters for each specific problem.**

Previously I have read in the large scale paper on FEniCS that the number of iterations required to reach solution increases with increasing the number of DoF in the system and keeping everything same. In ideal case it should be the same. 

From my own testing also I have found out that the residual norm increases with incresing numbers of DoF and thus requires more iterations to converge. 

Also, it should be noted that the problem should be well defined. The loading on the system should not be exorbitantly large.

My system was not converging with the standard choice of CG + ILU so based on [this post](https://scicomp.stackexchange.com/questions/2369/what-is-a-robust-iterative-solver-for-large-3-d-linear-elastic-problems?noredirect=1&lq=1) I changed it to CG + AMG and it converged in 20 iteration. 

> a direct solver is equivalent to an iterative method that stops after one preconditioner application. [from here](https://pages.tacc.utexas.edu/~eijkhout/pcse/html/petsc-solver.html#:~:text=32.2%20Direct%20solvers&text=PETSc%20has%20some%20support%20for,stops%20after%20one%20preconditioner%20application.)

---

## Summary

- The problem with direct solver is that they become prohibitive for large problems because of memory requirement.
- For that they are using sparse matrices with integer indexing to reduce the memory expense
- In case of PETSc, it cannot handle large problems with direct solvers(ndof>50,000) and in the documentation itself suggests to go for iterative solvers.
- Iterative solvers need to iterate to reduce the residual.
- The residual is dependent on initial guess as well as the size of the system
- For the same differential equation, iterative solvers require more iterations to converge with increasing number of DoF's.
- Moreover they are dependent on the pre-conditioners. For a simple cantilever beam problem with 1 Million elements no pre-conditioner other than AMG was producing results in less than 100 iterations. AMG took 16 iterations to converge to a solution.
- For a large system the condition number is bad for the same differential equation, thus we need pre-conditioner to condition the system. The amount of iteration that are required to solve the system depends on the solver and pre-conditioner used. This combination could result in better performance as compared to direct solvers for very large problem which are prohibitively expensive in terms of memory requirement.
- For symmetric positive definite systems that we can decompose and store in the memory, direct solvers will perform better. (There are many other factors that go into this but for our case this is the only thing.)

---

## References

- [Direct v/s iterative](https://caendkoelsch.wordpress.com/2018/11/29/direct-vs-iterative-solvers-in-fem/#:~:text=Direct%20Solver%3A,for%20computationally%20less%20expensive%20problems.)
- [MATLAB direct vs iterative](https://in.mathworks.com/help/matlab/math/iterative-methods-for-linear-systems.html#:~:text=MATLAB%20implements%20direct%20methods%20through,a%20finite%20number%20of%20steps.)
- [Nico's query on a similar problem to mine](https://scicomp.stackexchange.com/questions/5600/best-choice-of-solver-for-a-large-sparse-symmetric-but-not-positive-definite-s)
- [Which solver to choose](https://scicomp.stackexchange.com/questions/81/what-guidelines-should-i-follow-when-choosing-a-sparse-linear-system-solver)
- [Matlab vs python sparse solvers ](https://stackoverflow.com/questions/64401503/is-there-a-way-to-further-improve-sparse-solution-times-using-python)

