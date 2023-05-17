---
layout: post
title: "How I sped up my finite element code 100x by moving parentheses."
description: "The order of operation of multiplying three or more matrices could significantly affect the run-time of the algorithm."
categories:  [coding]
tags: [hpc, matrix, python, numpy, multiplication,home]
typora-root-url: ../../../../website
---
Matrix multiplication is one of the most fundamental operations in the finite element method. Recently there has been a lot of talk about  [AI that came up with a faster way to multiply matrices](https://www.deepmind.com/blog/discovering-novel-algorithms-with-alphatensor). And this is a big deal for the finite element world too. Faster matrix multiplication could lead to shorter run times for existing algorithms on existing hardware. In this post, I want to talk to you about one of my recent discoveries during a debugging session.

**But first, let me define the problem.**

We made a code to perform an eigen analysis of a finite element model. The first step for such analysis is to find the eigenvalues and eigenvectors of the system. The next step is to find the modal mass and the mass participation factor of the system, which is given by
 $$M_{eff}=\left(\frac{[\phi]^{T}[M]\{T\}}{\sqrt{m}}\right)^{2}$$

In the above equation, $\phi$ is the eigenvector matrix with a size of $(n\times m)$ (where $n$ is the number of degrees of freedom of the structure), $M$ is the mass matrix with a size of $(n\times n)$ , and $T$ is the rigid body vector with a size of $(n\times 1)$. There are two ways in which we can perform the multiplication in the numerator

1. $([\phi]^{T}[M])\{T\}$
2. $[\phi]^{T}([M]\{T\})$

For small-scale systems, these options will have the same time-to-solution. But, for large-scale systems, we get a massive difference in the time-to-solution between the two approaches. Let us boil down the problem to its basics and understand what is happening here.

## How many operations can a computer processor perform per second?

The speed of the computer processor is usually marketed in gigahertz (GHz). A standard Intel-i7 processor can achieve speeds up to 5GHz. This means that the processor can perform $5 \times 10^9$ cycles per second [(Ref)](https://stackoverflow.com/questions/43651954/what-is-a-clock-cycle-and-clock-speed). There is a whole lot of theoretical foundation that is required to understand what the computer can do in a single cycle. Still, for the sake of simplicity, we will assume that the computer can perform a single addition or multiplication in a single cycle. [(Ref)](https://qr.ae/pvQsE4)

Thus, theoretically, with our gross simplification, it should take a computer equipped with an  Intel-i7 processor one second to perform $5 \times 10^9$ operations. 

## How many operations are there in schoolbook matrix multiplication?
Multiplying two vectors of $n\times 1$ size takes  $n$ multiplications and $n-1$ additions.

$$
\left[\begin{array}{ll}
a & b 
\end{array}\right]\left[\begin{array}{l}
e \\
g 
\end{array}\right]=ae+bg
$$

Thus the total number of operations in a vector—vector dot product is : 

$$n+(n-1) = 2n-1$$

The resulting matrix from a $n \times n$ matrix-matrix dot product has $n^2$ elements. 

$$\begin{equation}
\left[\begin{array}{ll}
a & b \\
c & d
\end{array}\right] \times \left[\begin{array}{c c}
e & f \\
g & h
\end{array}\right]=\left[\begin{array}{ll}
a e+b g & a f+b h \\
c e+d g & c f+d h
\end{array}\right]
\end{equation}$$ 

Thus, the total number of operations in the matrix-matrix dot product is [(Ref)](https://math.stackexchange.com/questions/484661/calculating-the-number-of-operations-in-matrix-multiplication#:~:text=Thus%20the%20total%20number%20of,%3DO(n3).)

$$n^2 \times (2n-1) = 2n^3-n^2 = O(n^3).$$

Similarly, the total number of operations in a matrix-vector dot product is
$$n \times (2n-1) = 2n^2-n = O(n^2).$$

| size     | operations | time (seconds) |
| -------- | ---------- | -------------- |
| 10 | 1.90E+03  | 3.80E-07       |
| 1E+2  | 1.99E+06   | 3.98E-04       |
| 1E+4 | 2.00E+12   | 4.00E+02       |
| 1E+6 | 2.00E+18   | 4.00E+08 $\approx$ 13 years      |
| 1E+9 | 2.00E+24   | 4.00E+14 $\approx$ 13 million years |

As you can see from the above table, the cubic order increase ([Ref.](https://www.freecodecamp.org/news/big-o-notation-why-it-matters-and-why-it-doesnt-1674cfa8a23c/)) in operations leads to a huge increase in time-to-solution.  


## What is the associative property of multiplication?

It is a mathematical rule that states that the order in which matrices are grouped in a multiplication problem does not change the product, i.e., if we have three matrices with the correct size, then,

$$(AB)C=A(BC)$$

As is true for most of the concepts we learn in high school, every rule comes with a set of assumptions and limitations. In the realm of small scale, the multiplication order would not create much difference, but once we start working with huge matrices, we have to give special consideration to the simplest of operations.

## The finding
In one of the previous posts, I described [sparse matrices and why every researcher should understand and use them](https://abhigupta.io/2022/01/31/please-use-sparse-matrices.html). That handled one of the big problems of RAM requirement. The next step is to speed up the computation by developing the code based on good programming practices. I compiled most of my findings in this presentation: [Run-time from 300 years to 300 min: Lessons learned in large-scale modeling in FEniCS](https://www.researchgate.net/publication/352643174_Run-time_from_300_years_to_300_min_Lessons_learned_in_large-scale_modeling_in_FEniCS). By properly profiling the code for bottlenecks, we can figure out ways to increase its speed.

Let's create a test problem to understand

```python
from scipy.sparse import csr_matrix
import numpy as np

n=100000
m=40

t = np.ones(n)
rows_n=np.arange(n)
cols_n=np.arange(n)
val_n = np.ones(n)
mat_n = csr_matrix((val_n, (rows_n, cols_n)), shape=(n, n), dtype=int)
mat_m=csr_matrix(np.random.randint(10, size=(m, n)))/100
```

| Option 1 --- $([\phi]^{T}[M])\{T\}$ | Option 2 --- $[\phi]^{T}([M]\{T\})$ |
| -------- | -------- |
|Here the multiplication order is<br> $[m \times n]\cdot[n\times n] = [m\times n]$<br>$[m \times n]\cdot[n\times 1] = [m\times 1]$          |     Here the multiplication order is<br> $[n \times n]\cdot[n\times1] = [n\times 1]$<br>$[m \times n]\cdot[n\times1] = [m\times 1]$     |
|`(mat_m.dot(mat_n)).dot(t)` | `mat_m.dot(mat_n.dot(t))`|
|Wall time: 57.5 ms | Wall time: 6.77 ms |

The simple example presented above shows that there is almost a $7\times$ speed boost by just changing the multiplication order. This increase is because, in the first case, we multiply two matrices resulting in a matrix. Then we multiply the matrix with a vector. In the second case, we multiply the matrix to a vector, resulting in a vector, and then multiply the resulting vector to a matrix. Since matrix-matrix multiplication scales at $O(n^3)$ and matrix-vector product scales at $O(n^2)$, we get substantial savings in _Option-2_.

In the actual problem, I achieved a speed boost of around $100\times$. The run time for a single iteration reduces from 2 seconds to around 0.02 seconds. The simulation that took 2-3 hours to complete is completed in around 2mins. 🥳

Think of the kind of analysis we can perform if our code can handle a billion degrees of freedom.
