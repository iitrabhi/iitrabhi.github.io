---
layout: post
title: "Matrix multiplication is associative, but you have to be careful."
description: "The order of operation of multiplying three or more matrices could significantly affect the run-time of the algorithm."
categories:  [coding]
tags: [hpc, matrix, python]
typora-root-url: ../../../../website
---
Matrix multiplication is one of the most fundamental operations in the finite element method. Recently there has been a lot of talk about  [AI that came up with a faster way to multiply matrices](https://www.deepmind.com/blog/discovering-novel-algorithms-with-alphatensor). And this is a big deal for the finite element world too. Faster matrix multiplication could potentially lead to faster run times for existing algorithms on existing hardware. In this post, I talk about the associative property of matrices and why it is important to have the correct order of operation for multiplying three or more matrices.

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

| size     | operations | time (seconds) |
| -------- | ---------- | -------------- |
| 10 | 1.90E+03  | 3.80E-07       |
| 1E+2  | 1.99E+06   | 3.98E-04       |
| 1E+4 | 2.00E+12   | 4.00E+02       |
| 1E+6 | 2.00E+18   | 4.00E+08       |
| 1E+9 | 2.00E+24   | 4.00E+14 $\approx$ 13 million years |

As you can see from the above table, the cubic order increase in operations leads to a huge increase in time-to-solution.  


## What is the associative property of multiplication?

It is a mathematical rule that states that the order in which matrices are grouped in a multiplication problem does not change the product, i.e. if we have three matrices with the correct size, then,

$$(AB)C=A(BC)$$




As is true for most of the concepts we learn in high school, every rule comes with a set of assumptions and limitations. In the realm of small scale, the multiplication order would not create much difference, but once we start working with huge matrices, we have to give special consideration to the minutest things.


## References
- [First extension of AlphaZero to mathematics unlocks new possibilities for research](https://www.deepmind.com/blog/discovering-novel-algorithms-with-alphatensor)
- [Big O notation](https://www.freecodecamp.org/news/big-o-notation-why-it-matters-and-why-it-doesnt-1674cfa8a23c/)