---
layout: post
title: "How to apply equality volume constraint with MMA for topology optimisation"
description: "This post has the trick get a constant volume with MMA"
categories: [simulation]
tags: [fem, mma, volume constraint, equality constraint, topology optimization]
typora-root-url: ../../../../website
---

Method of moving asymptotes (`MMA`) is one of the most popular optimization algorithms for topology optimization. The standard definition of  topology optimization is:

$$\begin{array}{ll}
\text { Minimize } & c(\mathbf{x})=\sum_{e=1}^{N_{e}}\left(x_{e}\right)^{p}\left\{\mathbf{u}_{e}\right\}^{\mathrm{T}}\left[\mathbf{k}_{e}\right]\left\{\mathbf{u}_{e}\right\} \\
\text { subject to } & \frac{V(\mathbf{x})}{f \times V_{0}}-1=0 \\
& {[\mathbf{K}]\{\mathbf{U}\}=\{\mathbf{F}\}} \\
& \mathbf{x}_{\min } \leq \mathbf{x} \leq \mathbf{x}_{\max }
\end{array}$$

whereas the standard definition of `MMA` is:

$$\begin{array}{ll}
\text { Minimize } & f_0(x) + a_0 \times z + \sum( c_i \times y_i + 0.5 \times d_i \times (y_i)^2 ) \\
\text { subject to } & f_i(x) - a_i \times z - y_i <= 0,  i = 1,...,m \\
& z >= 0,   y_i >= 0,         i = 1,...,m \\
&  x^{min}_j \leq x_j \leq x^{max}_j,    j = 1,...,N_e
\end{array}$$

## The problem

The algorithm of `MMA` is designed with inquality constraint ($$\leq$$) but the porblem of topology optimization is defined with an equality constriant. We want the volume of the design domain to achieve a target value ($$\text{vol frac}\times \text{initial vol}$$) and not go below than that. 

## The solution

What we want is:

$$\text{current volume} = \text{target volume} \implies \frac{\text{current volume} }{ \text{target volume}}-1=0$$

To implement the above equality constraint with `MMA` we have to apply two constraints:

$$+\frac{\text{current volume} }{ \text{target volume}+0.01}-1 \leq 0 \implies \text{current volume} \leq \text{target volume}+0.01  $$

$$-\frac{\text{current volume} }{ \text{target volume} - 0.01}+1 \leq 0 \implies \text{target volume} - 0.01 \leq \text{current volume}$$

With the above two constraints applied we can force `MMA` to give us the design volume equal to the target volume.

$$\text{target volume} - 0.01 \leq \text{current volume}\leq \text{target volume}+0.01$$

# References

- [Github repo with implementation](https://github.com/iitrabhi/topopt2d_MMA)

