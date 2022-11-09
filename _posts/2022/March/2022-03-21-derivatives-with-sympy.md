---
layout: post
title: "How to take derivatives with sympy."
description: "Generate or validate your derivatives with symbolic calculation in python."
categories: 
  - [coding]
typora-root-url: ../../../../website
---

## The problem

Sometimes we need some validation for the derivatives of a function we are taking. I had a discussion regarding this with Tushar, and he suggested using the `sympy` library to verify the derivatives. We can not rely entirely on these symbolic libraries for the accuracy of the solution as they are still under active development. Still, I think they can give us some confidence in the solution.

## The solution

Define the derivative in symbolic language and use the `diff` function of sympy. In this example, my function to derivate is the p-norm of another function.

$$I = \left(\int{f(x)^p dx}\right)^{1/p}$$

```python
from sympy import *
from IPython.display import display, Latex
x,p = symbols("x p");
f = symbols("f", cls=Function);
I = integrate(f(x)**p,x)**(1/p)
sol = simplify(I.diff(f(x)));
print(sol)
print(latex(sol))
```
**Output** 

```mk
Integral(f(x)**p, x)**((1 - p)/p)*Integral(f(x)**(p - 1), x)
\left(\left(\int f^{p}{\left(x \right)}\, dx\right)^{\frac{1 - p}{p}}\right) \int f^{p - 1}{\left(x \right)}\, dx
```

We can then see this integral in latex by using the command.

```python
display(Latex('$'+latex(sol)+'$'))
```

**Output**

$$\left(\left(\int f^{p}{\left(x \right)}\, dx\right)^{\frac{1 - p}{p}}\right) \int f^{p - 1}{\left(x \right)}\, dx$$

