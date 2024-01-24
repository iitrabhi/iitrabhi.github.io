---
layout: post
title: First steps to understanding ice flow simulation with IcePack
description: What is Stokes equation? How to simulate Stokes equation? What is ice flow?
categories:
  - simulation
tags:
  - fenics
---
I know nothing about ice simulations, but I know a lot about variational mechanics and FEniCS. This post will record my learnings as I try to understand the subject with the Firedrake/FEniCS-based ice simulation package IcePack.

---

In 1848, J.D. Forbes was the first Western scientist to identify viscous deformation as the reason glaciers flow correctly. In his original paper on the subject, Forbes has a wonderful quote about this realization:

> There is something pleasing to the imagination in the unexpected analogies presented by a torrent of fiery lava and the icy stream of a glacier.

Both lava and ice flow can be described by the same mathematics, and that mathematics is the _Stokes equations._

---

The Stokes equations have three parts: a conservation law, a constitutive relation, and boundary conditions.

### Conservation law
$$\nabla\cdot\tau - \nabla p + f = 0.$$
### Constitutive relation
$$\dot\varepsilon(u) \equiv \frac{1}{2}\left(\nabla u + \nabla u^\top\right)$$
The deviatoric stress and strain rate tensor are linearly proportional for a plain old, viscous Newtonian fluid. But glacier ice is not a Newtonian fluid!
$$\dot\varepsilon = A|\tau|^{n - 1}\tau$$
### Boundary conditions
At the ice surface, there is effectively zero stress:
$$(\tau - pI)\cdot\nu|_{z = s} = -p_0\nu$$
Things get much more interesting at the ice base because there are different boundary conditions in the normal and tangential directions. In the normal direction, the ice velocity has to equal the rate of basal melting:
$$u\cdot\nu|_{z = b} = \dot m.$$
In the tangential direction, frictional contact with the bed creates resistive stresses. The **sliding law** contains the relationship between resistive stresses and the ice velocity and other fields. Weertman sliding is a power-law relation between stress and sliding speed:
$$(\tau - pI)\cdot\nu|_{z = b} = -C|u|^{1/m - 1}u,$$
## Action principle
Instead of dealing with the full complexity of the Stokes equations, this principle simplifies the problem into finding the critical point of a functional called the action.
$$J = \int_\Omega\left(\frac{n}{n + 1}A^{-1/n}|\dot\varepsilon|^{1/n + 1} - p\nabla\cdot u - f\cdot u\right)dx + \frac{m}{m + 1}\int_{\Gamma_b}C|u|^{1/m + 1}ds.$$
This method is advantageous because it allows the use of more effective numerical methods for solving optimization problems. The action itself is a mathematical expression that incorporates various aspects of ice movement and behavior, including viscosity, pressure, and sliding at the base of the ice. This approach is not only more efficient computationally but also more concise in its mathematical formulation.

