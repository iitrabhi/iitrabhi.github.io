---
layout: post
title: "Stokes Flow Simplified: A Beginner's Guide to Ice Movement"
description: Dive into the basics of Stokes flow, a key concept in understanding how glaciers and ice sheets move.
categories:
  - simulation
tags:
  - fenics
---
  
As part of my postdoctoral work with NASA, I've started exploring a concept that's relatively new to me: Stokes flow, especially as it relates to ice movement. I aim to understand how massive ice formations like glaciers can flow and change shape over time. This is a vital piece of the puzzle in studying Earth's climate and environmental changes. In this article, I'll be sharing what I learn as I delve into the intriguing science behind ice's slow yet significant movement.

## What is Stokes's flow?
Stokes flow, also known as creeping flow or low Reynolds number flow, refers to fluid motion where viscous forces are dominant over inertial forces. This is typically the case when dealing with very slow fluid motion, small-length scales, or highly viscous fluids. The Reynolds number in these situations is much less than one.

In the context of ice, Stokes flow can be particularly relevant in glaciology, where it's used to model the movement of ice sheets and glaciers. Despite being a solid, ice flows under its own weight over `geological timescales.` This flow is often very slow, and the viscous forces within the ice are much more significant than inertial forces, making Stokes flow a suitable model for describing it.

> To draw parallels with solid mechanics, consider how a metal beam bends under weight. The beam deforms slowly, showing strain akin to how ice deforms under its weight and external forces. However, unlike the instantaneous response of the metal beam, ice movement is a much slower process, taking years to manifest noticeable changes. This difference in timescale is a key aspect of Stokes flow in glaciology.

## What is Stokes equation?
The Stokes equation, also known as the Stokes flow or the creeping flow equation, is a simplified form of the Navier-Stokes equation. It describes the **motion of a viscous fluid at very low Reynolds numbers**, where the inertial forces are negligible compared to viscous forces. 

The Reynolds number (Re) is a dimensionless parameter that characterizes the relative importance of inertia to viscosity in a fluid flow. **Low Reynolds Number (Re << 1):** the fluid behaves as if it has high viscosity relative to its density and velocity.
$$\mathcal{R}_e = \frac{\rho V L}{\mu}$$
where Re is the Reynolds number, $\rho$ is the density of the fluid, V is the velocity of the fluid, L is a characteristic length scale (such as the diameter of a pipe or the radius of a sphere), $\mu$ is the dynamic viscosity of the fluid.
$$\begin{split}- \nabla \cdot (\nabla u + p I) &= f \quad {\rm in} \ \Omega, \\                \nabla \cdot u &= 0 \quad {\rm in} \ \Omega. \\\end{split}$$
Where **u** is the velocity vector field of the fluid, **p** is the pressure field, **I** is the identity tensor, and **f** represents the external forces acting on the fluid.

---
These equations are the steady-state, incompressible Stokes equations, and they describe the behavior of a viscous, incompressible fluid within a domain $\Omega$. Let's break down what each equation is trying to convey:

1. The first equation is the momentum equation (also called the Stokes momentum equation) and can be written as:

$$\nabla \cdot (\nabla u + pI) = f$$
This equation describes the balance of forces within the fluid. Here's what it's trying to say:

- $\nabla u$ represents the gradient of the velocity vector field $u$. It describes how the velocity changes from point to point within the fluid.
- $pI$ represents the pressure field multiplied by the identity tensor $I$. This term represents pressure forces acting within the fluid.
- $f$ represents external forces applied to the fluid (e.g., gravitational forces or other forces).

The equation essentially states that the divergence of the sum of the velocity gradient ($\nabla u$) and the pressure  ($p$) is equal to the sum of external forces ($f$) acting within the fluid. This equation represents the balance of forces in the fluid. It is a simplified version of the Navier-Stokes momentum equation, appropriate for low Reynolds number flows where inertial forces are negligible.

2. The second equation is the continuity equation, and it can be written as:
$$\nabla \cdot u = 0$$

This equation enforces the condition of incompressibility within the fluid. Here's what it's trying to say:
- $\nabla \cdot u$ represents the divergence of the velocity field $u$. It quantifies how much the fluid expands or contracts at a given point.
- The equation asserts that the divergence of the velocity field $u$ must be zero everywhere within the fluid domain $\Omega$.

In other words, the continuity equation ensures that the fluid is incompressible, meaning that the volume of fluid elements remains constant as they move within the flow. It is a fundamental constraint for incompressible fluid flows.

In summary, the Stokes equations describe the behavior of a viscous, incompressible fluid by expressing the balance of forces (momentum equation) and ensuring incompressibility (continuity equation) within a given domain $\Omega$. They are particularly applicable to low Reynolds number flows, such as slow and viscous flows, and they are a simplified version of the more comprehensive Navier-Stokes equations.

---
## Where is it used?
The Stokes equations model fluid flow in situations where viscous forces dominate, and inertial forces are negligible. These equations are particularly applicable to low Reynolds number flows. Here are some practical examples of where the Stokes equations might be used:

1. **Microfluidics:** In microfluidic devices and systems, the dimensions are often on the micrometer scale, and flow velocities are relatively low. Examples include lab-on-a-chip devices, where small volumes of fluids are manipulated for applications like chemical analysis and medical diagnostics. The Stokes equations accurately describe and predict fluid behavior in microfluidic channels.
2. **Biological Fluid Dynamics:** The movement of fluids within biological systems, such as blood flow in capillaries or the movement of cilia and flagella in cells, often occurs at very low Reynolds numbers. The Stokes equations can be used to study and understand these biological fluid flows.
3. **Suspensions of Small Particles:** When small solid particles are suspended in a fluid, viscous forces strongly influence their motion. The Stokes equations are used to model the motion of these particles in applications like colloidal suspensions, where understanding particle behavior is important in fields like chemistry and materials science.
4. **Boundary Layer Flows:** In situations where a fluid flows near a solid boundary, such as the flow of a viscous fluid over a flat plate (Stokes boundary layer), the flow can be analyzed using the Stokes equations. This is important in aerodynamics and fluid mechanics for understanding boundary layer behavior and drag forces.
5. **Rheology of Viscous Materials:** When studying the behavior of highly viscous materials like specific polymers or lubricants, the Stokes equations may be used to model their flow characteristics. This is important in materials processing, manufacturing, and lubrication engineering industries.
6. **Lubrication Problems:** In applications involving lubrication, where a thin film of viscous fluid separates moving surfaces (e.g., in journal bearings or slider bearings), the Stokes equations can be used to model the lubricant flow and study the lubrication properties.

The full Navier-Stokes equations account for viscous and inertial effects in many practical situations, especially in everyday engineering and fluid dynamics. However, the Stokes equations serve as a valid approximation in situations where viscous forces dominate.

---
## Why is the Stokes equation not time-dependent?
The steady-state, incompressible Stokes equations are not time-dependent because they describe fluid flow in a steady-state condition, where the flow properties stay the same with time. 

There are several reasons why these equations are often used in situations where time-dependent effects are not a primary concern:

1. **Simplified Modeling:** In many practical scenarios, especially when dealing with highly viscous fluids like honey, the primary interest is understanding the fluid's steady-state behavior. Time-dependent effects, such as fluid transients or unsteady flows, may be minor or of secondary importance.
2. **Practical Applicability:** Certain applications, like the flow of honey in food processing or the behavior of lubricants in machinery, often operate in steady-state conditions. Modeling these systems as steady-state simplifies the analysis without significantly compromising accuracy.
3. **Mathematical Convenience:** Steady-state problems can be easier mathematically than time-dependent problems. The absence of time-dependent terms makes the equations more straightforward to handle and often leads to analytical or numerical solutions that are computationally less intensive.
4. **Physical Realism:** In some situations, it is a reasonable approximation to assume that the flow has reached a quasi-steady state quickly, and the time it takes for transient effects to dissipate is much shorter than the time scale of interest. This makes the steady-state assumption valid for practical purposes.

That said, the full Navier-Stokes equations, which include time-dependent terms, can be used when it's necessary to account for transient effects and unsteady flows or when analyzing systems with rapidly changing conditions. 

The choice between steady-state and time-dependent models depends on the specific problem, the time scales of interest, and the desired level of accuracy in the analysis. 

In some cases, transient effects may be significant and cannot be neglected, and in those cases, time-dependent versions of the equations would be more appropriate.

## Solving the Stokes problem
The Stokes equation leads to a Saddle Point problem in its discretized form. This means that the resulting matrix from the discretization is neither symmetric nor positive definite. Such matrices are challenging to solve efficiently due to their unique properties. In the context of ice flow, this challenge is compounded by the large scale and complexity of glacial systems.

- Prof. Gilbert Strang talks about this problem in great detail in his [lecture](https://math.mit.edu/classes/18.086/2006/am65.pdf).
- We also have a very good tutorial on the solution of the Stokes equation using FEM on the [FEniCS website.](https://fenicsproject.org/olddocs/dolfin/1.6.0/python/demo/documented/stokes-iterative/python/documentation.html)
- [This post](https://math.stackexchange.com/questions/397644/existence-and-uniqueness-of-stokes-flow) on stack overflow is also good for understanding the problem and the solution process.

Because of a [0] on the block diagonal, the conditioning of the system of equations is really bad. A special way of handling this using specially designed pre-conditioners is discussed in the FEniCS tutorial.

> We could use a direct solver such as MUMPS to get a solution.

