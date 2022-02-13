---
layout: post
title: "Non-linear analysis of RCC structures for seismic loading"
description: "Understanding Newton Raphson and Arc Length method using ABAQUS."
langs: [ABAQUS, python]
year: "March 2012"
typora-root-url: ../../../../website
---

##  Objective of the study

- Studying different material models available in ABAQUS v.6.9.
- Modelling structural elements and seismic condition in ABAQUS v6.9 for conducting the non-linear analysis. 
- Studying the effect of material model parameters and mesh on analysis results.
- Conducting a comparative study of the material models of concrete and proposing the best suited model of concrete for the purpose of seismic analysis.

## Effect of mesh density on displacement field

![image-20220213213900583](/assets/images/image-20220213213900583.png)

| Displacement          | Ele_size = 300mm | 200mm    | 100mm   |
| --------------------- | ---------------- | -------- | ------- |
| Linear Hexahedral     | 3.509mm          | 3.501mm  | 3.069mm |
| Quadratic Hexahedral  | 2.761mm          | 3.051 mm | 3.058mm |
| Linear Tetrahedral    | 1.855mm          | 2.918mm  | 2.987mm |
| Quadratic Tetrahedral | 2.798mm          | 3.065mm  | 3.086mm |

Now, let's make an important observation: By meshing a part with a certain type of element and using a certain size and shape, we impose additional constraints on the part. The meshed part must conform to applied loads and restraints. But being meshed, it must also conform to constraints imposed by meshing. In other words, deformation must satisfy loads and restraints and be piecewise linear. Because the meshed part has the additional constraints, the process of meshing makes it stiffer.

The amount of additional stiffness depends on the element and their size. First-order elements require that the displacement field be piecewise linear. This is more restrictive than in second-order elements where the displacement field must be piecewise parabolic. Larger elements add more stiffness than small ones. However, the effect of added stiffness (call it artificial stiffness) always accompanies finite-element models. The effect of artificial stiffness is small but demonstrable in most cases, even with a reasonably well-refined mesh.
