---
layout: post
title: "Adaptive isogeometric topology optimization"
description: "Adaptive refinement for TO using PHT splines."
langs: [MATLAB, python]
year: "January 2022"
date: 2022-01-01
typora-root-url: ../../../../website
---

The goal of my research during the Ph.D. was to solve the problem of topology optimization using isogeometric analysis. This project completed the goal, and finally, I achieved topology optimization using adaptive refinement based on PHT-splines.

## Highlights

- A novel adaptive isogeometric topology optimization methodology using PHT-splines is presented.
- Complex geometries could be optimized with the help of multi-patch NURBS discretization.
- Variation of the element density is used to mark the sub-domains for refinement.
- First neighbourhood filter is utilized to achieve mesh independence.
- The proposed approach achieves 30 − 60% reduction in degree’s-of-freedom and 80 − 90% reduction in CPU time.

---

![image-20220215155746693](/assets/images/image-20220215155746693.png)

<figcaption>Hierarchical local refinement of the knot mesh and the resulting generation of PHT-spline basis. (a) PHT-basis on the initial level T1 mesh (b) Refined basis on the level T2 mesh.</figcaption>

![image-20220215155701639](/assets/images/image-20220215155701639.png)

<figcaption>The GIFT framework. (a) The initial coarse mesh is discretized with the help of NURBS basis function. (b) The parametric mesh for both analysis and geometric description is based on the initial coarse mesh. Only the analysis mesh is refined with PHT-basis to properly approximate the field variable. (c) The field variable is approximated over the initial coarse geometry but with refined PHT-basis.</figcaption>

![image-20220215155957651](/assets/images/image-20220215155957651.png)

<figcaption>Illustration of Adaptive meshing strategy for 2D PHT-Spline elements. The elements along the boundary are marked and refined to achieve a smooth boundary in the optimized solution.</figcaption>

![image-20220215160111698](/assets/images/image-20220215160111698.png)

<figcaption>Topology of the final optimized design of cantilever problem using APHT-ITO, super- imposed over the adaptively refined mesh.</figcaption>

![image-20220215160156412](/assets/images/image-20220215160156412.png)

<figcaption>Topology of final optimized design of the L-shaped geometry using APHT-ITO super- imposed over the adaptively refined mesh.</figcaption>

![image-20220215160519441](/assets/images/image-20220215160519441.png)

<figcaption>Topology of final optimized design of the plate-with-a-hole geometry using APHT-ITO superimposed over the adaptively refined mesh.</figcaption>

![image-20220215160601608](/assets/images/image-20220215160601608.png)

<figcaption>Topology of final optimized design of the 3-D cantilever problem, after having applied three levels of adaptive refinement, using the APHT-ITO algorithm.</figcaption>

## Conclusion

This work presents a novel adaptive refinement framework for topology optimization (TO) using isogeo- metric analysis (IGA). Within the context of IGA, the recently developed Geometry Independent Field approximaTion (GIFT) framework for Polynomial splines over Hierarchical T-meshes (PHT)-splines is utilized. The geometry is modeled with multi-patch Non-Uniform Rational B-Splines (NURBS), and PHT-splines are utilized to model the solution field in an adaptive setting. The adaptivity algorithm is controlled by the iterative variation of the density variable over each element. First neighbourhood filter is utilized to achieve mesh independence. The benefit of the hierarchical tree structure of PHT splines is made to store and access the variables of the TO problem.

Numerical studies conducted on four different design domains validated the effectiveness of the proposed method in solving 2-D and 3-D compliance minimization problems of TO. Compared to non-adaptive isogeometric topology optimization using a global refinement strategy of the same level, a reduction in degree’s-of-freedom between 30−60% was achieved. In addition, a considerable improvement in run-time in the range of 80 − 90% was also observed; this will further scale when the complexity, refinement levels, and iteration steps for the problem in consideration are increased. Therefore the proposed method is com- putationally efficient as compared to conventional ITO methods. Further, since the methodology is based on the GIFT framework, it can directly import multi-patch geometries from any professional computer- aided design (CAD) package, which helps establish a seamless CAD to engineering workflow.
