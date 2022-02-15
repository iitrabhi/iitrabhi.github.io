---
layout: post
title: "Adaptive isogeometric topology optimization"
description: "Adaptive refinement for TO using PHT splines."
langs: [MATLAB, python]
year: "January 2022"
typora-root-url: ../../../../website
---

## Highlights

- A novel adaptive isogeometric topology optimization methodology using PHT-splines is presented.
- Complex geometries could be optimized with the help of multi-patch NURBS discretization.
- Variation of the element density is used to mark the sub-domains for refinement.
- First neighbourhood filter is utilized to achieve mesh independence.
- The proposed approach achieves 30 − 60% reduction in degree’s-of-freedom and 80 − 90% reduction in CPU time.

---

![image-20220215145840919](/assets/images/image-20220215145840919.png)

<figcaption>Hierarchical local refinement of the knot mesh and the resulting generation of PHT-spline basis. (a) PHT-basis on the initial level T1 mesh (b) Refined basis on the level T2 mesh.</figcaption>

---

