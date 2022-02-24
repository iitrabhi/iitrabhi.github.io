---
layout: post
title: "Adaptive phase field fracture"
description: "Computational codes to carry out adaptive PFF in FEniCS."
langs: [python]
year: "April 2019"
date: 2019-04-01
typora-root-url: ../../../../website
---

Phase-field fracture(PFF) analysis was not part of my Ph.D. proposal. Still, since there existed a lot of literature on the phase-field method in the domain of fracture, I started learning it to use it further in topology optimization. I got a lot of help from my friend [Dr. Tushar K Mandal](tusharmandal.com), who was actively working in the phase-field method at the time. Also, since I was working in FEniCS and there existed an open-source code of PFF in FEniCS by [Dr. Hirshikesh](https://scholar.google.co.in/citations?user=2-tIkOcAAAAJ&hl=en) and [Dr. Emilio Martínez Pañeda](https://scholar.google.co.in/citations?user=DDVhQIcAAAAJ&hl=en), I opted to start learning about phase-field with PFF. One of my colleagues, [Meenu Krishnan](https://www.researchgate.net/profile/Meenu-Krishnan), also joined me in this effort. All of this led to the development of a generalized and modularized code of PFF in FEniCS with support for parallelization. The code can handle all the three most popular PFF models and supports checkpointing and restart analysis. We also developed time adaptivity and mesh adaptivity algorithms for PFF. This code is one of the significant achievements of my Ph.D.

![Ph.D.@2x](/assets/images/Ph.D.@2x.png)

# Outcome 1

[An auto-adaptive sub-stepping algorithm for phase-field modeling of brittle fracture](https://www.sciencedirect.com/science/article/pii/S0167844220301981)

In this work, we propose an auto-adaptive displacement stepping algorithm to auto-adaptively adjust the displacement step size while solving the quasi-static brittle fracture propagation problem with the phase-field method. There were multiple works on adaptive meshing for phase field fracture, but none on adaptive time stepping. Through this project we achieved adaptive time stepping in phase field fracture.

![image-20220217204215393](/assets/images/image-20220217204215393.png)

<figcaption>(a) Approximation of crack topology in PFM framework using damage parameter (*d*) and length scale parameter (l ); (b) distribution of the damage parameter over the domain.</figcaption>

![image-20220217204643579](/assets/images/image-20220217204643579.png)

<figcaption>Step-size change ratio as a function of error excess.</figcaption>

![image-20220217204754198](/assets/images/image-20220217204754198.png)

<figcaption>Proposed sub-stepping algorithm decides the step size based on slope of the energy curve.</figcaption>

![image-20220217204829797](/assets/images/image-20220217204829797.png)

<figcaption>Comparison of standard alternate minimization algorithm (AM) and auto-adaptive sub-stepping algorithm (AS) for the mode-1 problem: (a) reaction force versus displacement graph, (b) CPU time versus displacement graph. The proposed sub-stepping algorithm is able to accurately capture the response of the system with 90% reduction in CPU time to reach peak load and 62% reduction in CPU time for complete analysis.</figcaption>

![image-20220217204922697](/assets/images/image-20220217204922697.png)

![image-20220217210324812](/assets/images/image-20220217210324812.png)

![image-20220217210423964](/assets/images/image-20220217210423964.png)

## Conclusion

In this work, we proposed an auto-adaptive displacement stepping algorithm with a sub-stepping algorithm to achieve adaptivity in the displacement step for solving the coupled system of non-linear equations arising from the phase field formulation of brittle fracture. We achieved a reduction of 78–90% in the CPU time required to reach the peak reaction force. Also, there is a reduction of 56–78% in the CPU time corresponding to the final solution.

# Outcome 2

[A length scale insensitive phase field model for brittle fracture of hyperelastic solids](https://www.sciencedirect.com/science/article/pii/S0013794420307797)
