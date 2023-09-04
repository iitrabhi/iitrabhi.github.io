---
layout: post
title: Achieving significant computational gains with sparse matrices and proper multiplication order.
description: Changing the order of matrix multiplication can significantly impact run time. In this post, I discuss how I achieved a speed boost of around 100 times.
categories:
  - coding
tags:
  - sparse
  - matrix
  - hpc
---

## The reason

 In one of the previous posts, I described [sparse matrices and why every researcher should understand and use them](_posts/2022/Jan/2022-01-31-please-use-sparse-matrices.md). That handled one of the big problems of RAM requirement. The next step is to speed up the computation by developing the code based on good programming practices. I compiled most of my findings in this presentation: [Run-time from 300 years to 300 min: Lessons learned in large-scale modeling in FEniCS](https://www.researchgate.net/publication/352643174_Run-time_from_300_years_to_300_min_Lessons_learned_in_large-scale_modeling_in_FEniCS). By properly profiling the code for bottlenecks, we can figure out ways to increase its speed. In this post, I would like to share one such finding.

## The finding

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

Here the multiplication order is

1. $$[n \times n]\cdot[n\times1] = [n\times 1]$$
2. $$[m \times n]\cdot[n\times1] = [m\times 1]$$ 


```python
%%time
mat_m.dot(mat_n.dot(t))
```

    CPU times: user 6.53 ms, sys: 1.54 ms, total: 8.07 ms
    Wall time: 6.77 ms
    
    array([4495.66, 4500.96, 4515.95, 4497.58, 4495.07, 4500.8 , 4515.05,
           4511.28, 4487.11, 4492.63, 4504.96, 4502.63, 4501.69, 4505.42,
           4494.06, 4504.79, 4517.84, 4501.67, 4507.63, 4495.58, 4497.92,
           4495.33, 4495.05, 4506.82, 4509.54, 4510.65, 4497.09, 4506.66,
           4496.97, 4503.77, 4489.17, 4499.72, 4501.36, 4499.38, 4485.01,
           4490.23, 4502.94, 4505.94, 4503.82, 4499.95])

Here the multiplication order is

1. $$[m \times n]\cdot[n\times n] = [m\times n]$$ 
2. $$[m \times n]\cdot[n\times 1] = [m\times 1]$$ 


```python
%%time
(mat_m.dot(mat_n)).dot(t)
```

    CPU times: user 44.9 ms, sys: 5.98 ms, total: 50.9 ms
    Wall time: 57.5 ms
    
    array([4495.66, 4500.96, 4515.95, 4497.58, 4495.07, 4500.8 , 4515.05,
           4511.28, 4487.11, 4492.63, 4504.96, 4502.63, 4501.69, 4505.42,
           4494.06, 4504.79, 4517.84, 4501.67, 4507.63, 4495.58, 4497.92,
           4495.33, 4495.05, 4506.82, 4509.54, 4510.65, 4497.09, 4506.66,
           4496.97, 4503.77, 4489.17, 4499.72, 4501.36, 4499.38, 4485.01,
           4490.23, 4502.94, 4505.94, 4503.82, 4499.95])

From the simple example presented above, we can see that there is almost a $$7\times$$ speed boost by just changing the multiplication order. In the actual problem, I achieved a speed boost of around $$100\times$$. The run time for a single iteration reduces from 2 seconds to around 0.02 seconds. The simulation that was taking 2-3 hours to complete now gets completed in around 2mins. 🥳



