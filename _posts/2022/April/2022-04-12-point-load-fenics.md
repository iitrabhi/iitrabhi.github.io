---
layout: post
title: "How to apply point load in FEniCS."
description: "Applying point load in FEniCS is not straightforward but doable."
tag: 
  - [coding, fenics]
typora-root-url: ../../../../website
---

## The problem

I am trying to apply point load in FEniCS. It's a load term and thus couldnot be applied as `DirichletBC` and `pointwise` in FEniCS. The `PointSource` method also limits us to application of a single term and we need to hack our way out from vector terms. Also, from what I understand the `PointSource` is applied on assembled system. But, I would like to apply the load in the variational form.

## The solution

We can use `UserExpression` to define an expression that would be zero everywhere on the mesh and equal the the load at the point of application. This could be achieved with the following code snippet

```python
class PointLoad(UserExpression):
    def __init__(self, pt, vl,tol,**kwargs):
        super().__init__(**kwargs)
        self.point = pt
        self.value = vl
        self.tol = tol
    def eval(self, values, x):
        if near (x[0], self.point[0],self.tol) and near(x[1], self.point[1],self.tol):
            values[0] = self.value[0]
            values[1] = self.value[1]
        else:
            values[0] = 0
            values[1] = 0
    def value_shape(self):
        return (2,)
```

The class `PointLoad` is inherited from the `UserExpression` class. It is important to call the `super().__init__(**kwargs)` method to inherit all the functionality of the parent class. The `eval` method contains the logic for the evaluation of the expression over the mesh. It takes in the co-ordinate vector `x` and the values vector `values`.The expression is evaluated on each vertex of the mesh, and thus, if we specify the load vertices we can apply the point load as follows

```py
f0 = PointLoad(pt=(2.5,0), vl=(a0,b0), tol=1e-1,degree = 1)
l = dot(f0,v)*dx
```

## References

- [Point load question on discourse](https://fenicsproject.discourse.group/t/point-load-and-body-force/2556/7)
- [Point load on discourse with solution](https://fenicsproject.discourse.group/t/how-to-address-point-boundary-condition/3286/2)

