---
layout: post
title: "How to add images, latex and code to blog."
categories: misc
typora-root-url: ../../website
---
On a daily basis there are only three things that I have to put in my notes. My understanding of the research in well formatted text which I can easily write in markdown. The next important thing is to add latex based formulas. This is by default supported in markdown as well as Jekyll. Thus we can directly write what ever formula we wish to write inside double dollars and that will get parsed as latex formula:

$$R_{i, p}(\xi)=\frac{N_{i, p}(\xi) w_{i}}{W(\xi)}, \quad \text { with } \quad W(\xi)=\sum_{k=1}^{n} N_{k, p}(\xi) w_{k}$$

The next big thing is adding images to markdown. Now this is not supported by defualt in markdown since it is a text document, but by using [Typora](https://typora.io) we can directly paste our images to the markdown file. It will automatically add the image to the directory and paste a relative link to that file. 

For this method to work on typora you have to go to preferences → images and then change the copy action to `copy image to custom folder`

![image-20201122123717816](/assets/images/image-20201122123717816.png)

After that whenever you paste a file it will automatically save to the desired folder and referenced in your document. This way there is no need to link your image file manually in markdown. For example I can take a screenshot and directly paste the file here

![image-20201122115600774](/assets/images/image-20201122115600774.png)

In a similar fashion I can add a properly color-coded code block by putting it in 

> \`\`\`python ----- \`\`\`

```python
from dolfin import *
L = 25.0
H = 1.0
Nx = 250*6
Ny = 1000
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), Nx, Ny, "crossed")
def eps(v):
    return sym(grad(v))

E = Constant(1e5)
nu = Constant(0.3)

mu = E / 2 / (1 + nu)
lmbda = E * nu / (1 + nu) / (1 - 2 * nu)

def sigma(v):
    return lmbda * tr(eps(v)) * Identity(2) + 2.0 * mu * eps(v)

rho_g = 1e-3
f = Constant((0, -rho_g))

V = VectorFunctionSpace(mesh, "Lagrange", degree=1)
du = TrialFunction(V)
u_ = TestFunction(V)
a = inner(sigma(du), eps(u_)) * dx
l = inner(f, u_) * dx
k = assemble(a)
```

