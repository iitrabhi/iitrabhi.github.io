---
layout: post
title: "How to save a surface plot generated with matplotlib as STL to open in ParaView."
description: "This blog post details the steps to convert a surface plot into STL."
categories: [coding, presentation]

typora-root-url: ../../../../website
---

## The problem

I need to convert a surface plot generated in python to STL format for better visualization and modification in ParaView.

## The solution

Use [this library](https://github.com/iitrabhi/surf2stl-python) to convert your data into STL format. The process is quite simple. Surface plots in python follow a mesh-grid structure. Thus we can just feed the same data to the `surf2stl.write()` function to write our data to an STL format. Once written we can visualize and modify it in ParaView. Copy and paste the source file from [this repo](https://github.com/iitrabhi/surf2stl-python) into your work directory.

```python
import numpy as np
import surf2stl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
x = np.linspace(-6, 6, 30)
y = np.linspace(-6, 6, 30)
X, Y = np.meshgrid(x, y)
Z = np.sin(np.sqrt(X ** 2 + Y ** 2))
```

### Matplotlib plot

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface.
ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)
plt.show()
```

### Write to STL

```python
surf2stl.write('3d-sinusoidal.stl', X, Y, Z)
```

## Output

![image-20220304210941963](/assets/images/image-20220304210941963.png)

<figcaption>`left` matplotlib, `right` ParaView</figcaption>

## Note

If the data does not obey mesh-grid structure then we need to perform triangulation on the point cloud and generate a triangulated surface. Post triangulation either we can use [pygmsh](https://github.com/nschloe/pygmsh) or [pyvista](https://docs.pyvista.org/) to write the mesh in XDMF or VTU.