---
layout: post
title: "Creating interactive plots with python and plotly."
description: "This post contains some example interactive plots created with python and plotly."
categories:  [coding, presentation]
tags: [plotly, python, plot] 
typora-root-url: ../../../../website
---

## Installation

```
pip install plotly
pip install cufflinks
```

## Example

```python
import plotly.express as px
import numpy as np

t = np.linspace(0, 2*np.pi, 100)

fig = px.line(x=t, y=np.sin(t), labels={'x':'t', 'y':'cos(t)'})
fig.show()
```

{% include_relative assets/sin.html %}

---

```python
import plotly.graph_objs as go
layout = go.Layout(
  margin=go.layout.Margin(
        l=0, #left margin
        r=0, #right margin
        b=0, #bottom margin
        t=0, #top margin
    )
)
t = np.linspace(0, 2*np.pi, 100)
X,Y = np.meshgrid(t,t)
Z = np.sin(X)*np.sin(Y)
fig = go.Figure(data=[go.Surface(z=Z, x=X, y=Y, colorscale ='rdbu')])
fig.update_layout(layout)
fig.show()
```

{% include_relative assets/sin3D.html %}