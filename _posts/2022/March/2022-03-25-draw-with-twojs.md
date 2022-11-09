---
layout: post
title: "How to draw curves with two js."
description: "This blog post details my understaning of drawing with Two js."
tag: 
  - [coding, presentation]
typora-root-url: ../../../../website

---
## The motivation
I am trying to learn Twojs to make interactive widgets to share my knowledge of FEM and IGA. The simplest case is to draw a line and point on the canvas.

```js
var ob1 = document.getElementById('canvas_1');
var two1 = new Two();
two1.appendTo(ob1);

var curve = two1.makeCurve(0, 0, 220, 150, 240, 250, 260, 150, 280, 250, 290, 200, true);
curve.linewidth = 2;
curve.scale = 1.0;
// curve.rotation = Math.PI / 2;
curve.noFill();
curve.stroke = 'rgba(255, 0, 0, 0.5)';

var line = two1.makeLine(150, 150, 350, 350);
line.linewidth = 2;
// line.fill = "#881111";
line.noFill();
line.stroke = "rgba(255, 0, 0, 0.5)";
two1.update();
```

## References

- https://code.tutsplus.com/tutorials/drawing-with-twojs--net-32024