---
layout: post
title: "GUI desktop app for isogeometric analysis"
description: "The output of my first year of PhD coding IGA."
langs: [python, webgl, PYQT]
year: "November 2018"
date: 2018-11-01
typora-root-url: ../../../../website
---

Iso-geometric analysis (IGA) is a relatively new method that is a logical extension and generalisation of the classical finite element method (FEM). The idea of IGA was to merge the gap between computer aided design (CAD) and computer aided engineering (CAE) into one model by using a unified geometric representation.. When the researchers started exploring the concept, the primary intention was to perform finite element simulation in the 3D modeling package. But, when I began exploring software packages that allowed IGA, I was disappointed to find out that the only package that supported IGA was not active anymore. Thus, I started development of my own package.

To facilitate rapid prototyping and testing of new IGA based elements a GUI based framework to carry out IGA has been developed. This work is still under progress. The aim with this project is to create a scalable object oriented implementation of IGA. The figure below presents working of the in-house developed GUI based framework to carry out IGA for a plate bending problem.

![IITR-IGA1](/assets/images/IITR-IGA1.png)

![IITR-IGA2](/assets/images/IITR-IGA2.png)

![IITR-IGA3](/assets/images/IITR-IGA3.png)

![IITR-IGA4](/assets/images/IITR-IGA4.png)

The following software, technologies and libraries have been utilized for the development of this framework.

- **Programming Language:** Python

- **Libraries:** PyQT, Numpy, Scipy, Matplotlib, Pandas, Plotly, Cufflinks, Fenics

- **Version Control System:** Git, GitHub, Bitbucket

- **Documentation System:** Sphinx, Markdown

- **Development Environment:** PyCharm

- **Visualization:** Paraview, WebGL (Web Graphics Library)

The package under development is cross-platform which implies that it will run on Microsoft Windows, Linux, and MacOS.
