---
layout: post
title: "How to make tutorials for your class."
categories: [productivity, presentation]
tags: [python, teaching, jupyter, notebooks]
typora-root-url: ../../../website
---

## The problem 

I have always found the standard way of teaching just via textbooks soo limiting. When I was trying to understand the principles behind shape functions in FEM, 3D visualization tools were of great help. I wish to teach the students in the way I have learned, but I also understand that it is a very daunting task for students - who are not well versed in the field of programming - to code stuff on their own. That is why I find Jupyter notebooks really interesting and a solution to my problem. 

I am writing this post as I am creating [this repository](https://github.com/iitrabhi/fem) to help a student with FEM. 

## Benifits of Jupyter notebook

- Allows to write the theory of the problem with actual code execution. [See here.](https://nbviewer.jupyter.org/github/waltherg/notebooks/blob/master/2013-12-03-Crank_Nicolson.ipynb)
- Is easy to work with in Windows, Linux or MAC.
- Could be rendered online and thus shared with the whole class.
- Supports 2D as well as 3D graphs. Also, supports interactive graphs.
- Since it is based on python, learning in Jupyter opens up a whole new world of applications for the student.

## How to make a tutorial in Jupyter

- Make the repository in Github in which you wish to keep all your tutorials.

- Add a `README.ipynb` file to the repo and keep it in sync with the actual `README.md` file.

- Load your `README.ipynb`  in [nbviewer](https://nbviewer.jupyter.org/).

- Get the link to your rendered notebook and paste it in the Github repo `README.md`.

- Now, keep both the `README.ipynb` and `README.md` files updated with the source code. Think of it as an index for your tutorials.

- For every problem that you add to the repo, update the `Contents` with a link to that file. You need to adjust the link as follows


  **Github URL: -** https://github.com/iitrabhi/fem/blob/main/README.ipynb

  **Nbviewer URL:-**https://nbviewer.jupyter.org/github/iitrabhi/fem/blob/main/README.ipynb

  In your GitHub url, you need to just change `github.com` to `nbviewer.jupyter.org/github`

- Update both the `README.md` as well as `README.ipynb`

## References

- [Inspiration](https://github.com/mscroggs/bempp-acoustic-tutorials)
- [nbviewer](https://nbviewer.jupyter.org/github/waltherg/notebooks/blob/master/2013-12-03-Crank_Nicolson.ipynb)

