---
layout: post
title: "How to add nomenclature to latex document."
categories: [latex, presentation]
typora-root-url: ../../../website
---

## The problem 

Even after installing the `nomencl.sty` , I am not getting nomenclature in my pdf document. The reson for this is that to print nomenclature, it requires use to run the following command in the Terminal.

```bash
makeindex paper.nlo -s nomencl.ist -o paper.nls
```

![image-20211108111443455](/assets/images/image-20211108111443455.png)

Once the `.nls` is build we can run `pdflatex` normally to get nomenclature. So the proper order of commands to run is:

- `pdflatex paper.tex`
- `makeindex paper.nlo -s nomencl.ist -o paper.nls`
- `pdflatex paper.tex`

*Note: All of the above commands are to be run inside the terminal. First navigate to the directory where you have you source latex file (`paper.tex`) and then run the above commands.*

## References

- [stackexchange](https://tex.stackexchange.com/questions/62061/problem-with-the-nomenclature)

