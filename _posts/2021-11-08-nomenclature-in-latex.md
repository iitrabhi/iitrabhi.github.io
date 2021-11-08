---
layout: post
title: "How to add nomenclature to latex document."
tag: 
  - tools
typora-root-url: ../../website
---

## The problem 

Even after installing the `nomencl.sty` , I am not getting nomenclature in my pdf document. The reson for this is that to print nomenclature, it requires use to run the following command.

```bash
makeindex paper.nlo -s nomencl.ist -o paper.nls
```

Once the `.nls` is build we can run `pdflatex` normally to get nomenclature. So the proper order of commands to run is:

- `pdflatex paper.tex`
- `makeindex paper.nlo -s nomencl.ist -o paper.nls`
- `pdflatex paper.tex`

## References

- [stackexchange](https://tex.stackexchange.com/questions/62061/problem-with-the-nomenclature)

