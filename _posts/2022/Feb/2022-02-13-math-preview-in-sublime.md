---
layout: post
title: "Previewing math equations in Sublime with custom commands."
description: "The following latex package will allow you to do so."
tag: 
  - latex
typora-root-url: ../../../../website

---

## The problem

When working on a LaTeX project in sublime, if I have defined custom commands with the `\newcommand` command, then the math preview of `LaTeXTools` fails to preview them in sublime. I would like to have the full preview of the equations with custom commands in sublime itself.

## The solution

- Install the package [LaTeXYZ](https://packagecontrol.io/packages/LaTeXYZ).

- GoTo user prefrences of LaTeXYZ and set the following command to true.

  ```json
  {"auto_set_preview_math_template_preamble": false}
  ```

## Reference

[Question on Stack overflow](https://stackoverflow.com/questions/42284544/st3-latex-loading-locally-defined-commands-for-in-line-live-preview-of-math)