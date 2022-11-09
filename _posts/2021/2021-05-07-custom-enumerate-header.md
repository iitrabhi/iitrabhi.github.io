---
layout: post
title: "Setting custom header for list of items | Step 1:"
categories: 
  - [latex, presentation]
typora-root-url: ../../../website
---

Use the following code to get custom headers

```latex
\begin{enumerate}[label=\bfseries Step \arabic*:,leftmargin=*]
\item Discretize the domain and define the finite element problem with interval variables.
\end{enumerate}
```

The above command would result it

![image-20210507113158212](/assets/images/image-20210507113158212.png)

Another version

```latex
\begin{enumerate}[{Step }1:] 
\item Discretize the domain and define the finite element problem with interval variables.
\end{enumerate}
```



