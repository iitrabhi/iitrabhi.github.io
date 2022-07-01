---
layout: post
title: "Different math environments available in LaTeX"
description: "A picture of all the math environments available in latex"
tag: 
  - latex
typora-root-url: ../../../../website
---

![image-20220430145635931](/assets/images/image-20220430145635931.png)

```latex
\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{bm}
\usepackage{amsfonts}
\usepackage{mathrsfs}
\usepackage[margin=3cm]{geometry}

\title{math_symbols}
\author{}
\date{}

\begin{document}
\section{lowercase}

\begin{tabular}{r|l}
    mathbf & $\mathbf{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathit & $\mathit{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathnormal & $\mathnormal{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathrm & $\mathrm{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathsf & $\mathsf{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathtt & $\mathtt{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    boldsymbol & $\boldsymbol{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathfrak & $\mathfrak{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathbb & $\mathbb{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathscr & $\mathscr{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
    mathcal & $\mathcal{a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z}$ \\
\end{tabular}

\section{UPPERCASE}

\begin{tabular}{r|l}
    mathbf & $\mathbf{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathit & $\mathit{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathnormal & $\mathnormal{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathrm & $\mathrm{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathsf & $\mathsf{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathtt & $\mathtt{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    boldsymbol & $\boldsymbol{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathfrak & $\mathfrak{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathbb & $\mathbb{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathscr & $\mathscr{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
    mathcal & $\mathcal{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}$ \\
\end{tabular}

\end{document}

```



