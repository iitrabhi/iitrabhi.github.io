---
layout: post
title: "How to set sublime as text editor in gmsh."
categories: [productivity]
tags: [gmsh, sublime]
typora-root-url: ../../../website
---

Go to `Tools` → `options` → `General` → `Advanced`

![image-20210122154513524](/assets/images/image-20210122154513524.png)

`Windows` Paste the following in Text Editor command

```
C:\Program Files\Sublime Text 3\subl.exe %s 
```

on a `Mac` type the following command

```
/Applications/Sublime\ Text.app/Contents/SharedSupport/bin/subl %s
```

On `Linux` type

```
subl %s
```

Goto sublime and open package manager and install `gmsh-Tools`

![image-20210309103738751](/assets/images/image-20210309103738751.png)

After installation open a geo file and on the bottom right click on `Plain text `and go to Open all with current extension as ... → gmsh.![image-20210309104045408](/assets/images/image-20210309104045408.png)

Now you will get color coded gmsh files

![image-20210309104220484](/assets/images/image-20210309104220484.png)

