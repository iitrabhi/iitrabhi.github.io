---
layout: post
title: "What to do when linux can not find an installed python package"
description: "This post explains a little bit about linux $PATH variable and its use for python packages."
tag: 
  - [coding,home]
typora-root-url: ../../../../website
---

## The problem

Linux is tricky, and many times the `$PATH` gets corrupted because of some new installation. So, I would like to keep a reminder here on the website on the process of rectifying this issue.

## The Solution

First we need to understand the working of python and pip. pip is the package manager for Python packages. This is like an app store for python modules. Just like we go on the app store to find any app of our liking and install in from there, for python packages we go to [pip-store](https://pypi.org/) and install the packages with the command 

``` 
pip install package_name==version_number
```

In Linux, the package gets installed at the following location

- `/usr/local/lib/python3.6/site-packages/` or
- `/usr/local/lib/python3.6/dist-packages/`

You have to update `python3.6` in the above locations to the version number of python on your system. Suppose you have `python3.9` on your system then the installation path will be

- `/usr/local/lib/python3.9/site-packages/` or
- `/usr/local/lib/python3.9/dist-packages/`

Now the computer needs a way to find this package on the system. For this in LINUX we have the `$PATH` variable which is a system variable. You can think of it as the app drawer or the start menu in windows where you can find all the installed apps. 

`$PATH` variable contains the location of all the directories where system can find an installed package. You can see the content of `$PATH` variable by typing in terminal the following command

```
$PATH
```

On my system I get the following output

```mk
bash: /home/fenics/local/lib:/usr/local/lib/:/usr/lib/ccache:/usr/local/gmsh-4.3.0-Linux64/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
```

As you can see from above, I donot have the pip install directories on my `$PATH` variable. Thus, my system is unable to find the installed package. To rectify it, I just need to run the following command

``` 
export PATH=$PATH:/usr/local/lib/python3.6/dist-packages/
```

Now I can see my `$PATH` variable which gives the output

```mk
bash: /home/fenics/local/lib:/usr/local/lib/:/usr/lib/ccache:/usr/local/gmsh-4.3.0-Linux64/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/lib/python3.6/dist-packages/
```

And now I can use the packages installed via `pip`.

## References

- [Update path variable](https://opensource.com/article/17/6/set-path-linux)
- [Where does pip install its packages?](https://stackoverflow.com/questions/29980798/where-does-pip-install-its-packages)
