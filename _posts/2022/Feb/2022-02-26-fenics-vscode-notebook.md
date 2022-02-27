---
layout: post
title: "FEniCS code development in VSCode with jupyter notebook and Docker."
description: "This post describes my setup for development with FEniCS and Docker in VSCode."
tag: 
  - code
typora-root-url: ../../../../website
---

## The problem

I like to develop using a jupyter notebook for FEniCS. Right now my development workflow involves loading up my custom image in docker and starting a notebook server. Then I can use my web browser to open the notebook. Everything is fine with this workflow, but it has two primary limitations.

1. I do not like the theme of Jupyter Lab. I love working in the sublime text and am accustomed to the `Mariana` color scheme. 
2. We can not directly open `*.ipnby` files from the file browser. Sometimes I wish to just browse different notebooks in a folder without launching the container. This is right now not possible with my current workflow.

![image-20220226170047101](/assets/images/image-20220226170047101.png)

## The solution

I found out that VSCode supports the opening of `*.ipnby` files directly from the file browser. Once I searched for different color schemes available in VSCode, I was happy to find [Mariana Pro](https://marketplace.visualstudio.com/items?itemName=rickynormandeau.mariana-pro) theme for VSCode. This in itself was a great find and I love the outcome. Here is the same file in VSCode with Mariana theme.

![image-20220226170639588](/assets/images/image-20220226170639588.png)

So far so good. I can now directly open my notebooks with a double click to view their contents. But, what if I am able to attach the remote server to it and run the code 🤩.

---

Yes, we can attach the notebook to the remote server and, the process is quite simple. 

- First, we launch the notebook server as usual with the following command.

  ```
  docker run -p 8888:8888 -v D:\Codes\:/root/ -w /root/ iitrabhi/fenics_notebook
  ```

- Once launched we will get the remote server address inside the terminal. Just copy that address

- Now open VSCode and, on the bottom right of the application window click on the `Jupyter server` button

  ![image-20220226171418154](/assets/images/image-20220226171418154.png)

- This will open a command pallet. Click on the `Existing` option and then paste the remote address. 

  ![image-20220226171540958](/assets/images/image-20220226171540958.png)

- Now click on the `Select Kernel` option on the top right and select the first python from it.

  ![image-20220226171748655](/assets/images/image-20220226171748655.png)

- That is it. Now you can run your code by using the remote server's python.

Even though I am able to use the remote server, I am facing a small issue. In Jupyter Lab, when we open a notebook, it makes the directory of the notebook as current and will execute the code from that location. But, in the case of VSCode, it keeps the project directory as the root directory. This is problematic as to perform read and write operations from the current directory we need the latter.