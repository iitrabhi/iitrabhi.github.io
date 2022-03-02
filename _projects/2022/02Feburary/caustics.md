---
layout: post
title: "Caustics in FEniCS [WIP]"
description: "Recreating the winning entry of 3B1B SoME challenge."
langs: [python]
year: "Feburary 2022"
date: 2022-02-01
typora-root-url: ../../../../website
---

I participated in the [SoME one challenge 2020](https://www.youtube.com/watch?v=v_DqwOVr3Vw&feature=emb_title). It was a great experience to make a full video with python in Blender finally. One of the winning entries of the challenge was [this blog post](https://mattferraro.dev/posts/caustics-engineering) on caustics by [Matt Ferraro](). It was the first time that I heard of the term caustics. The blog post is well written, and I love how he has described the problem and the solution. 

Immediately after reading the blog post, I wanted to recreate the analysis portion in FEniCS and use the CNC machine in our lab to create the object. Moreover, after going through the [issues](https://github.com/MattFerraro/causticsEngineering/issues/10), I learned that we could model caustics in Blender. And thus, we can test out different results of the analysis run without even creating the thing in the real world.

The first step was to just run the code and test the generated model in Blender. Here is the 3D model setup. It's crude right now but I did get the results.

![image-20220217234453361](/assets/images/image-20220217234453361.png)

The image generated from this experiment is a good first try. The next step is to get a clean image in Blender with the parameters presented in the blog post.

![image-20220218000257522](/assets/images/image-20220218000257522.png)

`[Update 25/2/22]` The initial machining failed to produce any image with caustics. After reading the original blog post, I figured out one key difference, the toolpath for roughing and finishing operations must be parallel. Also, we had scaled the model to fit the dimensions. The next try will be with actual dimensions.

`[Update 2/3/22]` With some tweaks, I got the following image from the blender simulation. This simulation is with the 100mm $$\times$$ 100mm size and 5mm thickness of the acrylic sheet. The total cut depth is 3mm. Will try to machine it tomorrow and update the results.

![image-20220302214214166](/assets/images/image-20220302214214166.png)