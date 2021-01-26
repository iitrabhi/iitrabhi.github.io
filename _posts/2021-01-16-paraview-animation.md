---
layout: post
title: "How to export animation from paraview."
tag: 
	-  tools
typora-root-url: ../../website
---
Here are my steps to create a video file of the animation from Paraview.

- Select the proper color code and position of the plot in Paraview.  Paraview will export the animation as it is presented on the application
  ![image-20210115171028235](/assets/images/image-20210115171028235.png)

- Now goto File  →  Save Animation and save the animation  into a folder of your choice. Save the animation as `.png` images.

  ![image-20210116143746320](/assets/images/image-20210116143746320.png)

- Next open `Video Editor` app that comes bundled with Windows and click on `New Project`

  ![image-20210116143919343](/assets/images/image-20210116143919343.png)

- This will create a new project on your computer. Next, drag and drop all the images that you have exported from Paraview into the `Project library` area.

  ![image-20210116144145988](/assets/images/image-20210116144145988.png)

- Drag all the images from project library to timeline.
  ![image-20210116144238957](/assets/images/image-20210116144238957.png)

- Select all the images in timeline by pressing `Ctrl + A` on your keyboard twice. Now, click on the duration button and change the duration of each image to $$\frac{1}{frame-rate} seconds$$. That is if you want the frame-rate to be 20 fps then the image duration will be 1/20. ![image-20210116144612334](/assets/images/image-20210116144612334.png) 

- That is it. Now you can click on `Finish video` on top right and export the video.![image-20210116144720280](/assets/images/image-20210116144720280.png)