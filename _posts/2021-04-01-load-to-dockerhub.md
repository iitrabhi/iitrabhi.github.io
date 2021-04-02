---
layout: post
title: "Load the working copy of your docker image to docker hub for reuse"
tag: 
  - code
typora-root-url: ../../website
---

I have been working with a team on few projects based on FEniCS. Time and again we have to build our images on different systems and that requires for all the lab members to have an updated version of the Dockerfile with them. I have placed the docker file on github [here](https://github.com/iitrabhi/fenics-docker). This procedure was still creating a lot of problems and thus I have now decided to load the images to dockerhub so that my team mates can pull them easily in just one line. The documentation provides us with all the steps that are necessary to load an image to dockerhub. This was everyone from the team will have the same version of all the software packages installed on their systems. The steps are as follows:

- First we need to build the image in our own system using the command

  ``` 
  docker build --target base -t fenics .
  ```

- Then we have to login to our dockerhub account using

  ```
  docker login --username=yourhubusername
  ```

- Next we have to tag our image with a name

  ```
  docker tag bb38976d03cf iitrabhi/fenics:latest
  ```

- Then we have to push the image to docker hub

  ```
  docker push iitrabhi/fenics
  ```

That is it. These are the steps that are necessary for uploading the image to dockerhub. Now anyone can pull this image with

``` 
docker pull iitrabhi/fenics:latest
```

---

## Summary

The process of uploading an image that you have build on your system is quite simple as is described in this post. Next, I would like to learn about the auto-builds that can happen when we push a code to GitHub.

---

## References

- [Official documentation](https://ropenscilabs.github.io/r-docker-tutorial/04-Dockerhub.html)

