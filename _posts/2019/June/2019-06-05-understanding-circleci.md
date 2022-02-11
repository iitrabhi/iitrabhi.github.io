---
layout: post
title: "Understanding CircleCI"
description: "Continuous testing of the code can only guarantee its correctness."
tag: 
  - gsoc
typora-root-url: ../../../../website
---

The official [about us page of CircleCI](https://circleci.com/about/) states the following:

> CircleCI allows teams to rapidly build quality projects, at scale. Our mission is to give people everywhere the power to build and deliver software at the speed of imagination.

Their package allows us as a developer to develop on the principles of [continuous integration](https://circleci.com/continuous-integration/).

> Continuous integration is a software development strategy that increases the speed of development while ensuring the quality of the code. Developers continually commit small increments of code (at least daily, or even several times a day), which is then automatically built and tested before it is merged with the shared repository.

The official repo of FEniCS uses CircleCI to test the commits by different developers working on the repository before merging them into the master branch. In this post, I would like to introduce you to the concept of CircleCI and how I have incorporated it into [my GSoC repository](https://github.com/iitrabhi/GSoC2019).

The first step is to create a “.circleci” folder in the root directory of your project. Inside that directory, we have to create a “config.yml” file. Inside “config.yml” we would load our commands that would then be utilized by CircleCI to build and test our code. After configuring our “config.yml” we will connect our repository with our CircleCi account. More on it [here](https://circleci.com/docs/2.0/getting-started/). That’s it. Now every time we commit our code to our Github repo, CircleCI would automatically build and run tests on our code.

------

CircleCI can build your code inside a docker image. You can read more about docker image [here](https://computationalmechanics.in/understanding-the-dockerfile-of-dolfin-x/). When you ask CircleCI to build and test your repository it essentially performs the following steps:

- Loads the specified docker image.
- Loads the contents of your repository into the root of the docker image.
- Build and run your repository inside docker image.

------

First, we start with specifying which image we would like to use for our build. For my own repository, I am using the official dolfinx image. This command alone would load the docker image and my repository into the root of that image.

```
version: 2
jobs:
  build:
    docker: 
      - image: quay.io/fenicsproject/dolfinx:real 
```

Next, we could just simply specify the commands that we wish to execute inside the image. Here, we first get inside the python directory inside directory Scripts and then run the program “understanding_mesh_workflow.py”.

```
steps:
      - checkout # check out the code in the project directory
      - run: |
             cd Scripts/python
             python3 ./understanding_mesh_workflow.py
```

This simple config file is all that was required to test the code using the concepts of CI. You can see the build output [here](https://circleci.com/gh/iitrabhi/workflows/GSoC2019).
