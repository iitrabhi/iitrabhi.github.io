---
layout: post
title: "Load the project to PyPi"
tag: 
  - code
typora-root-url: ../../../website
---

I have developed so many tools in the lab (in python) for easing my workflow. Now I wish to convert those to opensource packages so that I can easily install them in the future from PyPi and maybe they will become useful for some other people. 

- The process to convert my script to installable format starts with finding a good name for the script. For that I need to go to the [following url](https://pypi.org/search/) and search whether that name is already taken or not. Fortunately for me the name `meshx` was not taken so I am using the same. 

- Next, we need to make a GitHub repository with the same name. [Here is mine](https://github.com/iitrabhi/meshx)

- Then we have to clone the repo to our development machine and start adding our code. PyPi expects our code to be in a particular format that is available at this [Real python tutorial](https://realpython.com/pypi-publish-python-package/)

- To test the code install it in editable format

  ```bash
  pip install -e .
  ```

- Once the code is thoroughly tested, it is time to push the code to PyPi

- In the root folder of project type the following

  ```bash
  pip install twine
  python setup.py sdist bdist_wheel
  twine upload dist/*
  ```

- You can now install the project by using the following command

  ```bash
  sudo -H pip install meshx
  ```

  

---

## Summary

Creating a repository on PyPi is quite straight forward if we follow the steps correctly. I have added the references below that I have used to make this project. This project is quite small and anyone can read through all the files to get a feel for making there own project on PyPi.

---

## References

- [Real python tutorial](https://realpython.com/pypi-publish-python-package/)
- [Developing command line tools with argparse](https://realpython.com/command-line-interfaces-python-argparse/)

