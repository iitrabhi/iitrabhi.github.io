---
layout: post
title: "How to set up your first JavaScript project."
categories: [coding]
tags: [js, javascript]
typora-root-url: ../../../website
---

 

- Create a new project in GitHub and clone that to your desktop. In my case that is the `nurbs` project.

- Create an `index.html` file inside the project. Usually this is the file that is used my many to store the primary code of the website or the entry point code of the website.

- An HTML file is made of two parts. The `head` and the `body`. From what I understand the `head` is the part where all the `scripts` and `style` live and the `body` is the part where the primary layout of the website live. 

  ```html
  <!DOCTYPE html>
  <html>
  <head>
      <!-- Load the Paper.js library -->
      <script type="text/javascript" src="node_modules/paper/dist/paper-full.js"></script>
      <script type="text/paperscript" src="line.js" canvas="myCanvas">
      </script>
      <style type="text/css">
      body {height: 100%;}
      /* Scale canvas with resize attribute to full size */
      canvas[resize] {
          width: 100%;
          height: 100%;
      }
      </style>
  </head>
  <body>
      <canvas id="myCanvas" resize></canvas>
  </body>
  </html>
  ```

  

- In this example I wish to learn how to run a script on the website and also how to use `paper.js`.

- The first step is to open the directory in `cmder` and run the following command

  ```bash
  npm install paper
  ```

  This will install paper in the repository. This simply means that it will copy the source code of the paper.js into the repository which we can then access locally.

- Now once the initial setup is done the folder should look something like this
  ![image-20210214184907365](/assets/images/image-20210214184907365.png)

- The directory `node_modules` contains the code for  `paper.js`. This is also referenced in the head of `index.html`. This reference is somewhat similar to `import` statement in python, as in it imports all the methods from `paper.js` for our use.

  ```html
  <!-- Load the Paper.js library -->
      <script type="text/javascript" src="node_modules/paper/dist/paper-full.js"></script>
  ```

- Next, we need to tell the browser where to look for the primary source code for our JavaScript. In my case I have added the JS code to `line.js`

  ````html
  <script type="text/paperscript" src="line.js" canvas="myCanvas">
  ````

  paper.js will draw everything on a canvas object that we will create inside the `body` of our page. In the above line the script tag will link the PaperScript code to the canvas with the given ID and produces a `Project` and `View` object for it on the fly. So here we are just telling the browser that use the canvas `myCanvas` to draw using the `line.js` script.

- Finally in the body of the HTML we will create a canvas object and give it the id `myCanvas`

  ```html
  <body>
      <canvas id="myCanvas" resize></canvas>
  </body>
  ```

- Ok, so now the page is ready. To start adding paper object to the code we can simply update our `line.js` file. Initially I will just try out the example given on their website.

  ```javascript
  // Create a circle shaped path with its center at the center
  // of the view and a radius of 30:
  var path = new Path.Circle({
      center: view.center,
      radius: 30,
      strokeColor: 'black'
  });
  
  function onResize(event) {
      // Whenever the window is resized, recenter the path:
      path.position = view.center;
  }
  ```

- Now to finally view the page I need to setup a local server [because read this](https://developer.mozilla.org/en-US/docs/Learn/Common_questions/set_up_a_local_testing_server). The way to do is it simply go to the project root folder and run

  ```bash
  python -m http.server
  ```

- This will create a local server that I can access at `http://localhost:8000/`

- Now I can actively debug the project. To debug I have to press `ctrl + shift + i` and go to console tab.