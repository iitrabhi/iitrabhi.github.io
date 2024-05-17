---
layout: post
title: "Learning javascript basics: getElementById"
description: My notes on javascript
categories:
  - coding
tags:
  - javascript
---


The `getElementById` method is a widely used function in web programming, specifically in JavaScript. It is part of the Document Object Model (DOM) and is used to access and manipulate elements within an HTML document. The method returns the element that has the ID attribute with the specified value. If no such element exists, it returns `null`.

Here are some common uses and methods associated with an element obtained using `getElementById`:

1. **Modifying Text Content**
   - **innerText**: Sets or returns the text content of the specified element.
   - **textContent**: Sets or returns the text content of the specified element, including all whitespace and HTML tags.

   ```javascript
   var element = document.getElementById("demo");
   element.innerText = "Hello, World!";
   ```

2. **Changing Style**
   - **style**: Allows you to change the style of an element by setting properties like `color`, `backgroundColor`, `fontSize`, etc.

   ```javascript
   var element = document.getElementById("demo");
   element.style.color = "red";
   ```

3. **Modifying HTML Content**
   - **innerHTML**: Sets or returns the HTML content (inner HTML) of an element.

   ```javascript
   var element = document.getElementById("demo");
   element.innerHTML = "<strong>Updated Content</strong>";
   ```

4. **Adding Event Handlers**
   - **addEventListener**: Adds an event listener to the specified element.

   ```javascript
   var element = document.getElementById("demo");
   element.addEventListener("click", function() {
       alert("Element Clicked!");
   });
   ```

5. **Attributes Manipulation**
   - **setAttribute**: Sets a new value for an attribute on the specified element.
   - **getAttribute**: Returns the value of a specified attribute on the element.

   ```javascript
   var element = document.getElementById("demo");
   element.setAttribute("data-custom", "12345");
   console.log(element.getAttribute("data-custom"));
   ```

6. **Class Manipulation**
   - **classList**: Provides methods like `add`, `remove`, `toggle`, and `contains` for class manipulation.

   ```javascript
   var element = document.getElementById("demo");
   element.classList.add("new-class");
   element.classList.remove("old-class");
   element.classList.toggle("active");
   ```

These methods and properties provide a robust interface for interacting with and manipulating the elements of a webpage dynamically using JavaScript.
