---
layout: post
title: "How to use SVG images on your website with light and dark mode."
description: "This post contains simple examples to showcase the way I use SVG in html."
tag: 
  - code
typora-root-url: ../../../../website
---

I love vector-based drawing since it offers infinite zoom and is free from any kind of pixelation artifact. Moreover, for simple scientific drawings, they are substantially smaller in size as compared to their raster counterpart. I have been using vector format for all of my drawings in my research. Now I would like to do the same on the website. Adding a vector-based image to the website is quite easy. You just wrap the image in the `<object>`tag.

```html
<object data="{{ site.url }}/assets/svg/abhi.svg" width="80%" style="pointer-events: none;"></object>
```

This code renders the SVG image on my website. Currently, I am using Jekyll with Liquid and thus, the syntax for the path contains the variable- `site.url`. Setting `pointer-events` to `none` will allow setting `href` to the image. You can find more about it on the internet.

---

The above code snippet successfully sets my desired image at my desired location. Now, I would like the image to change with the theme. That means in the dark mode, the image must be light, and in the light mode, the image must be dark. This could be achieved by using the `filter` method.

We can use the `data-theme` property to set the dark and light filters. In the dark mode, I want the filter to invert the image from black to white. Thus, I will set the invert intensity to 100%. In the light mode, I want my image to be its standard color, i.e. black. Thus I will set the filter intensity to 0%. Then we can define a class that will switch the filter based on the variable and assign that class to the object. This way I will have light and dark icons on the website.

```css
html[data-theme="dark"] {--invert-intensity: invert(1);}
html,html[data-theme="light"] {--invert-intensity:invert(0);} 

.abhi-icon{
    filter: var(--invert-intensity);
 }
```

```html
<object class="abhi-icon" data="{{ site.url }}/assets/svg/abhi.svg" width="80%" style="pointer-events: none;"></object>
```

