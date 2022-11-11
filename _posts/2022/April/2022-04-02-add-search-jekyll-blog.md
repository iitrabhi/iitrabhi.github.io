---
layout: post
title: "How to add search functionallity to jekyll blog."
description: "Use this simple method to add search functionallity to your jekyll blog."
categories: [coding]
tags: [jekyll, search, liquid]
typora-root-url: ../../../../website
---

## The problem

I have curated a lot of information on this website and would like to have search functionality to make it easy to retrieve the data.

## The solution

The first step is to create an index of all the posts on the website. This could be done automatically with the help of a small `Liquid` script. The full process is detailed on [this blogpost](https://blog.webjeda.com/instant-jekyll-search/). Copy and paste the following to the root directory of the project in a file named `search.json`.

```json
{% raw %}---
---
[
  {% for post in site.posts %}
    {

      "title"    : "{{ post.title | strip_html | escape }}",
      "url"      : "{{ site.baseurl }}{{ post.url }}",
      "category" : "{{post.categories | join: ', '}}",
      "tags"     : "{{ post.tags | join: ', ' }}",
      "date"     : "{{ post.date }}",
      "discription" : "{{post.description | strip_html | strip_newlines | escape }}"

    } {% unless forloop.last %},{% endunless %}
  {% endfor %}
]{% endraw%}
```

The above script will automatically create an index of all the posts on the website.

Next, we need to make a search box somewhere on the website and connect that to the `search.json` index file generated. For this paste the code to the page on which you wish to use search functionallity.

```html
<!-- Html Elements for Search -->
<div class="form-group" id="search-container">
        <input class="form-field" type="text" id="search-input"  placeholder="Search">
        <!-- <span>Search</span> -->
</div>
<div>
    <ul id="results-container"></ul>
</div>

<!-- Script pointing to search-script.js -->
<script src="{{ base.url | prepend: site.url }}/assets/js/search_script.js" type="text/javascript"></script>

<!-- Configuration -->
<script>
SimpleJekyllSearch({
  searchInput: document.getElementById('search-input'),
  resultsContainer: document.getElementById('results-container'),
  json: '/search.json'
})
</script>
```

Add the CSS styling for the form to the `assets/css` folder and link it to the webpage. Finally add the `search_script.js` file to the `assets/js` folder. You can see the script in action [here](https://abhigupta.io/tags).

## References

- [This blog has detailed instructions](https://blog.webjeda.com/instant-jekyll-search/)
- [Search Script](https://raw.githubusercontent.com/christian-fei/Simple-Jekyll-Search/master/dest/simple-jekyll-search.min.js)
