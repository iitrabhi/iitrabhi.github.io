---
layout: post
title: "How to properly add bib to your latex project using zotero."
categories: [latex, presentation]
typora-root-url: ../../../website
---

# The problem 

I have 100's of references in my zotero library. When I bring them to my `bib` file it ends up in a mess. There are many unwanted keys and there is no order to the references therein. This is the workflow that I have developed to handle the problem

# Workflow

- **Better bibtex:** Use the better-BibTeX plugin to clean up the export in Zotero. After installation of the plugin, you can set the fields to omit. Goto `Preferences ->Better Bibtex -> Export ` . I remove the following fields

  ```
  abstract, url, language, urldate, file
  ```

  ![image-20210701085733630](/assets/images/image-20210701085733630.png)

- **BibTeX Tidy:** This website is excellent to beautify your bib file and clean up further. I use it to have a proper indentation in my bib file and also to sort the entries by year.

  ![image-20210701090006256](/assets/images/image-20210701090006256.png)

# Before and After

![image-20210701090627770](/assets/images/image-20210701090627770.png)

# References

- [Better Bibtex](https://retorque.re/zotero-better-bibtex/)
- [Bibtex Tidy](https://flamingtempura.github.io/bibtex-tidy/)

