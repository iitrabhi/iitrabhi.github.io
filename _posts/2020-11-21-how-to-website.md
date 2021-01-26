---
layout: post
title: "How to setup the Moving theme for your Jekyll blog."
tag: 
	- misc
---

Github pages is a very simple way to convert your markdown files into static website. The workflow to publish your pages to the web is similar to the way we push our code to Github. That's why, once setup, we can work locally on our laptop in markdown and then push the code to Github, all the formatting will be taken care by the theme. Getting started with Github pages is quite simple, we have to just follow the simple steps written here:

> https://pages.github.com/

Once you have created your repository named '*username***.github.io**' you can go to that url to see your website. Right now that would be empty, but still following these simple steps you have your website up and running on the internet. Next, we have to setup **Jekyll** on our laptop. For this again follow the simple steps given here:

> https://jekyllrb.com

This will install **Jekyll** on your laptop and then we are ready to build our website. Clone the repository that you have just created on Github to your laptop and then *cd* to that in Terminal. Now here either you can setup the default Jekyll website or clone some existing website. In my case I cloned the website [Moving](https://github.com/huangyz0918/personal-page-blog) to the root of my repository and made the following changes.

- Deleted all the posts from the **_posts** folder.
- Updated the **Gemfile** and **Gemfile.lock**. You can directly use the files from my repository.
- Update the **_config.yml** file. Here you have to make the blog your own. Fill in your details. Remember to change the url to your own website, which in my case it "https://iitrabhi.github.io". This step is very important, without it Github will not render your website.

That's it. Now you can add posts to the **_posts** folder in markdown and then push them to your website. You have to follow the [naming convention and header syntax](https://jekyllrb.com/docs/posts/)  for the posts to be formatted properly.