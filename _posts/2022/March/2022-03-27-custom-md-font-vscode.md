---
layout: post
title: "How to have a custom markdown preview font-style and size for VSCode."
description: "This blog post details the settings required to have a custom VSCode fonts in markdown preview."
categories: 
  - [coding]
typora-root-url: ../../../../website
---
## The problem

Recently I have shifted to VSCode from sublime, and I like most things here. It's much more powerful than sublime but seems as light as sublime. The inbuild markdown previewer is very useful, and even though I use Typora for writing the blog posts for my website, I love this nifty feature. The only minor thing I want is to have the same fonts that I use on my website and Typora to be rendered in the markdown preview of VSCode.

## The solution

The first thing to do is to make a `CSS` file that would detail the layout properties of `MD`. You can download my [current settings](https://github.com/iitrabhi/iitrabhi.github.io/blob/master/markdown.css) and modify them to your liking. Add this file to the root directory of your project.

The next thing to do is save your project as a workspace in VSCode. That would create the `*.code-workspace` file in your root directory. Open that and paste the following to control the font size of the markdown preview.

```json
{
	"folders": [
		{
			"path": "."
		}
	],
	"settings": {"markdown.styles": [
		"markdown.css"
	],
	"markdown.preview.fontSize":16}
}
```

These settings will change the font size of the markdown for the current workspace and set the styling as per the `markdown.css` file. Here is the output with 'Bitter' font for the headings and paragraphs.

![Screenshot 2022-03-27 at 12.13.04 PM](/assets/images/Screenshot%202022-03-27%20at%2012.13.04%20PM.png)

This kind of synchrony gives a sense of WYSIWYG while writing the blog posts and I like it. 

## References

- [This gave me an idea of implementation](https://github.com/raycon/vscode-markdown-style)
