---
layout: post
title: "How to have a clean latex folder"
tag: 
  - latex
typora-root-url: ../../website
---

I am using sublime as the latex editor of choice and I really like it. The only problem that I had with sublime was that it was creating all the output files and auxilary files in the root directory. This was creating a mess. I like to keep things clean and minimal. Today I found a solution. We just need to set the `aux_directory` and the `output_directory` in the sublime-project file. Sublime will then create all the respective folders and files. If you need to clean your directory you only need to delete the `cache` folder.

```json
{
	"folders":
	[
		{
			"file_exclude_patterns":
			[
				"*.cls",
				"*.md",
				"*.gitignore",
				"*.sublime-project"

			],
			"path": "."
		}
	],
	"settings":
	{
		"TEXroot": "thesis.tex",
		"ensure_newline_at_eof_on_save": true,
		"spell_check": true,
		"trim_trailing_white_space_on_save": true,
		"word_wrap": true,
		"wrap_width": 100,
		"aux_directory": "cache",
		"output_directory": "cache"
	}
}
```

