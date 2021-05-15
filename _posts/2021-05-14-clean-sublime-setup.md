---
layout: post
title: "How to have a clean latex folder"
tag: 
  - tools
typora-root-url: ../../website
---

# Download

- Download Sublime text here: https://www.sublimetext.com/
- Download Miktex here: https://miktex.org/download
- Download Sumatra PDF here: https://www.sumatrapdfreader.org/download-free-pdf-viewer

# Install

- Install Sumatra first and then Miktex. Finally install Sublime.

- Open the Command Palette: Press `Ctrl+Shift+P`

- Type ‘install’ in the Command Palette input box, which should autocomplete to ‘Install Package Control.’ Press Enter to select it.

- Sublime Text 3 will start installing Package Control. This may take a short while. Once installed, a pop-up will display the message: Package Control was successfully installed.

- Go to `Preferences  → Package Settings  →  Package Control  → Settings` and paste the following

  ```json
  {
  	"bootstrapped": true,
  	"in_process_packages":
  	[
  	],
  	"installed_packages":
  	[
  		"A File Icon",
  		"Agila Theme",
  		"AutoPEP8",
  		"ayu",
  		"DistractionFreeWindow",
  		"Dracula Color Scheme",
  		"Fold Comments",
  		"ImageMagick",
  		"LaTeX Word Count",
  		"LaTeX-cwl",
  		"LaTeXTab",
  		"LaTeXTools",
  		"Non Text Files",
  		"Package Control",
  		"Python 3",
  		"SideBarEnhancements",
  	]
  }
  ```

- Save the file. This will automaticaly install all the packages necessary for the setup.  Wait for 5-10 mins for the installation to complete.

- Next open `Preferences  → Settings`  and paste the following there.

- ```json
  {
  	"auto_complete_triggers":
  	[
  		{
  			"characters": ".",
  			"selector": "source.python - string - comment - constant.numeric"
  		},
  		{
  			"characters": "\\",
  			"selector": "text.tex.latex"
  		}
  	],
  	"color_scheme": "Packages/Color Scheme - Default/Mariana.sublime-color-scheme",
  	"default_line_ending": "unix",
  	"font_size": 14,
  	"ignored_packages":
  	[
  		"Vintage"
  	],
  	"open_externally_patterns":
  	[
  		"*.jpg",
  		"*.jpeg",
  		"*.png",
  		"*.gif",
  		"*.zip",
  		"*.pdf"
  	],
  	"rulers":
  	[
  		100
  	],
  	"tab_size": 4,
  	"theme": "Agila.sublime-theme",
  	"translate_tabs_to_spaces": true
  }
  ```

- Goto `Preferences→` and paste the following 

  ```json
  [
      { "keys": ["f1"], "command": "toggle_side_bar" },
      { "keys": ["f2"], "command": "distraction_free_window" },
  ]
  ```

- Now finally go to `Preferences→Package settings→Latex tools→Check system` to check whether everything is fine or not.

