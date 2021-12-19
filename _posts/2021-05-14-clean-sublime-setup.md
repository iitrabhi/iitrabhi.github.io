---
layout: post
title: "How to have a clean and minimal sublime setup + Dark Mode."
tag: 
  - tools
typora-root-url: ../../website
---

<iframe width="960" height="540" src="https://www.youtube.com/embed/SZEXdXL_P4w" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

# Download

- [Download Sublime text here.]( https://www.sublimetext.com/)
- [Download Miktex here.](https://miktex.org/download)
- [Download Sumatra PDF here.](https://www.sumatrapdfreader.org/download-free-pdf-viewer)
- [Download ImageMigick here.](https://imagemagick.org/script/download.php#windows)
- [Download source tex file here.](https://github.com/iitrabhi/paper-template)

# Install

- Install Sumatra first and then Miktex and then ImageMagick. Finally install Sublime.

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

- Save the file. This will automatically install all the packages necessary for the setup.  Wait for 5-10 mins for the installation to complete.

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

- Go to `Preferences→Key Bindings` and paste the following 

  ```json
  [
      { "keys": ["f1"], "command": "toggle_side_bar" },
      { "keys": ["f2"], "command": "distraction_free_window" },
      { "keys": ["f3"], "command": "fold" },
      { "keys": ["f4"], "command": "unfold" },
  ]
  ```

- Now finally go to `Preferences→Package settings→Latex tools→Check system` to check whether everything is fine or not.

# Dark Mode

In Sumatra PDF go to `Settings→Advanced Options`. The settings will open in a new text document. Change the `MainWindowBackground = #252a33` and replace the code in `FixedPageUI` with the following

```
FixedPageUI [
	# Light Mode
	# TextColor = #000000
	# BackgroundColor = #ffffff

	# Dark Mode
	TextColor = #ffffff
	BackgroundColor = #343d46

	SelectionColor = #f5fc0c
	WindowMargin = 2 4 2 4
	PageSpacing = 4 4
]
```

This will activate dark mode in Sumatra. To revert back to light mode just un-comment the lines under `Light Mode` and comment the lines under `Dark Mode`.

```
FixedPageUI [
	# Light Mode
	TextColor = #000000
	BackgroundColor = #ffffff

	# Dark Mode
	# TextColor = #ffffff
	# BackgroundColor = #343d46

	SelectionColor = #f5fc0c
	WindowMargin = 2 4 2 4
	PageSpacing = 4 4
]
```

# Learn Latex

- [Overleaf](https://www.overleaf.com/learn/latex/Learn_LaTeX_in_30_minutes)
- [Latex Tutorial](https://latex-tutorial.com/quick-start/)