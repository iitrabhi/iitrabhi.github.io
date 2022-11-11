---
layout: post
title: "How to have a clean latex folder"
categories: [latex, presentation]
typora-root-url: ../../../website
---

I am using sublime as the latex editor of choice, and I like it. You can learn more about it [here](https://abhigupta.io/2021/05/14/clean-sublime-setup.html). The only problem with the current workflow is that it creates all the output files and auxiliary files in the root directory. These additional files make a mess. I desire to keep things clean and minimal. 

Today I found a solution for my problem. We need to set the `aux_directory` and the `output_directory` keys in the `.sublime-project` file. Sublime will then create all the additional folders and files in that directory. If you need to clean your LaTeX project, you only need to delete the cache folder.

```json
{
   "folders":[
      {
         "file_exclude_patterns":[
            "*.cls",
            "*.md",
            "*.gitignore",
            "*.sublime-project"
         ],
         "path":"."
      }
   ],
   "settings":{
      "TEXroot":"thesis.tex",
      "ensure_newline_at_eof_on_save":true,
      "spell_check":true,
      "trim_trailing_white_space_on_save":true,
      "word_wrap":true,
      "wrap_width":100,
      "aux_directory":"cache",
      "output_directory":"cache"
   }
}
```

This procedure only works with TexLive in windows. On MAC, it works with MacTex. Unfortunately, it does not work with MikeTex.
