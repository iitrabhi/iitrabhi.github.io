## Welcome to my academic blog. 
This platform is being developed with Jekyll, which evolves with my progress.

Jekyll is a simple, blog-aware, static site generator perfect for personal, project, or organization sites. It takes text written in your favorite markup language and uses layouts to create a static website. You can tweak the site’s look and feel, URLs, the data displayed on the page, and more.

GitHub Pages is a service offered by GitHub that allows users to host websites directly from their GitHub repositories. It's particularly well-suited for Jekyll-based blogs or sites because GitHub Pages automatically builds and serves Jekyll sites, making the deployment process very efficient.

Below is a brief overview of the typical folders and files structure in a Jekyll site, organized in a table format:

| **Folder/File** | **Description** |
|-----------------|-----------------|
| `_posts` | Contains the blog posts written in Markdown. Each file is named with the date and title of the post. |
| `_layouts` | Contains templates that wrap around your content. These can be used to define the layout of your webpage. |
| `_includes` | Contains snippets of code that can be included in layouts and posts using Jekyll's `{% include file.ext %}` tag. |
| `_config.yml` | Stores configuration data for the Jekyll site. I define the title, description, URL, and other settings here. |
| `assets` | Typically used to store static files like images, JavaScript, and CSS used in the layouts and posts. |

To build a Jekyll website, several technologies are involved:

1. **Ruby**: Jekyll is written in Ruby, so you need Ruby installed to run Jekyll.
2. **Markdown/HTML**: Posts and pages can be written in Markdown or HTML, which Jekyll converts into static files.
3. **Liquid Templating Language**: Used for templates in Jekyll, allowing the integration of dynamic content in the site.
4. **YAML**: Used for configuration files and front matter, which is metadata for posts and pages.
5. **Sass/CSS**: Jekyll supports Sass out of the box for styling the site, and the stylesheets can be processed as part of the site build.

These technologies work together to efficiently create static sites that can be hosted anywhere, but they are particularly well-suited for hosting on GitHub Pages.

In Jekyll, the naming convention of folders and files is significant, especially with the use of underscores at the beginning of folder names. Here's the distinction:

1. **Folders with underscores:**
   - **Purpose:** These folders are special to Jekyll. They're used for processing and not included in the output site when built. They contain files that Jekyll processes to generate the static site. The contents are typically not served directly as they are; instead, they are combined, transformed, or read to generate the final site.
   - **Common Examples:**
       - `/_posts`: Contains the blog posts that are processed into the site's content.
       - `/_layouts`: Includes templates for pages and posts.
       - `/_includes`: Houses snippets of code to be included in layouts and posts.
       - `/_drafts`: Holds draft posts that aren't ready to be published.
       - `/_data`: Contains data files in .yml, .yaml, .json, or .csv format that can be used in the site.
   - **Behavior:** When Jekyll builds the site, it processes these directories specially. For example, files in `_posts` are turned into individual blog posts, and `_layouts` are used as templates for generating HTML files.

2. **Folders without underscores:**
   - **Purpose:** These are regular directories that Jekyll copies directly to the output site as-is. They're not processed in any special way but are included in the site generation. These folders typically contain assets that you want to be used in the final site as they are, like images, CSS, or JavaScript files.
   - **Common Examples:**
       - `/assets`: Usually contains static files like CSS, JavaScript, and images.
       - `/images`: A common folder for storing image files used in the site.
   - **Behavior:** The contents of these folders are typically referenced in your site's pages, posts, or layouts and are copied as-is into the `_site` directory (the default output directory of Jekyll) when the site is built.

In summary, folders prefixed with an underscore have a special meaning in Jekyll and are used for content that needs to be processed, while folders without an underscore are for content that should be copied as-is to the final site. 
## Building
GitHub Actions can automatically build and deploy your Jekyll site whenever you push changes to your repository. This means your live site updates automatically with each commit.

```bash
bundle exec jekyll serve --livereload --drafts
```

The command `bundle exec jekyll serve --livereload --drafts` is used for local development. It builds your site, serves it on your local machine, automatically refreshes the browser when changes are made (`--livereload`), and includes draft posts in the build (`--drafts`). This allows you to see and test changes in real-time as you develop your site.
## Categories
- coding
- simulation
- presentation
- productivity
- finance

## Guidelines

|Entity|Description|
|:--|:--|
|Image size|800 px wide|
|  |  |
|  |  |
|  |  |
