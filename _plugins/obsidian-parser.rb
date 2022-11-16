Jekyll::Hooks.register :posts, :pre_render do |post, payload|
  docExt = post.extname.tr('.', '')
  # only process if we deal with a markdown file
  if payload['site']['markdown_ext'].include? docExt
    content = post.content.gsub(/_posts\//, '{% post_url ')
    content = content.gsub(/\.md/, ' %}')
    post.content = content
  end
end










# module Jekyll
# class ObsidianParser < Converter
#   safe true
#   priority :highest

#   def matches(ext)
#     ext =~ /^\.(md|markdown)$/i
#   end
#   def output_ext(ext)
#     ".html"
#   end
#   def convert(content)
#     # content = content.gsub(/_posts\//, '% post_url ')
#     # content = content.gsub(/\.md/, ' %')

#     # Now call the standard Markdown converter
#     content
#     # site = Jekyll::Site.new(@config)
#     # mkconverter = site.find_converter_instance(Jekyll::Converters::Markdown)
#     # mkconverter.convert(content)
#   end
# end
# end

# https://stackoverflow.com/questions/35614552/with-jekyll-3-can-i-transform-a-posts-markdown-before-actual-markdown-parsing
# https://github.com/wildlyinaccurate/jekyll-responsive-image/issues/15