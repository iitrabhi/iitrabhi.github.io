---
layout: post
title: "May 27th – May 31st, First week of coding"
description: "My first week as a GSoC2019 participant."
tag: 
  - [gsoc,coding]
typora-root-url: ../../../../website
---

The aim of the project is to preserve the named metadata in `.geo` files through to the tagged regions (MeshValueCollection) in dolfinx. This would help the end user to easily model with a large number of tagged regions. After a thorough discussion with my mentor Dr.Jack S. Hale and the author of “meshio”, [Dr. Nico Schlömer](https://github.com/nschloe), it was decided that the right way to proceed with the project is to use the `<Information>` tag to store a key value collection that contains the mapping between the user friendly strings and the numerical id.

This would require changes to be made in the source code of meshio and dolfinx. It was decided that I would initially work on dolfinx with a ‘mock-up’ of an XDMF file with the required metadata and then come back to meshio.

My plan was to initially start with the python version of dolfinx and then work on cpp version. Thus, I made a simple unit square in gmsh, exported it via meshio and then manually added the required data in the Information tag of XDMF. All the four boundaries were marked with string names and stored in information tag.

```xml
<Information Name="BoundaryTags">
       <![CDATA[
       <main>
       <map key="bottom">1</map>
       <map key="top">2</map>
       <map key="left">3</map>
       <map key="right">4</map>
       </main>
       ]]>
</Information>
```

The following code successfully read and created a dictionary of the map stored in the XDMF file.

```python
from lxml import etree

with open(r'Models/tag_all.xdmf') as fobj:
       xml = fobj.read()

root = etree.fromstring(xml)
for domain in root.getchildren():
       for elem in domain.getchildren():
           if elem.tag=='Information':
               cdata = etree.fromstring(elem.text)

tags = {}
for main in cdata.getchildren():
       tags[main.values()[0]]=int(main.text)
```

Output:

```
{‘bottom’: 1, ‘top’: 2, ‘left’: 3, ‘right’: 4}
```

In the next post, I will try to touch upon the concept of [docker](https://docs.docker.com/get-started/) and continuous integration (CI) using [CircleCI](https://circleci.com/continuous-integration/).
