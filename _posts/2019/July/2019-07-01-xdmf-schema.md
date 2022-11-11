---
layout: post
title: "FEniCS: The XDMF schema"
description: "Understanding the XDMF schema as a part of GSoC 2019"
categories: [gsoc,coding]
typora-root-url: ../../../../website
---

XDMF stands for eXtensible Data Model and Format. This type of file format is used to store and exchange a large amount of scientific data between high-performance computing codes. FEniCS supports ‘io’ operations with XDMF file format which is further supported by the visualization tool of choice PARAVIEW. In this post, I would describe a bit about the XDMF and HDF5 file format.

### HDF5

HDF5 stands for Hierarchical Data Format version 5. This type of file format is used to store large volumes of data in a structured format so that we could easily retrieve the data when required. For a much more detailed overview of HDF5, the reader is referred to [Python and HDF5 by Andrew Collette](https://www.oreilly.com/library/view/python-and-hdf5/9781491944981/ch01.html). To view tables stored in an HDF file we use [HDFView](https://www.hdfgroup.org/downloads/hdfview/). You can think of HDF5 file as an excel spreadsheet with multiple tables. You can access the tables by their name and they have a hierarchical structure so that we could group them according to their use. The idea would get much clear with the discussion in the next section.

### XDMF

XDMF is based on XML file structure and thus the file format is highly readable by humans. XDMF file format is used to describe the data model of the data stored in the HDF5 file. This makes it easy for the end-user to read and understand the hierarchy of data stored in an HDF5 file. In addition to HDF5, XDMF can also store data in XML itself. Thus it is advisable to store heavy data (GB or TB) in XDMF and light data (KB or MB) in XML.

The code below presents a simple XDMF file where the data is stored in an HFD5 file named `poisson.h5`

```xml
<Xdmf Version="3.0">
  <Domain>
    <Grid Name="mesh" GridType="Uniform">
      <Topology NumberOfElements="2048" TopologyType="Triangle" NodesPerElement="3">
        <DataItem Dimensions="2048 3" NumberType="UInt" Format="HDF">poisson.h5:/Mesh/mesh/topology</DataItem>
      </Topology>
      <Geometry GeometryType="XY">
        <DataItem Dimensions="1089 2" Format="HDF">poisson.h5:/Mesh/mesh/geometry</DataItem>
      </Geometry>
      <Attribute Name="u" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="1089 1" Format="HDF">poisson.h5:/VisualisationVector/0</DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
```

### Understanding the code

- `<Xdmf>`: This is the root node of the tree. Every XML file needs a root node under which all the other nodes are housed.
- `<Domain>`: [Computational domain] could be one or more.
- `<Grid>`: Domains can have one or more `<Grid>`
- `<Geometry>`: Specify the location of grid nodes.
- `<Topology>`: Specify connectivity between nodes defined by `<Geometry>`.
- `<Attribute>`: It is used to specify values such as scalar or vectors, that are located at the node, edge, face, cell, center, or grid center.
- `<DataItem>`: It is used to specify the actual values for `<Geometry>`, `<Topology>`, or `<Attribute>`. The values are usually stored in an HDF5 file. For example, the values defining the location of nodes is stored in the file named poisson.h at location `/Mesh/mesh/geometry`.

You can read more about the various commands and parameters used in XDMF [here](http://www.xdmf.org/index.php/XDMF_Model_and_Format).
