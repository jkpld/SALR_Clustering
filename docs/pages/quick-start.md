---
layout: page
sidebar: left
title:  "Quick start guide"
subheadline:  "Installation and requirements"
teaser: ""
categories:
    - information
tags:
    - documentation
    - help
permalink: /quick-start/
header:
    image_fullwidth: syd-wachs-120737sm.jpg
---

## Installation

1. Download the [GitHub repository](https://github.com/jkpld/seed_point_detection.git).
2. Open Matlab on your computer, and change the current folder to the downloaded repository.
3. Run `setup.m`. This function will
    * Add the necessary files to your Matlab path.
    * Attempt to compile two C functions that will increase the program speed. If a C compiler is not found on your computer, then the code will still run, but it could be slower.
4. Start by browsing the examples located in the `examples/` folder. These same examples can also be seed on this website >here<.


#### Installation Notes
* The GitHub repository could take a few minutes to download as the repository is ~200 MB. The majority of this space comes from the example data included.
* `setup.m` will not permanently add the files to your Matlab path. If you want the files to remain on your Matlab path, the run the command `savepath` in the Matlab command window after running `setup.m`.

## Requirements

* Matlab version R2016b or higher. _The reason for this requirement is that R2016b is the first version that allows for implicit expansion, which is a feature heavily used in the code._
* The Image Processing Toolbox. (?)
* The Statistics and Machine Learning Toolbox. (?)
* The Parallel Computing Toolbox - This toolbox is optional, but if you have it, then you will be able to run computations in parallel.
