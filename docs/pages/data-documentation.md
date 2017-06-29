---
layout: documentation
sidebar: left
title:  "Data documentation"
subheadline:  "How is the data stored?"
teaser: "This page provides a description of the data files included as well as the format of the truth data provided."
categories:
    - information
tags:
    - documentation
permalink: /data-documentation/
header:
    image_fullwidth: syd-wachs-120737sm.jpg
---

_The data files discussed here can be found in the [example data](https://github.com/jkpld/SALR_Clustering/tree/master/exampleData) folder of the [GitHub repository](https://github.com/jkpld/SALR_Clustering), or in the **Supplemental Data** of the manuscript[^1]._

## Images
The images provided with this work are composite images formed by combining many nuclei clumps from a full-slide image. The image names follow the following convention

`(image|mask)_(clumpSize)[_normalized].tif`

where
* The `image` files contain the actual images of the nuclei clumps without modification.
* The `mask` files contain the mask of each nuclei clump as obtained using an adaptive log-weighted Otsu thresholding on the full-slide image.
* `clumpSize` is one of the following values, `2`, `3`, `4`, `5`, or `67`.
* The optional flag `_normalized` at the end gives the images of the nuclei clumps after intensity normalization. _These are the images used with the previous methods we compared to in the manuscript._

Example:
<div class="row">
	<div class="medium-4 columns t30">
        <strong>image_2.tif</strong>
        <img src="{{ site.urlimg }}image_2.png">
    </div>
    <div class="medium-4 columns t30">
        <strong>mask_2.tif</strong>
        <img src="{{ site.urlimg }}mask_2.png">
    </div>
    <div class="medium-4 columns t30">
        <strong>image_2_normalized.tif</strong>
        <img src="{{ site.urlimg }}image_2_normalized.png">
    </div>
</div>

## Truth data
Each nucleus from the `testImage`s has been labeled. The data is stored in the five matlab files

`markedCenters_(clumpSize).mat`

where `clumpSize` is the same as above. Each `.mat` file contains a single `Nx3` array. Each row in the array corresponds to a labeled nuclei center. The first column gives the object number to which the nuclei belongs, the second column gives the y-coordinate (vertical axis), and the third column gives the x-coordinate (horizontal axis). The objects are numbered starting from in the upper left corner of the image and increasing down the columns. _This is the same numbering as obtained using Matlab's `bwlabel()`, or a similar function._

Here is an example of the object numbering and format of the truth data:

<div class="row">
<div class="medium-6 columns t30">
    <img src="{{ site.urlimg }}truth_example.png">
</div>
<div class="medium-6 columns t30" markdown="1">

| Obj # | y | x |
| :---- | :-: | :-: |
| 1 | 27 | 18 |
| 1 | 23 | 32 |
| 2 | 79 | 19 |
| 2 | 71 | 31 |
| ... | ... | ... |
| 8 | 81 | 120 |
| 8 | 69 | 130 |
| 9 | 127 | 132 |
| 9 | 123 | 118 |

</div>
</div>

## 5D data set
The 5-dimensional data set used in the manuscript is located in the Matlab file, `damage_data_5D.mat`. This file contains a single `Nx5` array where each row is a new point and each column is a dimension. This data file is approximately 100 MB.

[^1]: J. Kapaldo et al. **(submitted)**.
