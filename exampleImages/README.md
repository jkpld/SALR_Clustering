# Data Description
This document provides a description of the data files included as well as the format of the truth data provided.

## Images
The images provided with this work are composite images formed by combining many nuclei clumps from a full-slide image. The image names follow the following convention

`(image|mask)_(clumpSize)[_normalized].tif`

where
* The `image` files contain the actual images of the nuclei clumps without modification.
* The `mask` files contain the mask of each nuclei clump as obtained using an adaptive log-weighted Otsu thresholding on the full-slide image.
* `clumpSize` is one of the following values, `2`, `3`, `4`, `5`, or `67`.
* The optional flag `_normalized` at the end gives the images of the nuclei clumps after intensity normalization. _These are the images used with the previous methods we compared to in the manuscript._

Example:

| image_2.tif | mask_2.tif | image_2_normalized.tif |
| :-------------: | :--: | :--: |
| ![image](/docs/images/image_2.png) | ![mask](/docs/images/mask_2.png) | ![normalized image](/docs/images/image_2_normalized.png) |


## Truth data
Each nucleus from the `testImage`s has been labeled. The data is stored in the five matlab files

`markedCenters_(clumpSize).mat`

where `clumpSize` is the same as above. Each `.mat` file contains a single `Nx3` array. Each row in the array corresponds to a labeled nuclei center. The first column gives the object number to which the nuclei belongs, the second column gives the y-coordinate (vertical axis), and the third column gives the x-coordinate (horizontal axis). The objects are numbered starting from in the upper left corner of the image and increasing down the columns. _This is the same numbering as obtained using Matlab's `bwlabel()`, or a similar function._

Here is an example of the object numbering and format of the truth data:

![truth example](/docs/images/truth_example.png)

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

## 5D data set
The 5-dimensional data set used in the manuscript is located in the Matlab file, `damage_data_5D.mat`. This file contains a single `Nx5` array where each row is a new point and each column is a dimension. This data file is approximately 100 MB.
