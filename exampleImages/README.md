# Data Description
This document provides a description of the data files included as well as the format of the truth data provided.

## Images
The images provided with this work are composite images formed by combining many nuclei clumps from a full-slide image. The image names follow the following convention

`testImage_(image|mask)_LD(clumpSize)P24[_normalized].tif`

where
* The `image` files contain the actual images of the nuclei clumps without modification.
* The `mask` files contain the mask of each nuclei clump as obtained using an adaptive log-weighted Otsu thresholding on the full-slide image.
* `clumpSize` is one of the following values, `2`, `3`, `4`, `5`, or `67`.
* The optional flag `_normalized` at the end gives the images of the nuclei clumps after intensity normalization. _These are the images used with the previous methods we compared to in the manuscript._

## Truth data
