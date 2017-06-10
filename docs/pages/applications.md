---
layout: code-example
title:  "Applications and validation"
subheadline:  "What can SALR clustering be used for?"
teaser: "Locating object and distribution centers is a key step in both image processing and unsupervised machine learning. Both of these disciplines are used in the automated analysis of biological images. Here, SALR particle clustering is used to located the centers of overlapping cell nuclei as well as to locate clusters in scatter point data describing nuclei damage."
categories:
    - information
tags:
    - data clustering
    - seed-point detection
    - machine learning
    - salr clustering
    - short-range attractive long-range repulsive
permalink: /applications/
header:
    image_fullwidth: page_header.svg
---


## Locating cell nuclei centers

<div class="row">
<div class="medium-8 columns t30" markdown="1">

When automatically analyzing a biological image it is necessary to separate the objects of interest (cells/nuclei) from the background; this is called &raquo;[segmenting](https://en.wikipedia.org/wiki/Image_segmentation)&laquo;. Some of the best methods of segmenting (repulsive level-set[^Qi] or active contours[^Zhang]) require &raquo;seed-points&laquo; before they can be used. A seed-point is a point known to be inside the object of interest. In order to obtain good segmentation results, it is important that there only be one seed-point per object located as close as possible to the object's center.

</div><!-- /.medium-7.columns -->
<div class="medium-4 columns t30">
<img src="/images/nuclei.png">
<figcaption class="text-right">
Example of overlapping nuclei. Red points represent the manually labeled centers of the nuclei.
</figcaption>
</div><!-- /.medium-5.columns -->
</div><!-- /.row -->


<div class="row">
<div class="medium-7 columns t30" markdown="1">

### Results

We applied SALR clustering to a set of ~2500 images of nuclei clumps (examples shown above) and compared the results with those of several of the leading and popular methods:  MPV[^Parvin], SPV<sub>Qi</sub>[^Qi], SPV<sub>Xu</sub>[^Xu_spv], SBF[^Quelhas], gLoG[^Xu_glog], MSER[^Matas], the distance transform local maximum.  The comparison metric used is the [F<sub>1</sub> score](https://en.wikipedia.org/wiki/F1_score), F<sub>1</sub> = 2TP / (2TP + FN + FP), where TP, FN, and FP are the number of true positives, false negatives, and false positives.

The F<sub>1</sub> score (F<sub>1</sub>=1: good, F<sub>1</sub>=0: bad) depends on how close a calculated seed-point must be to a true nuclei center before it is considered a TP point. For this reason, the F<sub>1</sub> score is computed over a range, starting with the requirement that the seed-point and true center be very close to each other (&delta;r=3 pixels), and going to the requirement that they need not be very close to each other (&delta;r=10 pixels). _Note, the average nuclei radius is 19 pixels._

The results of the comparison are in the figure above; SALR clustering is 8.2% better than the next best method at &delta;r=3, and 3.0% better than the next best average F<sub>1</sub> score (across &delta;r).

</div><!-- /.medium-7.columns -->
<div class="medium-5 columns t90">
<img src="/images/resultsComparison.svg">
<figcaption class="text-right">
Comparison of F<sub>1</sub> versus &delta;r. The truth data used in computing F<sub>1</sub> was created by manually labeling the center of each nuclei in the images (~7800 nuclei).
</figcaption>
</div><!-- /.medium-5.columns -->
</div><!-- /.row -->


<div class="row">
<div class="medium-7 medium-push-5 columns t30" markdown="1">

### Method

The steps of applying SALR clustering to locating the centers of cell nuclei are briefly described below and a graphical animation of the process can be seen to the right. _For details, please see the manuscript[^1]._

**1.** An initial segmentation of the image is made with the goal of creating a binary mask of the objects, but not trying to separate the overlapping objects.[^BW]

**2.** The confining potential is created using V(r)=1/dt(~BW(r)), where dt() is the distance transform and ~ represents the not operator. This confining potential can be improved to be scale invariant by first scaling the distance transform to a range of [1, &lambda;<sub>max</sub>] before computing the confining potential.

**3.** The particles are initialized. There are two important quantities worth mentioning, the particle charge _q_ and the particle's initial positions. The charge is set to _q_=_N_<sup>-&beta;</sup>, where _N_ is the number of particle used in the simulation and &beta;=1/3; this helps to normalize the interaction strength of large particle clusters. There are several methods of initializing the particle positions. The method we use first computes a set of possible positions using centers of curvature of each boundary vertex from the binary mask. We then overlay a hexagonal lattice over the binary mask and select one of the possible positions from each hexagon. The size of the hexagons, defined by the particles Wigner-Seitz radius _r_<sub>s</sub>, controls the number of particles in the simulation.

**4.** The particles dynamics are modeled and the centers of each cluster of particles are used as the final seed-points.

</div><!-- /.medium-7.columns -->
<div class="medium-5 medium-pull-7 columns t90">
<img src="/images/animated2d.gif">
</div><!-- /.medium-5.columns -->
</div><!-- /.row -->



## Data clustering

<div class="row">
<div class="medium-8 columns t30" markdown="1">

In unsupervised machine learning and data mining, data often comes in the form of scatter point data where each point is an observation and each dimension is a different measured quantity. The easiest way _(but not the only way)_ to apply SALR clustering to scatter point data is to bin the data: we overlay a grid on the data and then count the number of points in each bin (grid element). The result of this can be thought of as a (multi-dimensional) image.

</div><!-- /.medium-7.columns -->
<div class="medium-4 columns t30">
<img src="/images/data_binning.svg">
<figcaption class="text-right">
Example of binning 2D data.
</figcaption>
</div><!-- /.medium-5.columns -->
</div><!-- /.row -->

<div class="row">
<div class="medium-7 medium-push-5 columns t30" markdown="1">

### Results

We test our method on a 5D data set that represents healthy and damaged cells. To the right we show the count isosurfaces (the count is directly proportional to the data point density) of a projection of the data to 3d. You can see the data set has one dense region (which are the healthy cells) with three primary low density regions extending from it. Our method is able to locate all of the regions well, but k-means++ is not.

### Method

If the different locally convex regions of the image are approximately the same size along any direction (or they can be made to have the same size with a simple scaling of the axis), then the distance transform can again be used to create the confining potential. This could be the case for 3d images of cell nuclei; however, in general, this requirement is quite stringent and the distance transform cannot be used.

Instead, one over the point density (the value of each pixel in our binned data) can be used as the confining potential. However, it is very important that the magnitude of the potential gradient \|∇V\| be scaled so that it is approximately of the same order of magnitude as the particle repulsion 1/rₐ².

</div><!-- /.medium-7.columns -->
<div class="medium-5 medium-pull-7 columns t90">
<img src="/images/animated5d.gif">
</div><!-- /.medium-5.columns -->
</div><!-- /.row -->



[^1]: J. Kapaldo, X. Han, D. Mary, and S. Ptasinska. Nature Methods.
[^Qi]:  X. Qi, F. Xing, D. J. Foran, and L. Yang, [IEEE Transactions on Biomedical Engineering 59, 754 (2012)](http://dx.doi.org/ 10.1109/TBME.2011.2179298).
[^Zhang]: X. Zhang, H. Su, L. Yang, and S. Zhang, [Proceedings of the IEEE Computer Society Conference on Computer Vision and Pattern Recognition 07-12-June, 5361 (2015)](http://dx.doi.org/ 10.1109/CVPR.2015.7299174).
[^Parvin]: B. Parvin, Q. Yang, J. Han, H. Chang, B. Rydberg, and M. H. Barcellos-Hoff, [IEEE Transactions on Image Processing 16, 615 (2007)](http://dx.doi.org/ 10.1109/TIP.2007.891154).
[^Xu_spv]:  H. Xu, C. Lu, and M. Mandal, [IEEE Journal of Biomedical and Health Informatics 18, 1729 (2014)](http://dx.doi.org/10.1109/JBHI.2013.2297030).
[^Quelhas]:  P. Quelhas, M. Marcuzzo, A. M. Mendon¸ ca, and A. Campilho, [IEEE Transactions on Medical Imaging 29, 1463 (2010)](http://dx.doi.org/10.1109/TMI.2010.2048253).
[^Xu_glog]: H. Xu, C. Lu, R. Berendt, N. Jha, M. Mandal, and S. Member, [Journal of Biomedical and Health Informatics XX, 1 (2016)](http://dx.doi.org/10.1109/JBHI.2016.2544245).
[^Matas]: J. Matas, O. Chum, M. Urban, and T. Pajdla, [In British Machine Vision Conference , 384 (2002)](http://dx.doi.org/ 10.5244/C.16.36).
[^Otsu]: N. Otsu, IEEE transactions on systems, man, and cybernetics 9, 62 (1979).
[^BW]: SALR clustering does not handle this step as it can be very different depending on the type of image. However, the binary mask used in creating the results were generated using an adaptive Otsu threshold[^Otsu], and the manuscript[^1] shows that excellent results can also be obtained using the binary mask created by gLoG filtering followed by thresholding.

*[MPV]: Multi-pass voting
*[SPV]: Single-pass voting
*[gLoG]: Generalized Laplacian of Gaussian
*[MSER]: Maximally stable extremal regions
*[SBF]: Sliding-band filter
*[SALR]: Short-range attractive long-range repulsive
