---
layout: page-fullwidth
title:  "What is SALR particle clustering?"
subheadline:  "Data clustering"
teaser: "SALR clustering »<em>short-range attractive long-range repulsive particle clustering</em>« is a <a href='https://en.wikipedia.org/wiki/Cluster_analysis'>data clustering</a> technique to locate the centers of partially overlapping convex objects or distributions, and it overcomes several common problems of clustering methods."
categories:
    - information
tags:
    - data clustering
    - seed-point detection
    - machine learning
    - salr clustering
    - short-range attractive long-range repulsive
permalink: /what-is-it/
header:
    image_fullwidth: page_header.svg
    caption: "Particles close to each other are attracted together and particles far from each other are repulsed apart."
---

<div class="row">
	<div class="medium-7 columns t30">
    <h2>Data clustering</h2>
    <blockquote>
    Is there structure to the data?<br>
    Is the data made from distinct/overlapping groups?
    </blockquote>
    <p>These are often the first questions asked when analyzing data, and data clustering helps to answer them by searching for regions in the data where a measure of the difference between the points in each group is minimized. An example of data clustering is shown in this image, where the center of seven different data clusters have been found. Each data points can then be assigned to the cluster that has the smallest difference with itself.</p>
	</div><!-- /.medium-7.columns -->

	<div class="medium-5 columns t30">
    <img src="{{ site.urlimg }}clustering1.png">
    <figcaption class="text-right">
    A distribution of data where the center of seven clusters has been found. This image is looking down on a 3D point distribution.
    </figcaption>
	</div><!-- /.medium-5.columns -->
</div><!-- /.row -->




<div class="row">
	<div class="medium-5 columns t30">
    <img src="{{ site.urlimg }}what_is_it.gif">
    <figcaption class="text-right">
    » Data distribution » confining potential (dark blue: hill, light green: valley) » model SALR particles »  particle clusters used for <a href="https://en.wikipedia.org/wiki/Voronoi_diagram">Voronoi clustering</a> »
    </figcaption>


	</div><!-- /.medium-7.columns -->

	<div class="medium-7 columns t30">
    <h2>SALR particle clustering</h2>
    <p>SALR particle clustering locates cluster centers by first creating an intermediate <q>surface</q> called the »confining potential«. The confining potential can be thought of as being like <abbr title="maxima of the confining potential">hills</abbr> and <abbr title="local minima of the confining potential">valleys</abbr>, with valleys near the locations we think the data clusters could be and hills near the locations where they could not be. The data clusters are then determined by modeling how a set of particles move in the confining potential and extracting the final positions of the particles. The key aspect of SALR clustering is the short-range attractive long-range repulsive interaction of the particles.</p>


	</div><!-- /.medium-5.columns -->
</div><!-- /.row -->

<div class="row">
	<div class="medium-7 columns t30">
    <h3>Why short-range attractive long-range repulsive interaction?</h3>

    <p>Consider a set of repulsive particles (e.g. electrons) confined to a region. These particles will try to move as far apart from each other as possible. Image <b>a</b> shows this for 10 particles confined to a region with one valley. If the region the particles are confined to has several valleys, then the particles will try to be near these valleys. Assuming our goal is to have particles at each valley and nowhere else, we must use the same number of particles as valleys (since two particles cannot be near each other), see images <b>b</b> and <b>c</b>. This is a problem as the number of valleys is normally not known before hand, but it can be solved by modifying how the particles interact with each other, so particles near each other are attracted to each other (not repulsed). Now, we can use (many) more particles than we expect to be valleys, and we can have a cluster of particles located at each valley, see image <b>d</b>.</p>
	</div><!-- /.medium-7.columns -->

	<div class="medium-5 columns t30">
    <img src="{{ site.urlimg }}theoryIntro.svg">
	</div><!-- /.medium-5.columns -->
</div><!-- /.row -->

## What are the advantages of SALR clustering?

There are three primary advantages of using SALR clustering over other clustering/seed-point detection methods.

* The number of clusters does not need to be known or guessed beforehand, only the approximate size of the clusters must be known.
* The cluster locations are not strongly biased to regions with very high data density, this makes it useful in cases with rare clusters (low density clusters next to high density clusters).
* The cluster locations do not need to have local maximum in the data density, as opposed to mean-shift clustering[^2], for example.
* SALR clustering can be very fast: one repetition on a 5D data set with 2.7 million points only takes about 2 seconds, as opposed to k-means clustering[^3] where one repetition takes about 15 seconds.[^1]
* SALR clustering achieves an F1 score that is 8.1% higher than any other method when used for locating the centers of overlapping nuclei, and SALR clustering is able to better locate rare clusters without density maxima than k-means clustering.[^1]

> SALR clustering can represent a significant improvement in locating the centers of overlapping convex objects: it locates the correct number of nuclei more often and the nuclei centers more accurately than standard and leading methods; it can significantly improve the performance of previous methods; and it is able to determine, not only the number of clusters, but the correct position of the cluster centers in data clustering while not required a cluster to have a local density maximum.
<cite>J. Kapaldo et al., **(submitted)**[^1]</cite>

*[SALR]: short-range attractive long-range repulsive
*[repulsive particles]: Particles that do not want to be near each other and push each other apart.


[^2]: D. Comaniciu and P. Meer, <a href='http://dx.doi.org/10.1109/34.1000236'>IEEE Transactions on Pattern Analysis and Machine Intelligence 24, 603 (2002).</a>
[^3]: J. MacQueen, <a href="http://www-m9.ma.tum.de/foswiki/pub/WS2010/CombOptSem/kMeans.pdf">Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability 1, 281 (1967).</a>
[^1]: J. Kapaldo, X. Han, and D. Mary. **(submitted)**



<div d="post-nav" class="row">
<div class="small-12 text-center columns">
<a class="button small radius next" href="{{ site.baseurl }}/theory/">The theory &raquo;</a>
</div>
</div>
<div d="post-nav" class="row">
<div class="small-12 text-center columns">
<a class="button small radius next" href="{{ site.baseurl }}/applications/">Applications &raquo;</a>
</div>
</div>
<div d="post-nav" class="row">
<div class="small-12 text-center columns">
<a class="button small radius next" href="{{ site.baseurl }}/quick-start/">Getting started &raquo;</a>
</div>
</div>
