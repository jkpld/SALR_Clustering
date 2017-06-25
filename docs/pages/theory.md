---
layout: page-fullwidth
title:  "The Theory"
subheadline:  "How SALR particle clustering works"
teaser: "Motivation of SALR clustering comes from the fundamental physics of modeling classical Wigner crystals<a href='http://dx.doi.org/10.1006/spmi.1993.1026'><abbr title='F. Bolton and U. Rössler, Superlattices and Microstructures 13, 139 (1993).'><sup>1</sup></abbr></a> and modeling the formation of clusters<a href='http://dx.doi.org/10.1021/la048554t'><abbr title='S. Mossa, F. Sciortino, P. Tartaglia, and E. Zaccarelli, Langmuir 20, 10756 (2004).'><sup>2</sup></abbr></a><sup>,</sup><a href='http://dx.doi.org/10.1103/PhysRevLett.93.055701'><abbr title='F. Sciortino, S. Mossa, E. Zaccarelli, and P. Tartaglia, Physical Review Letters 93, 5 (2004).'><sup>3</sup></abbr></a>. This page goes over the theoretical model driving the particle dynamics, and then discusses how/why/when SALR clustering is able to perform better than other clustering/seed-point detection methods."
header:
    image_fullwidth: page_header.svg
categories:
    - information
tags:
    - theory
    - data clustering
    - seed-point detection
    - salr clustering
    - short-range attractive long-range repulsive
permalink: /theory/
---

## How the particles are modeled

<div class="row">
<div class="medium-7 columns t30" markdown="1">

A set of N interacting charge particles can be described by the [Hamiltonian](1) (1), where _m_, **p**, **r**, and _q_ represent the mass, momentum, position, and charge of a particle, and \|**r**<sub>i</sub> - **r**<sub>j</sub>\| is the distance between **r**<sub>i</sub> and **r**<sub>j</sub>. _V_ is the confining potential, and _V_<sub>int</sub> is the short-range attractive long-range repulsive (SALR) interaction potential created with a Gaussian attractive term and a _1/r_ repulsive term, given by (2). An example of what this potential looks like is shown in the figure. The values _A_, _µ_, and _σ_ in (2) are set implicitly by giving the potential depth, _d_<sub>o</sub>; the location of the potential minimum, _r_<sub>o</sub>; and the location at which the potential switches from being attractive to repulsive, referred to as the attractive extent, _r_<sub>a</sub>.

</div>
<div class="medium-5 columns t30">
<img src="{{ site.urlimg }}eq1.png">
<img src="{{ site.urlimg }}eq2.png">
<img src="{{ site.urlimg }}eq34.png">
</div>
</div>

<div class="row">
<div class="medium-5 columns t30">
<img src="{{ site.urlimg }}interactionPotential.png">
</div>

<div class="medium-7 columns t30" markdown="1">

The equations actually describing the motion of the particles can be found from the Hamiltonian, and are given by (3) and (4). And extra term has been added to (3), -&alpha;(t) / _m_<sub>i</sub> · **p**<sub>i</sub>, that does not come from the Hamiltonian; this term slows the particles down over time so that they eventually come to a stop.

Equations (3) and (4) give a set of 2N differential equations that are modeled until the solution converges. The solution is said to converge when the average particle speed over the last few solution steps is below a pre-defined threshold.

</div>
</div><!-- /.row -->



## How, why, and when SALR clustering works?

The answer comes by analyzing the force on the particles in the simulation: the confining potential puts a force on the particles with a magnitude of |∇_V_| (which is the slope of the confining potential, like a ball rolling down a hill), and the particles repel each other with a force of about _r_<sub>a</sub><sup>-2</sup>. If the force from the confining potential is much larger than the repulsive force, |∇_V_|≫_r_<sub>a</sub><sup>-2</sup>, then the particles will only find the local minimum of the confining potential---this is what most/all other methods do in finding the local extrema of a voting landscape or finding the local maxima of the point density (mean-shift[^5]). Alternatively, if the repulsive force is much larger than the confining force, _r_<sub>a</sub><sup>-2</sup>≫|∇_V_|, then the particles will ignore the confining potential minima and spread out as far apart from each other as possible---this will result is clusters of particles, separated by at least a distance _r_<sub>a</sub>, that uniformly cover the region. In between these two regimes, when the two forces are approximately the same order of magnitude, |∇_V_|&#126;_r_<sub>a</sub><sup>-2</sup>, there is an interaction between the confining potential and the particle repulsion, and this is the regime that leads to SALR clustering's improved performance.
<cite>J. Kapaldo et al., **(submitted)**[^1]</cite>



#### For additional details, please consult the publication [^1]

{% include next-prev next-text='Applications' next-url='/applications/' prev-text='What is SALR clustering' prev-url='/what-is-it/' %}


[1]: https://en.wikipedia.org/wiki/Hamiltonian_mechanics

[^1]: J. Kapaldo, X. Han, and D. Mary. **(submitted)**
[^2]: F. Bolton and U. Rössler, [Superlattices and Microstructures 13, 139 (1993)](http://dx.doi.org/10.1006/spmi.1993.1026).
[^3]: S. Mossa, F. Sciortino, P. Tartaglia, and E. Zaccarelli, [Langmuir 20, 10756 (2004)](http://dx.doi.org/10.1021/la048554t).
[^4]: F. Sciortino, S. Mossa, E. Zaccarelli, and P. Tartaglia, [Physical Review Letters 93, 5 (2004)](http://dx.doi.org/10.1103/PhysRevLett.93.055701).
[^5]: D. Comaniciu and P. Meer, [IEEE Transactions on Pattern Analysis and Machine Intelligence 24, 603 (2002)](http://dx.doi.org/10.1109/34.1000236).

*[SALR]: short-range attractive long-range repulsive
