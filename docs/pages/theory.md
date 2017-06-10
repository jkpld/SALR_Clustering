---
layout: page-fullwidth
title:  "The Theory"
subheadline:  "How SALR particle clustering works"
teaser: "Motivation of SALR clustering comes from the fundamental physics of modeling classical Wigner crystals<a href='http://dx.doi.org/10.1006/spmi.1993.1026'><abbr title='F. Bolton and U. Rössler, Superlattices and Microstructures 13, 139 (1993).'><sup>1</sup></abbr></a> and modeling the formation of clusters<a href='http://dx.doi.org/10.1021/la048554t'><abbr title='S. Mossa, F. Sciortino, P. Tartaglia, and E. Zaccarelli, Langmuir 20, 10756 (2004).'><sup>2</sup></abbr></a><sup>,</sup><a href='http://dx.doi.org/10.1103/PhysRevLett.93.055701'><abbr title='F. Sciortino, S. Mossa, E. Zaccarelli, and P. Tartaglia, Physical Review Letters 93, 5 (2004).'><sup>3</sup></abbr></a>. This page goes over the theoretical model driving the particle dynamics, and then discusses how/why/when SALR clustering is able to perform better than other clustering/seed-point detection methods."
header:
    image_fullwidth: syd-wachs-120737sm.jpg
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
	<div class="medium-7 columns t30">
    <p>A set of N interacting charge particles can be described by the <a href="https://en.wikipedia.org/wiki/Hamiltonian_mechanics">Hamiltonian</a> (1), where <em>m</em>, <b>p</b>, <b>r</b>, and <em>q</em> represent the mass, momentum, position, and charge of a particle, and |<b>r</b><sub>i</sub> - <b>r</b><sub>j</sub>| is the distance between <b>r</b><sub>i</sub> and <b>r</b><sub>j</sub>. <em>V</em> is the confining potential, and <em>V</em><sub>int</sub> is the short-range attractive long-range repulsive (SALR) interaction potential created with a Gaussian attractive term and a <em>1/r</em> repulsive term, given by (2). An example of what this potential looks like is shown in the figure. The values <em>A</em>, <em>µ</em>, and <em>σ</em> in (2) are set implicitly by giving the potential depth, <em>d<sub>o</sub></em>; the location of the potential minimum, <em>r<sub>o</sub></em>; and the location at which the potential switches from being attractive to repulsive, referred to as the attractive extent, <em>r<sub>a</sub></em>.</p>
    </div>

    <div class="medium-5 columns t30">
    <img src="/images/eq1.png">
    <img src="/images/eq2.png">
    <img src="/images/eq34.png">
	</div>
</div>

<div class="row">
    <div class="medium-5 columns t30">
    <img src="/images/interactionPotential.png">
    </div>

	<div class="medium-7 columns t30">
    <p>The equations actually describing the motion of the particles can be found from the Hamiltonian, and are given by (3) and (4). And extra term has been added to (3), -α(t) / <em>m<sub>i</sub></em> · <b>p</b><sub>i</sub>, that does not come from the Hamiltonian; this term slows the particles down over time so that they eventually come to a stop.</p>

    <p>Equations (3) and (4) give a set of 2N differential equations that are modeled until the solution converges. The solution is said to converge when the average particle speed over the last few solution steps is below a pre-defined threshold.</p>
	</div>

</div><!-- /.row -->



## How, why, and when SALR clustering works?

The answer comes by analyzing the force on the particles in the simulation: the confining potential puts a force on the particles with a magnitude of |∇<em>V</em>| (which is the slope of the confining potential, like a ball rolling down a hill), and the particles repel each other with a force of about <em>r</em><sub>a</sub><sup>-2</sup>. If the force from the confining potential is much larger than the repulsive force, |∇<em>V</em>|≫<em>r</em><sub>a</sub><sup>-2</sup>, then the particles will only find the local minimum of the confining potential---this is the same thing that the previous methods do in finding the local extrema of the voting landscape or finding the local maxima of the point density (mean-shift). Alternatively, if the repulsive force is much larger than the confining force, <em>r</em><sub>a</sub><sup>-2</sup>≫|∇_V_|, then the particles will ignore the confining potential minima and spread out as far apart from each other as possible---this will result is clusters of particles, separated by at least a distance _r_<sub>a</sub>, that uniformly cover the region. In between these two regimes, when the two forces are approximately the same order of magnitude, |∇_V_|&#126;r<sub>a</sub><sup>-2</sup>, there is an interaction between the confining potential and the particle repulsion, and this is the regime that leads to SALR clustering's improved performance.
<cite>J. Kapaldo et al., Nature Methods[^1]</cite>


#### For additional details, please consult the publication [1]

<!-- [1]: J. Kapaldo, X. Han, D. Mary, and S. Ptasinska. Nature Methods. (submitted) -->

<br>
<br>
<html>_______</html>
[1] J. Kapaldo, X. Han, D. Mary, and S. Ptasinska. Nature Methods. (submitted)
