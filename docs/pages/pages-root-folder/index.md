---
#
# Use the widgets beneath and the content will be
# inserted automagically in the webpage. To make
# this work, you have to use › layout: frontpage
#
layout: frontpage
header:
  image_fullwidth: page_header.svg
  caption: "The content on this site summarizes work presented in <strong>J. Kapaldo et al., (submitted)</strong>."
widget2:
  title: "Applications"
  url: '/applications/'
  image: applications_img_sm.gif
  text: '<em>SALR particle clustering</em> can be applied to any problem where the center of overlapping objects or distributions needs to be found. This covers a broad range of fields, from unsupervised maching learning and data clustering to locating the centers of overlapping nuclei in biological images.'
widget1:
  title: "What is SALR particle clustering?"
  url: '/what-is-it/'
  image: img1.png
  text: 'Data clustering is prevalent in almost all areas technology and science, from identifying documents that are similar to analyzing biological images. <em>SALR particle clustering</em> is a data clustering technique for locating the centers of partially overlapping objects or distributions.'
widget3:
  title: "Getting started"
  url: '/quick-start'
  image: widget-github-303x182.jpg
  text: '<em>SALR particle clustering</em> is written in Matlab. Setup is as simple as downloading the <a href="https://github.com/jkpld/SALR_Clustering">repository</a> and running the included setup functions. The code is fully documented with several examples to get you started using it quickly.'
#
# Use the call for action to show a button on the frontpage
#
# To make internal links, just use a permalink like this
# url: /getting-started/
#
# To style the button in different colors, use no value
# to use the main color or success, alert or secondary.
# To change colors see sass/_01_settings_colors.scss
#
callforaction:
  url: '/images/manuscript.pdf'
  text: View the original paper ›
  style: alert
permalink: /index.html
#
# This is a nasty hack to make the navigation highlight
# this page as active in the topbar navigation
#
homepage: true
---
