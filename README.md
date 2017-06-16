


# SALR particle clustering
### A method for locating the centers of partially overlapping convex objects and distributions.
_This project is apart of the manuscript **J. Kapaldo et al. Nature Methods (submitted)**._


SALR particle clustering is a simple technique for finding the centers of objects. It can be used to locate the centers of overlapping cell nuclei (see to the right), and it can be used to locate the cluster centers in scatter point data. Further details with examples can be seed at **website** and the manuscript.

<table style="width:100%">
  <tr>
    <td valign="top" width="45%" align="center"><img src="docs/images/animated2d.gif"></td>
    <td valign="top" width="45.8%" align="center"><img src="docs/images/animated5d.gif"></td>
  </tr>
</table>

## Installation

1. Download the repository.
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

## Example use

*[SALR]: short-range attractive long-range repulsive
