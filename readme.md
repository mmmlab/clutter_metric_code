# Instructions on How to Use the Clutter Metric Code

*************
### Overview 

This repo contains a MATLAB implementation of the clutter metric described in Semizer and Michel (2018). 
This metric is a modified version of the clutter metric described in Bravo and Farid (2008).

*************

### How to run the code

This code uses a segmentation algorithm by Felzenszwalb and Huttenlocher (2004). The current implementation of the code requires a MATLAB wrapper for the original segmentation algorithm, which can be found [here](https://github.com/cvjena/Felzenszwalb-Segmentation). 

After downloading the wrapper, call `getClutter(image_dir)` function in **getClutter.m** to compute clutter as follows:

`[nr_regions,clutter_distro,clutter_ranks,seg_images,fileNames] = getClutter(image_dir)`

The input argument `image_dir` represents the directory which contains images. The function returns the following variables. 

Variable Name | Explanation
---------- | ----------------------------------------------------------------------------------------
`nr_regions` |  number of regions in the segmented images at each scale of segmentation
`clutter_distro` | distribution of clutter values across images
`clutter_ranks` | distribution of ranks of images based on clutter values
`seg_images` | example segmented images at single scale of segmentation
`fileNames` | names of image files

************* 

Copyright &copy; 2018 by Yelda Semizer <yelda.semizer@rutgers.edu> \& Melchi M. Michel <melchi.michel@rutgers.edu>
 
Licensed under GNU General Public License 3.0 or later. 

Some rights reserved. See COPYING, AUTHORS @license GPL-3.0+ <http://spdx.org/licenses/GPL-3.0+>
