Camille Couprie, 11 December 2012

/* ------------------ GENERAL INFORMATION ----------------------*/

Implementation of the segmentation algorithm described in:
 "Causal graph-based video segmentation", 
C. Couprie, C. Farabet, Y. Lecun, 2013

The program takes a sequence of images extracted from a video and
produces a temporally consistent segmentation with a random color color
assigned at each region.


/* ------------- SOFTWARES USED IN THIS PACKAGE --------------- */

The files convolve.h, filter.h, image.h, imutil.h, pnmfile.h,
 imconv.h, disjoint-set.h, segment-image.h, segment-graph.h, misc.h of
 this package come from the segmentation package implementing

Efficient Graph-Based Image Segmentation
Pedro F. Felzenszwalb and Daniel P. Huttenlocher
International Journal of Computer Vision, 59(2) September 2004.

available at http://www.cs.brown.edu/~pff/segment/

The code in the files basic-morpho.h, list_utils.h is extracted and modified code from 
the Pink image library available at http://pinkhq.com/.

The rest of the files (scene-labeling.cpp, semantic-segment.h, graph-matching.h,
 label-components.h, segment-video.h) contain specific code for the
 video segmentation described in "Causal graph-based video
 segmentation", C. Couprie, C. Farabet, Y. Lecun, 2013

/* -------------------- USAGE INSTRUCTIONS --------------------- */

The program takes a sequence of color images (PPM format) and produces a segmentation
with a random color assigned to each region.

1) Type "make" to compile "segment".

2) Run "segment k min_size sigma input (without .ppm) output (without .ppm) nbimages nb_classes static_mode <OPT labeling (without .pgm)>".

The parameters are: (see the paper for details)

* sigma: Used to smooth the input image before segmenting it.  
* k: Value for the threshold function. 
* min_size: Minimum component size enforced by post-processing.
* input: prefix name for the sequence of input image. We
 assume that the images are ordered and that their name ends
 with 5 digits and the extension ".ppm" Example: if the 
 sequence of images to segment are named "img-00001.ppm, 
 img-00002.ppm, img-00003.ppm", the input argument is "img-"
* output: Output prefix name to store the resulting image
 sequence. 
* nbimages: number of images in the sequence
* nb_classes: 0 if no semantic segmentation, otherwise, the 
 number of possible classes in the final semantic segmentation.
* static_mode: 0 if the camera is moving, 1 otherwise
* optional: if nb_classes>0, prefix name for the sequence of 
 input semantic predictions, computed image per image, to be 
 temporally smoothed   
* optional: maximum displacement parameter 

Typical parameters are sigma = 0.5, k = 400, min = 100.
Larger values for k result in larger components in the result.


/*--------------------- DEMO EXAMPLES ---------------------------*/


State of the art of semantic segmentation on NYU-Scenes dataset: 
(for evaluation, uncomment #define GROUND_TRUTH_AVAILABLE in misc.h)
./segment 1200 100 1.2  videos/NYU-Scenes/input- out 73 34 0 videos/NYU-Scenes/argmax- 70 

Segmentation on NYU-Scenes dataset: 
time ./segment 350 50 0.7  videos/NYU-Scenes/input- out 73 0 0  

360 video:
./segment 400 50 0.9  videos/NYC-360deg/input- out 50 34 0 videos/NYC-360deg/argmax- 15

Two women:
./segment 130 480 0.5 videos/two_women/two_women- out 90 0 1 



/* ------------ USEFUL PROGRAMS FOR VIDEO PROCESSING -----------*/

convert a video to images 
ffmpeg -i video.mp4 -s 320x240 image%05d.png 

create movies from a list of images : 
ffmpeg -r 5 -i out-%05d.png -vcodec libx264 movie.mp4


/* -------------------------- DATA ---------------------------- */

The folder "video" contains 3 video sequences: 
- the NYU-Scenes dataset was acquired by Clement Farabet 
- the NYC-360deg dataset was acquired by Marco Scoffier 
For both dataset described above, semantic predictions are given 
in the files "argmax-*.pgm" computed frame by frame using the scene parsing technique described in 
Clement Farabet, Camille Couprie, Laurent Najman, Yann LeCun 
"Learning Hierarchical Features for Scene Labeling" 
To appear in IEEE Transaction on Pattern Analysis and Machine Intelligence, 2013.

- the sequence two-women is extracted from Sylvain Paris' video
http://people.csail.mit.edu/sparis/publi/2008/eccv/Paris_08_Stream_Processing.mov 
