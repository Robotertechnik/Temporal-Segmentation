/*
Copyright (C) 2012 Camille Couprie

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#include <cstdio>
#include <cstdlib>
#include <image.h>
#include <misc.h>
#include <pnmfile.h>
#include "segment-image.h"
#include "segment-video.h"
#include "segment-graph.h"
#include "semantic-segment.h"
#include "graph-matching.h"
#include "color_flow.h"

#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
using namespace std;
#define FILESIZE 100

void Print_Confusion(int ** Confusion, int nb_class);

/* The program takes a sequence of images extracted from a video and
produces a temporally consistent segmentation with a random color color
assigned at each region. */
int main(int argc, char **argv) {

  if (argc < 8 ) {
    fprintf(stderr, "usage: %s  k min_size sigma input (without .ppm) output (without .ppm) nbimages nb_classes static_mode <OPT labeling (without .pgm)> <OPT max_dep> \n", argv[0]);
    return 1;
  }

  // (1) variables declarations

  char * imname = new char[FILESIZE];
  char * outname = new char[FILESIZE];
  char * outname2 = new char[FILESIZE];
  char * prediction_name = new char[FILESIZE];
  char * labeling_name = new char[FILESIZE];
  char * optical_name = new char[FILESIZE];
  char * groundtruth_name = new char[FILESIZE];
  char * appel = new char[1000];
  image<rgb> *input;
  image<rgb> *next_input;
  image<uchar> *predictions;
  int num_ccs, num_ccs2, i,j,jj; 
  int num_edges, nedges_bef_corr;  
  struct timeval tim; double t1,t2; 
  label_image l1, l2;
  graph G; 
  float ** weights;
  int start=1;
 
 // (2) reading arguments

  sprintf(prediction_name, "%s0%004d.pgm", argv[9], start);
  float k = atof(argv[1]);
  int min_size = atoi(argv[2]);
  float sigma = atof(argv[3]);
  sprintf(imname, "%s0%004d.ppm", argv[4], start);
  sprintf(groundtruth_name, "%s0%004d.txt", argv[4], start);
  sprintf(outname, "%s0%004d.ppm", argv[5], start);
  int nb_images = atoi(argv[6]);
  int nb_class = atoi(argv[7]);
  int static_mode = atoi(argv[8]);
  int max_dep;
  if (argc <= 10) max_dep = 0;
  else max_dep = atoi(argv[10]);
  printf("%s\n",imname);
  
  // (3) loading first image

  printf("loading input image.\n");
  input = loadPPM(imname);

  int width = input->width();
  int height = input->height();
  int N = width*height;
  edge *edges = new edge[N*6]();
  int *im2_sizes = new int[N]; 
  int *im1_sizes = new int[N]; 
  int **centroids1 = new int*[N]; 
  for (int y = 0; y < N; y++) centroids1[y] = new int[4];
  int **centroids2 = new int*[N];
  for (int y = 0; y < N; y++) centroids2[y] = new int[4];

  // centroid[][0]: coordinate x of centroid 
  // centroid[][1]: coordinate y of centroid
  // centroid[][2]: index centroid point
  // centroid[][3]: mean intensity of the region

  int num_vertices = 2*N;
  int *labels = new int[num_vertices];
  rgb *colors = new rgb[num_vertices];
  uchar *semantic = new uchar[num_vertices];
  int *labels1 = new int[N];
  int *labels2 = new int[N];
  image<uchar> *output = new image<uchar>(width, height);
  image<rgb> *seg = new image<rgb>(width, height);
  image<rgb> *seg2 = new image<rgb>(width, height);
  image<uchar> *semantic_result = new image<uchar>(width, height);

  // Set parameters
  int DIST_CENTROID_MIN_DIST = (30.0*N)/(320*240);
  int DIST_CENTROID_MAX_DIST = (75.0*N)/(320*240);
  int LARGE_SIZE = ((1250.0)*N)/(320*240);
  if (static_mode > 0) { 
    DIST_CENTROID_MIN_DIST = DIST_CENTROID_MIN_DIST/2.0;
    DIST_CENTROID_MAX_DIST = DIST_CENTROID_MAX_DIST/2.0;
    LARGE_SIZE = LARGE_SIZE/2.0;
  }

  // confusion matrix
  int ** Confusion = new int*[nb_class]();
  for (int y = 0; y < nb_class; y++) Confusion[y] = new int[nb_class]();
  for (int x = 0; x < nb_class; x++)
    for (int y = 0; y < nb_class; y++)
      Confusion[x][y]=0;
    
  int *optx = new int [N]();
  int *opty = new int [N]();

  #ifdef OPTICALFLOW
    image<rgb> *output_optical = new image<rgb>(width, height);
  #endif

  #ifdef CONVERT
    sprintf(appel, "convert %s0%004d.png %s  ;", argv[9], start, prediction_name);
    system(appel); 
  #endif
    if (nb_class > 0) 
     predictions = loadPGM(prediction_name);
 
 
  // (4) segment first image and returns the labeled nodes

  segment_image(input, sigma, k, min_size, &num_ccs, &edges, &labels1, &im1_sizes, &seg, &centroids1, predictions, nb_class); 
  // fprintf(stderr,"image 1 : got %d components\n", num_ccs);
  #ifdef SAVE_ALL
    savePPM(seg, outname);
  #endif

  //(4 bis) semantic segmentation
  if (nb_class > 0) {
    for (int y = 0; y < height; y++) 
      for (int x = 0; x < width; x++) {
	imRef(semantic_result,x,y) = 0;
      }
    sprintf(labeling_name, "labeled%s0%004d.ppm", argv[5], start);
    semantic_segment( input, predictions, labels1, num_ccs, nb_class+1, labeling_name, &semantic_result, groundtruth_name, &Confusion);
  }
 
  // (5) treat subsequent images 
  for (i=start+1;i<=nb_images;i++) {
    
    // (6) loading image at time t
    sprintf(prediction_name, "%s0%004d.pgm", argv[9], i);
    sprintf(imname, "%s0%004d.ppm", argv[4], i); 
    sprintf(outname, "%s0%004d.ppm", argv[5], i); 
    next_input = loadPPM(imname);
    
    if (nb_class > 0) {
      sprintf(labeling_name, "labeled%s0%004d.ppm", argv[5], i);
      predictions = loadPGM(prediction_name);
  }

    // (7) segment the following image at time t independently from image t-1.
    segment_image(next_input, sigma, k, min_size, &num_ccs2, &edges, &labels2, &im2_sizes, &seg2, &centroids2, predictions, nb_class);
    // printf("image 2' : got %d components\n", num_ccs2);
   
    #ifdef SAVE_ALL
      sprintf(outname2, "alone_segm%05d.ppm",  i); 
      savePPM(seg2, outname2);
      sprintf(appel, "convert %s %s.png  ;", outname2, outname2);
      system(appel);
      sprintf(appel, "rm %s  ;", outname2);
      system(appel);
    #endif
    
    // (8) compute a similarity ratio between regions of image 1 and 2. The higher the closest regions are.
    l1.I = labels1;  l1.width = width;  l1.height = height; 
    l1.nb_cpts = num_ccs; l1.colorI= seg; l1.sizes = im1_sizes; l1.semanticI= semantic_result;
    l2.I = labels2; l2.width = width;  l2.height = height; 
    l2.nb_cpts = num_ccs2; l2.colorI= seg2; l2.sizes = im2_sizes;
    
    G =  build_graph(l1, l2, centroids1, centroids2, DIST_CENTROID_MIN_DIST, DIST_CENTROID_MAX_DIST, LARGE_SIZE);
     
    for (int x=0;x<num_vertices;x++)
      labels[x]=-1;

    graph_matching(G, l1, l2, centroids1, centroids2,  &edges, &num_edges, &nedges_bef_corr, labels, colors, semantic, static_mode);

    // (11) compute the final segmentation of image at time t 
    segment_video(next_input, min_size, edges, num_edges, nedges_bef_corr, &labels1, labels, colors, &num_ccs, num_ccs2, &im1_sizes, &seg, &centroids1, &optx, &opty, semantic, &semantic_result, max_dep, nb_class);


   // (12) Optical flow
    #ifdef OPTICALFLOW
      sprintf(optical_name, "optical%s0%004d.pgm", argv[5], i);
      MotionToColor(optx,  opty, &output_optical);
      savePPM(output_optical, optical_name);
      sprintf(appel, "convert %s %s.png  ;", optical_name, optical_name);
      system(appel);
      sprintf(appel, "rm %s  ;", optical_name);
      system(appel);
    #endif

    input = next_input;

   // (13) semantic segmentation
   
    if (nb_class > 0) {
      sprintf(groundtruth_name, "%s0%004d.txt",argv[4], i); 
      printf("image %02d ", i);
      semantic_segment( input, predictions, labels1, num_ccs, nb_class+1, labeling_name, &semantic_result, groundtruth_name, &Confusion);
       #ifdef GROUND_TRUTH_AVAILABLE
      //Print_Confusion(Confusion, nb_class);
       #endif
    }
 
    savePPM(seg, outname);

    #ifdef SAVE_ALL 
      sprintf(appel, "convert %s %s.png  ;", outname, outname);
      system(appel);
      sprintf(appel, "rm %s  ;", outname);
      system(appel);
      //printf("image %d : got %d components\n", i, num_ccs);
    #endif
  }

  //free memory
  for (int y = 0; y < N; y++) 
    delete [] centroids1[y];
  delete [] centroids1;
  for (int y = 0; y < N; y++) 
    delete [] centroids2[y];
  delete [] centroids2;
  delete [] labels; 
  delete [] colors;
  delete [] labels1;
  delete [] labels2;
  delete [] im2_sizes;
  delete [] im1_sizes;
  delete [] edges;
 
  printf("done.\n");
  
  return 0;
}



void Print_Confusion(int ** Confusion, int nb_class) {
  int nb_false_pos, nb_false_neg;
  for (int x = 0; x < nb_class; x++) {
    nb_false_pos = 0;
    for (int y = 0; y < nb_class; y++) { 
      printf("%4d ", Confusion[x][y]);
      nb_false_pos =	nb_false_pos + Confusion[x][y];
      nb_false_neg = 0;
      for (int xx = 0; xx < nb_class; xx++)
	nb_false_neg = nb_false_neg + Confusion[xx][x];
    }
    if ((nb_false_pos+nb_false_neg)==0) printf("|   0  \n");
    else printf("| %2.3f \n",(2.0*Confusion[x][x])/(nb_false_pos+nb_false_neg));
    }
}
