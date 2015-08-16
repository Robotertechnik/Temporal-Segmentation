/*
Copyright (C) 2006 Pedro Felzenszwalb

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





int main(int argc, char **argv) {
  if (argc != 5) {
    fprintf(stderr, "usage: %s input (without .ppm) output (without .ppm) nbimages ratio \n", argv[0]);
    return 1;
  }
  // (1) variables declarations
  char * imname = new char[100];
  char * outname = new char[100];
  char * appel = new char[1000];

 // (2) reading arguments
  int nb_images = atoi(argv[3]);
  int start=1;
  int i;
  sprintf(imname, "%s0%004d.ppm", argv[1], start);
  sprintf(outname, "%s0%004d.ppm", argv[2], start);
  float ratio = atof(argv[4]);  
  printf("%s\n",imname);
  image<rgb> *input;  
  // (3) loading first image
  printf("loading input image.\n");
  input = loadPPM(imname);
  int width = input->width();
  int height = input->height();
  int N = width*height;
  
  
  
  image<rgb> *output= new image<rgb>((int)(width/ratio), (int)(height/ratio));
  
 for (i=start;i<=nb_images;i++) {
 
   sprintf(imname, "%s0%004d.ppm",argv[1], i); 
   sprintf(outname, "%s0%004d.ppm",argv[2], i);   
   input = loadPPM(imname);
   
   int moy_r, moy_g, moy_b, norm;
   for (int y = 0; y < height/2; y++) {
     for (int x = 0; x < width/2; x++) {
      moy_r = imRef(input, x*2,y*2).r;
      moy_g = imRef(input, x*2,y*2).g;
      moy_b = imRef(input, x*2,y*2).b;
      norm = 1;                      
      if (x*2+1<width) {
	moy_r += imRef(input, x*2+1,y*2).r;
	moy_g += imRef(input, x*2+1,y*2).g;
	moy_b += imRef(input, x*2+1,y*2).b;
	norm++;}                                  
      if (y*2+1<height) {
	moy_r += imRef(input, x*2,y*2+1).r; 
	moy_g += imRef(input, x*2,y*2+1).g; 
	moy_b += imRef(input, x*2,y*2+1).b; 
	norm++;}
      if ((x*2+1<width) && (y*2+1<height)) {
	moy_r +=  imRef(input, x*2+1,y*2+1).r;
	moy_g +=  imRef(input, x*2+1,y*2+1).g;
	moy_b +=  imRef(input, x*2+1,y*2+1).b;
	norm++; }
      imRef(output, x,y).r = moy_r/norm;
      imRef(output, x,y).g = moy_g/norm;
      imRef(output, x,y).b = moy_b/norm;
    }
  }
  
  savePPM(output, outname);

  }

  

  printf("done.\n");
  
  return 0;
}

 
