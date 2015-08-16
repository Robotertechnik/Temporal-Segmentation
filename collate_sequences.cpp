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
  if (argc < 6) {
    fprintf(stderr, "usage: %s  nbimages ratio_horizontal ratio_vertical input1 (without .ppm) input2 (without .ppm) \n", argv[0]);
    return 1;
  }
  // (1) variables declarations
  char * imname = new char[100];
  char * imname2 = new char[100];
  char * outname = new char[100];
  char * appel = new char[1000];

 // (2) reading arguments
  int nb_images = atoi(argv[1]);
  int start=1;
  int i;
  sprintf(imname, "%s0%004d.ppm", argv[4], start);
  sprintf(imname2, "%s0%004d.ppm", argv[5], start);
  int ratio_horizontal = atoi(argv[2]);  
  int ratio_vertical = atoi(argv[3]);  
  printf("%s\n",imname);
  image<rgb> *input; image<rgb> *input2;   
  // (3) loading first image
  printf("loading input image.\n");
  input = loadPPM(imname);
  int width = input->width();
  int height = input->height();
  int N = width*height;
  
  
  image<rgb> *big = new image<rgb>(width*ratio_horizontal, height*ratio_vertical);
  
  
 for (i=start;i<=nb_images;i++) {
 
   sprintf(imname, "%s0%004d.ppm",argv[4], i); 
   sprintf(imname2, "%s0%004d.ppm",argv[5], i);   
   sprintf(outname, "%sout-0%004d.ppm", argv[5], i);
   input = loadPPM(imname);
   input2 = loadPPM(imname2);

   if(ratio_vertical==2){
      for (int y = 0; y < height; y++) {
	for (int x = 0; x < width; x++) {
	  imRef(big,x,y) = imRef(input, x,y);
	  imRef(big,x,y+height) = imRef(input2, x,y);
	}
      }
   }
   else  if(ratio_horizontal==2){
      for (int y = 0; y < height; y++) {
	for (int x = 0; x < width; x++) {
	  imRef(big,x,y) = imRef(input, x,y);
	  imRef(big,x+width,y) = imRef(input2, x,y);
	}
      }
   }
  
  
   savePPM(big, outname);
   sprintf(appel, "convert %s %s.png  ;", outname, outname );
   system(appel);
  }

  printf("done.\n");
  
  return 0;
}

 
