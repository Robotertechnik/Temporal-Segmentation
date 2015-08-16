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

#ifndef SEMANTIC_SEGMENT
#define SEMANTIC_SEGMENT

#include <cstdlib>
#include <image.h>
#include <misc.h>
#include <math.h>
#include <iostream>
#include <fstream>


/*
 * Saves the semantic segmentation of an image and evaluate pixel classification if a ground truth file is available
 */
void semantic_segment(image<rgb> *input, // in : input image 
		      image <uchar>*predictions, // in : noisy semantic predictions 
		      int *labels, // in : input superpixels segmentation
		      int num_ccs, // in : nb of superpixels in labels
		      int nb_classes, // in : nb of different possible semantic classes +1
		      char * outname, // in : name of the output file for the resulting semantic segmentation 
		      image<uchar> **output, //in : previous semantic result; out : current semantic result 
		      char * gt_name, // in : ground truth name
		      int *** Confusion) { // in and out : confusion matrix 

  int width = input->width();
  int height = input->height();
  int N = width*height;
  int i,j, argmax;
  
  uchar * presence_class = new uchar [nb_classes]();
  float * colormap = new float [nb_classes];
  float * colormap_g = new float [nb_classes];
  float * colormap_b = new float [nb_classes];
  
  
  if (nb_classes==35) {
    /* SIFTFLOW COLOR MAP */
    colormap[1] = 0.0; colormap_g[1] =0.0; colormap_b[1]= 0.0;
    colormap[2] = 0.5; colormap_g[2]= 0.5; colormap_b[2]= 0.5; // awning
    colormap[3] = 0.9; colormap_g[3]= 0.3; colormap_b[3]= 0.3; // balcony
    colormap[4] = 0.8; colormap_g[4]= 0.3; colormap_b[4]= 0.2; // bird
    colormap[5] = 0.4; colormap_g[5]= 0.4; colormap_b[5]= 0.8; // boat
    colormap[6] = 0.5; colormap_g[6]= 0.9; colormap_b[6]= 0.9; // bridge
    colormap[7] = 0.7; colormap_g[7]= 0.7; colormap_b[7]= 0.3; // building 
    colormap[8] = 0.4; colormap_g[8]= 0.7; colormap_b[8]= 0.8; // bus
    colormap[9] = 0.4; colormap_g[9]= 0.4; colormap_b[9]= 0.8; // car  
    colormap[10]= 0.8; colormap_g[10]= 0.6; colormap_b[10]= 0.6; // cow
    colormap[11]= 0.9; colormap_g[11]= 0.7; colormap_b[11]= 0.9; // crosswalk
    colormap[12]= 0.9; colormap_g[12]= 0.9; colormap_b[12]= 0.5; // desert
    colormap[13]= 0.5; colormap_g[13]= 0.3; colormap_b[13]= 0.0; // door  
    colormap[14]= 0.6; colormap_g[14]= 0.5; colormap_b[14]= 0.1; // fence
    colormap[15]= 0.7; colormap_g[15]= 0.7; colormap_b[15]= 0.1; // field
    colormap[16]= 0.0; colormap_g[16]= 0.9; colormap_b[16]= 0.0; // grass
    colormap[17]= 0.0; colormap_g[17]= 0.2; colormap_b[17]= 0.2; // moon
    colormap[18]= 0.7; colormap_g[18]= 0.5; colormap_b[18]= 0.3; // mountain
    colormap[19]= 1.0; colormap_g[19]= 0.0; colormap_b[19]= 0.3; // person 
    colormap[20]= 0.3; colormap_g[20]= 0.7; colormap_b[20]= 0.1; // plant
    colormap[21]= 0.4; colormap_g[21]= 0.2; colormap_b[21]= 0.2; // pole 
    colormap[22]= 0.1; colormap_g[22]= 0.4; colormap_b[22]= 0.9; // river
    colormap[23]= 0.3; colormap_g[23]= 0.3; colormap_b[23]= 0.3; // road 
    colormap[24]= 0.5; colormap_g[24]= 0.4; colormap_b[24]= 0.2; // rock
    colormap[25]= 0.8; colormap_g[25]= 0.8; colormap_b[25]= 0.5; // sand
    colormap[26]= 0.1; colormap_g[26]= 0.1; colormap_b[26]= 0.9; // sea
    colormap[27]= 0.5; colormap_g[27]= 0.5; colormap_b[27]= 0.5; // sidewalk 
    colormap[28]= 1.0; colormap_g[28]= 0.1; colormap_b[28]= 0.1; // sign 
    colormap[29]= 0.0; colormap_g[29]= 0.7; colormap_b[29]= 0.9; // sky 
    colormap[30]= 0.9; colormap_g[30]= 0.4; colormap_b[30]= 0.3; // staircase
    colormap[31]= 0.1; colormap_g[31]= 1.0; colormap_b[31]= 0.1; // streetlight
    colormap[32]= 1.0; colormap_g[32]= 1.0; colormap_b[32]= 0.0; // sun
    colormap[33]= 0.2; colormap_g[33]= 0.8; colormap_b[33]= 0.1; // tree 
    colormap[34]= 0.1; colormap_g[34]= 0.6; colormap_b[34]= 1.0; // window 
    for (i=2;i<nb_classes;i++)
      presence_class[i]=i;
  }else {
    /* CAMVID COLOR MAP*/
    colormap[1] = 0.25; colormap_g[1]= 0.5; colormap_b[1]= 0.25; // animal
    colormap[2] = 0.75; colormap_g[2]= 0.0; colormap_b[2]= 0.5; // archway
    colormap[3] = 0; colormap_g[3]= 0.5; colormap_b[3]= 0.75;   presence_class[3] = 3;  // bicyclist 
    colormap[4] = 0; colormap_g[4]= 0.5; colormap_b[4]= 0.25; // bridge 
    colormap[5] = 0.5; colormap_g[5]= 0.0; colormap_b[5]= 0.0;presence_class[5] = 5; // building 
    colormap[6] = 0.25; colormap_g[6]= 0.0; colormap_b[6]= 0.5;presence_class[6] = 6; // car 
    colormap[7] = 0.25; colormap_g[7]= 0; colormap_b[7]= 0.75; presence_class[7] = 17; // cartluggagePram 
    colormap[8] = 0.75; colormap_g[8]= 0.5; colormap_b[8]= 0.25; presence_class[8] = 17;// child
    colormap[9] = 0.75; colormap_g[9]= 0.75; colormap_b[9]= 0.5; presence_class[9] = 9;//  column pole 
    colormap[10]= 0.25; colormap_g[10]= 0.25; colormap_b[10]= 0.5; presence_class[10] = 10; // fence 
    colormap[11]= 0.5; colormap_g[11]= 0.0; colormap_b[11]= 0.75; presence_class[11] = 18;// lane_MkgsDriv
    colormap[12]= 0.75; colormap_g[12]= 0.0; colormap_b[12]= 0.25; presence_class[12] = 18;// laneMkgsNonDriv
    colormap[13]= 0.5; colormap_g[13]= 0.5; colormap_b[13]= 0.25; presence_class[13] = 21;// Misc text 
    colormap[14]= 0.75; colormap_g[14]= 0.0; colormap_b[14]= 0.75; // Motorcycle Scooter
    colormap[15]= 0.5; colormap_g[15]= 0.25; colormap_b[15]= 0.25; // Other Mooving
    colormap[16]= 0.25; colormap_g[16]= 0.75; colormap_b[16]= 0.5; // Parking Block
    colormap[17]= 0.25; colormap_g[17]= 0.25; colormap_b[17]= 0.0; presence_class[17] = 17;// Pedestrian 
    colormap[18]= 0.5; colormap_g[18]= 0.25; colormap_b[18]= 0.5; presence_class[18] = 18; // Road 
    colormap[19]= 0.5; colormap_g[19]= 0.5; colormap_b[19]= 0.75; presence_class[19] = 18;// Road Shoulder
    colormap[20]= 0.0; colormap_g[20]= 0.0; colormap_b[20]= 0.75; presence_class[20] = 20;// Sidewalk  
    colormap[21]= 0.75; colormap_g[21]= 0.5; colormap_b[21]= 0.5; presence_class[21] = 21;// Sign Symbol  
    colormap[22]= 0.5; colormap_g[22]= 0.5; colormap_b[22]= 0.5; presence_class[22] = 22;// sky   
    colormap[23]= 0.25; colormap_g[23]= 0.5; colormap_b[23]= 0.75; presence_class[23] = 6 ;// SUV Pickup truck
    colormap[24]= 0.0; colormap_g[24]= 0; colormap_b[24]= 0.25; // Traffic Cone
    colormap[25]= 0.0; colormap_g[25]= 0.25; colormap_b[25]= 0.25; presence_class[25] = 21;// traffic light
    colormap[26]= 0.75; colormap_g[26]= 0.25; colormap_b[26]= 0.5; // train
    colormap[27]= 0.5; colormap_g[27]= 0.5; colormap_b[27]= 0.0; presence_class[27] = 27;// tree  
    colormap[28]= 0.75; colormap_g[28]= 0.5; colormap_b[28]= 0.75; presence_class[28] = 6 ; // Truck bus
    colormap[29]= 0.25; colormap_g[29]= 0; colormap_b[29]= 0.25; // tunnel
    colormap[30] = 0.75; colormap_g[30] =0.75; colormap_b[30]= 0.0; presence_class[30] = 27 ;//vegetation misc 
    colormap[31]= 0.0; colormap_g[31]= 0.0; colormap_b[31]= 0.0;// presence_class[31] = 1;// void
    colormap[32]= 0.25; colormap_g[32]= 0.75; colormap_b[32]= 0; presence_class[32] = 5;// wall
  }
    
  image<rgb> *output_color = new image<rgb>(width, height);
  int *indic_region_off = new int [num_ccs]();
  int **indic_class = new int* [nb_classes];
  
  for (i = 0; i < nb_classes; i++) 
    indic_class[i] = new int [num_ccs]();
  
  
    //(1) count the number of pixels of each class for each segment. This number is stored in indic_class
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      j=y*width+x;
      if ( imRef((*output),x,y)==0) 
	indic_region_off[labels[j]]++;
      indic_class[imRef(predictions, x, y)][labels[j]] ++;
    }
  }
    
  //(2) compute argmax for each ligne of the array indic_class, store it in the first cell of each line.
  for (int y = 0; y <num_ccs ; y++) {
    argmax = 0;
    // printf("%d \n", indic_region_off[y]);
    for (int x = 1; x < nb_classes; x++) {
      if (indic_class[x][y]> indic_class[0][y]){
        indic_class[0][y]=indic_class[x][y];
	argmax = x; 
      }
    }
    indic_class[0][y]=argmax;
  }
  
  //(3) produce the semantic color map
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      // case 1 : no temporal prediction
      //imRef((*output),x,y)=imRef(predictions,x,y)+1;
 
      // case 2 : add temporal consistency
        if (imRef((*output),x,y)==0) 
     	imRef((*output),x,y)= indic_class[0][labels[y * width + x]]+1;
    
      //else use propagated values already computed in output.
      imRef(output_color,x,y).r= colormap[imRef((*output),x,y)]*128+imRef(input,x,y).r/2;
      imRef(output_color,x,y).g= colormap_g[imRef((*output),x,y)]*128+imRef(input,x,y).g/2;
      imRef(output_color,x,y).b= colormap_b[imRef((*output),x,y)]*128+imRef(input,x,y).b/2;
    }
  }

 // (4) Evaluate results
  #ifdef GROUND_TRUTH_AVAILABLE
    int M = width*height;
    int *GT = new int [M];
    ifstream myfile;
    myfile.open(gt_name);
    //if(!myfile) 
    // cout << endl << "Failed to open file " << gt_name<<endl;
    if (myfile) {
    for (i=0;i<M ;i++){
    myfile >> GT[i];
  }
    myfile.close();
    
    for (int y = 0; y <height ; y++) {
      for (int x = 0; x < width; x++) { 
        if(presence_class[GT[y*width+x]+1] > 0)
	  (*Confusion)[presence_class[GT[y*width+x]+1]-1][imRef((*output),x,y)-1]++;
      }
    }
    
    int nb_classified_pix=0;
    int nb_correct_pix=0;
    for (int x = 0; x < nb_classes-1; x++) {
      for (int y = 0; y < nb_classes-1; y++) {
	  nb_classified_pix += (*Confusion)[x][y];
      }
	nb_correct_pix += (*Confusion)[x][x];
    }
    
    printf(" %f correct \n", (1.0*nb_correct_pix)/nb_classified_pix );
  }// if (myfile)

  delete [] GT;
  #endif

  // (5) Save results
  savePPM(output_color, outname);
  #ifdef SAVE_ALL
    char * appel = new char[1000];
    sprintf(appel, "convert %s %s.png  ;", outname, outname);
    system(appel);
    sprintf(appel, "rm %s  ;", outname);
    system(appel);
    delete [] appel;
  #endif

  // (6) free memory
 
  delete [] colormap;
  delete [] colormap_g;
  delete [] colormap_b;

  for (i = 0; i < nb_classes; i++) 
    delete [] indic_class[i];
  delete [] indic_class;
  
  output_color->image::~image(); 
}

#endif
