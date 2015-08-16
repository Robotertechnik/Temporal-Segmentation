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

#ifndef SEGMENT_VIDEO
#define SEGMENT_VIDEO

#include <cstdlib>
#include <image.h>
#include <misc.h>
#include <math.h>
#include <filter.h>
#include "label-components.h"
#include <assert.h>
#include <time.h>
#include <sys/time.h>

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )  
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define max3( a, b, c) ( ((max(a,b)) > (c)) ? (max(a,b)) : (c) )  

/*
 * Segment an image based on the previous one
 * Returns a color image representing the segmentation.
 */
void segment_video(image<rgb> *im, // in :  image to segment
		   int min_size, // in : minimum component size (enforced by post-processing stage)			  
		   edge *edges, // in : edges of the graph G defined on the image to segment + edges for markers 
		   int num_edges, // in : total nb of edges of the graph G
		   int nedges_bef_corr, // in : total nb of edges of G without counting the edges associated to "correcting markers" 
		   int ** labels1, // out: resulting labeling
		   int *  labels, // in : seeds/markers values
		   rgb *  colors, // in : colors of the seeds/markers 
		   int *  nlabels, // in : nb of regions in S_t : the segmentation of the previous frame
		                   //   out: nb_of regions in S_{t+1} resulting segmentation 
		   int    nlabels2, // in : nb of regions in S'_{t+1} : the independant segmentation of the current frame 
		   int ** image_size, // in-out: array containing each region size  
		   image<rgb> **output, // out : resulting color segmentation S_{t+1}
		   int *** centroids, // in : centroids of S_t,  out : centroids of S_{t+1} 
		   int ** optx, // out : optical flow displacements along the x axis  
		   int ** opty, // out : optical flow displacements along the y axis  
		   uchar * semantic, // in : previous semantic result
		   image<uchar> **semantic_result, // out : current semantic result 
		   int MAX_REGION_DISPLACEMENT, // in : threshold on optical flow values for propagation of semantic results
		   int nb_class // in : nb_classes in case of semantic segmentation
		   ) {
  int width = im->width();
  int height = im->height();
  int N = width*height;
  int num_vertices = N+(*nlabels)+nlabels2;
  int label, j;
  struct timeval tim; 
  int l1label; 
 
  universe * u = new universe(num_vertices);
  int num_flat_edges = (width-1)*height + width*(height-1) + 2*(width-1)*(height-1);

  // (1) Segmentation 
  


  // (a) for each correction edge, in non-decreasing weight order...
  for (int i = nedges_bef_corr; i < num_edges; i++) {
    // components connected by this edge
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if (a != b) 
      if (  !((labels[a]>=0) && (labels[b]>=0) && (labels[a]!=labels[b]))) {
	u->join(a, b);
	label = max(labels[a], labels[b]);
	labels[a] =label; 
	labels[b] =label; 
      }
  }

  // (b) for each seed edge, in non-decreasing weight order...
  for (int i = num_flat_edges; i < nedges_bef_corr; i++) {
    // components connected by this edge
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if (a != b) 
      if (  !((labels[a]>=0) && (labels[b]>=0) && (labels[a]!=labels[b]))) {
	u->join(a, b);
	label = max(labels[a], labels[b]);
	labels[a] =label; 
	labels[b] =label; 
      }
  }

   // (c) for each flat edge, in non-decreasing weight order...
  for (int i = 0; i < num_flat_edges; i++) {
    // components connected by this edge
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if (a != b) 
      if (  !((labels[a]>=0) && (labels[b]>=0) && (labels[a]!=labels[b]))) {
	u->join(a, b);
	label = max(labels[a], labels[b]);
	labels[a] =label; 
	labels[b] =label; 
      }
  }
 
  int * comp = new int[N];

  // (2) identify the color labels
  for (int i = N; i < num_vertices; i++) {
    colors[u->find(i)]=colors[i];
    semantic[u->find(i)]=semantic[i];
  }
  
  // Useful for computing the final segmentation of current frame
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      j = y * width + x;
      comp[j] = u->find(j);
    }
  }

  // (3) First step of optical flow and semantic propagation
  int l1nb_cpts = *nlabels;
   
  if (nb_class > 0) { 
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
	j = y * width + x;    
	l1label = labels[comp[j]]-N;
	if (l1label>=l1nb_cpts) { // new label
	  (*optx)[j] = 0;
	  (*opty)[j] = 0; 
	}
      else { //label propagation
	(*optx)[j] = - (*centroids)[l1label][0];
	(*opty)[j] = - (*centroids)[l1label][1];
      } 
      }  
    }
  }
  
  // (4) Clean-up segmentation
  remove_small_components(comp, width, height, min_size, colors, output, image_size, nlabels, labels1, semantic, semantic_result);
  
  
  // (5) Compute mean of regions and centroids
  int *im1_means_r = new int[N]();
  int *im1_means_g = new int[N]();
  int *im1_means_b = new int[N]();
  int *centroids_tmp = new int[2*N]();
  
  int jj;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      j = y * width + x;
      im1_means_r[(*labels1)[j]]+=imRef(im,x,y).r; 
      im1_means_g[(*labels1)[j]]+=imRef(im,x,y).g; 
      im1_means_b[(*labels1)[j]]+=imRef(im,x,y).b; 
      centroids_tmp[(*labels1)[j]]+=x; 
      centroids_tmp[N+(*labels1)[j]]+=y; 
      (*centroids)[(*labels1)[j]][2]=j;  
    }
  }
  
  for (jj = 0; jj < (*nlabels); jj++) {
    (*centroids)[jj][0]= (centroids_tmp[jj]/(*image_size)[jj]);
    (*centroids)[jj][1]= (centroids_tmp[jj+N]/(*image_size)[jj]);
    (*centroids)[jj][3]= max3(im1_means_r[jj], im1_means_g[jj], im1_means_b[jj])/(*image_size)[jj];
  }
   

  // (6) if semantic segmentation : 
  if (nb_class > 0) {

    // (a) 2nd step for Optical flow and semantic propagation
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
	j = y * width + x;
	l1label = labels[u->find(j)]-N;
	if (l1label>=l1nb_cpts) {
	  (*optx)[j] = 0;
	  (*opty)[j] = 0; 
	}
	else {
	  (*optx)[j] += (*centroids)[(*labels1)[j]][0];
	  (*opty)[j] += (*centroids)[(*labels1)[j]][1];
	  // optical flow too large: do not propagate semantic label 
	  if ((abs((*optx)[j])+abs((*opty)[j]) > MAX_REGION_DISPLACEMENT))
	    { 
	      (*optx)[j] = 0;
	      (*opty)[j] = 0;  
	    }
	} 
      }  
    }
    
    // (b) loop to make sure that semantic_results has one label per region 
    int *indic_region_off = new int [(*nlabels)]();
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
	j= y*width+x;
	if (( (*optx)[j]==0) && ((*opty)[j] == 0))
	  indic_region_off[(*labels1)[j]]++;
      }
    }
  
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
	j= y*width+x;
	if (((*optx)[j]==0) && ((*opty)[j] == 0)&& (indic_region_off[(*labels1)[j]]>min_size))
	  imRef((*semantic_result),x,y)= 0;
      }
    }
    delete [] indic_region_off;
  }
  
  delete [] centroids_tmp;
  delete [] im1_means_r;
  delete [] im1_means_g;
  delete [] im1_means_b;
  delete u;
}

#endif
