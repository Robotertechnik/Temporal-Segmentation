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

#ifndef GRAPH_MATCHING
#define GRAPH_MATCHING

#include <cstdlib>
#include <image.h>
#include <misc.h>
#include <imutil.h>
#include <math.h>
#include <filter.h>
#include <iostream>
#include "segment-graph.h"
#include "list_utils.h"
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include "basic-morpho.h"

int SMALL_SIZE; 
int SMALL_SIZE_NEW; // threshold for new region appearance 
int EROSION_PARAM;  
//#define VERBOSE



graph build_graph(label_image l1, // input : segmentation of image at time t
		  label_image l2, // input : independant segmentation of image at time t+1
		  int ** centroids1, // input : contains centroids coordinates of segmentation l1 
		  int ** centroids2, // input : contains centroids coordinates of segmentation l2 
		  int DIST_CENTROID_MIN_DIST,
		  int DIST_CENTROID_MAX_DIST,
		  int LARGE_SIZE) 
		  
{  
  // (1) Variables declarations and initialization
  int x, y, w, k, y1, y2, ny;
  int N = l1.height*l1.width;
  int ** Union = new int * [l1.nb_cpts];
  for (x = 0; x <l1.nb_cpts ; x++)
    Union[x] = new int [l2.nb_cpts];
  int ** Inter = new int * [l1.nb_cpts];
  for (x = 0; x <l1.nb_cpts ; x++)
    Inter[x] = new int [l2.nb_cpts]();

  int * Marked = new int[N]();
  graph G;  
  Lifo *LIFO= CreeLifoVide(N);
  
  // (2) add edges to the graph: region that are not overlaping but with close centroids
  
  //a) compute a distance matrix for centroids
  int ** DistanceCentroids = new int * [l1.nb_cpts];
  for (x = 0; x <l1.nb_cpts ; x++)
    DistanceCentroids[x] = new int [l2.nb_cpts];
  for (x = 0; x <l1.nb_cpts ; x++) {
    for (y = 0; y <l2.nb_cpts ; y++) {
      DistanceCentroids[x][y]= sqrt((centroids1[x][1]-centroids2[y][1])*(centroids1[x][1]-centroids2[y][1])+ (centroids1[x][0]-centroids2[y][0])*(centroids1[x][0]-centroids2[y][0]));
      Union[x][y]=2*(l2.sizes[y]+ l1.sizes[x]); 
    }
  }

  //b) if the distance is less than a threshold, compute an overlap between them
  int x_, y_;
  int cpt_Marked = 0;
  for (x_ = 0; x_ <l1.nb_cpts ; x_++) {
    for (y_ = 0; y_ <l2.nb_cpts ; y_++) {
      if ((DistanceCentroids[x_][y_]<DIST_CENTROID_MIN_DIST)||( (DistanceCentroids[x_][y_]<DIST_CENTROID_MAX_DIST) && (min(l1.sizes[x_],l2.sizes[y_])>LARGE_SIZE))) {
	// i- compute similarity (y_ , x_)
	cpt_Marked++;
	x = centroids2[y_][2];
	Marked[x]=cpt_Marked;
	LifoPush(LIFO,x);
	while (!LifoVide(LIFO)) {
	  w = LifoPop(LIFO);
	  for (k = 0; k < 8; k += 1) {
	    y = voisin(w, k, l1.width, N);
	    if ((y != -1) && (Marked[y] < cpt_Marked)) 
	      {
		// for a given segment of label l2[w] in image 2
		if  (l2.I[y] == l2.I[w])
		{
		  Marked[y]=cpt_Marked;
		  
		  if (  min(l1.sizes[x_],l2.sizes[y_]) <LARGE_SIZE){
		    // a) compute normalized coordinates
		    y1 = y/l1.width + centroids1[x_][1]-centroids2[l2.I[y]][1];
		    y2 = y%l1.width + centroids1[x_][0]-centroids2[l2.I[y]][0];
		    
		    // b) centering the regions before comparison
		    if ((y1>=0) && (y2>=0) && (y1<l1.height) && (y2<l1.width)) { 
		      ny = y1*l1.width+y2;
		      Inter[l1.I[ny]][l2.I[w]]++; 
		    Union[l1.I[ny]][l2.I[w]]--; 
		    }
		  }
		  LifoPush(LIFO,y);		  
		}
	      } // if y ...
	  } // for k ... 
	} // while (! LifoVide(LIFO)) 
	
	// ii- compute similarity (x_ , y_)
	cpt_Marked++;
	x = centroids1[x_][2];

	Marked[x]=cpt_Marked;
	LifoPush(LIFO,x);
	while (!LifoVide(LIFO)) {
	  w = LifoPop(LIFO);
	  for (k = 0; k < 8; k += 1) {
	    y = voisin(w, k, l1.width, N);
	    if ((y != -1) && (Marked[y] < cpt_Marked)) 
	      {
		// for a given segment of label l1[w] in image 1
		if  (l1.I[y] == l1.I[w])// to do : pas de decalage pour des regions de grosse taille 
		  {
		    Marked[y]=cpt_Marked;
		     if (  min(l1.sizes[x_],l2.sizes[y_]) <LARGE_SIZE){
		    // a) compute normalized coordinates
		    y1 = y/l1.width + centroids2[y_][1]-centroids1[l1.I[y]][1];
		    y2 = y%l1.width + centroids2[y_][0]-centroids1[l1.I[y]][0];
		    
		    // b) centering the regions before comparison
		    if ((y1>=0) && (y2>=0) && (y1<l1.height) && (y2<l1.width)) { 
		      ny = y1*l1.width+y2;
		      Inter[l1.I[w]][l2.I[ny]]++; 
		      Union[l1.I[w]][l2.I[ny]]--;
		    }}
		    else {//no shift for large size regions
		    Inter[l1.I[w]][l2.I[y]]++; 
		    Union[l1.I[w]][l2.I[y]]--;
		  }

		    LifoPush(LIFO,y);
		  }
	      } // if y ...
	  } // for k ... 
	} // while (! LifoVide(LIFO)) 
      }
    }
  }
  
  // (3) count the number of edges
  G.nbE = 0; 
  for (x = 0; x <l1.nb_cpts ; x++){
    for (y = 0; y <l2.nb_cpts ; y++){
      if (Inter[x][y]!=0) 
	G.nbE ++;
    }
  }
  
  // (4) build the edges
  G.nbV = l1.nb_cpts+l2.nb_cpts;
  G.edges = new edge[G.nbE];
 
  k=0;
  float a = 0;
  #ifdef VERBOSE
    printf("\n [%2d]  ",0); 
    for (x = 0; x <l2.nb_cpts ; x++)
      printf("[%4d] ", x+l1.nb_cpts);
    printf("\n [%2d]  ",0); 
  #endif
   for (x = 0; x <l1.nb_cpts ; x++){
     for (y = 0; y <l2.nb_cpts ; y++){
       if (Inter[x][y]!=0) {
 	 // if the segments are similar, the score is low
	 G.edges[k].w = (Union[x][y] *DistanceCentroids[x][y]/4.0 )/(Inter[x][y]+0.1); 
	 if (DistanceCentroids[x][y]>=DIST_CENTROID_MIN_DIST) 
	   G.edges[k].w = (Union[x][y] *DistanceCentroids[x][y]/10.0 )/(Inter[x][y]+0.1) + 5*fabs(centroids2[y][3]-centroids1[x][3])/255.0;
	 G.edges[k].a = x; 
	 G.edges[k].b = y+l1.nb_cpts;
	
         #ifdef VERBOSE
	   printf("%6.1f ", G.edges[k].w);
         #endif
	 k++;
       }
       #ifdef VERBOSE
       else printf("%6.1f ", a);
       #endif
    } //end for y
    #ifdef VERBOSE
      printf("\n [%2d]  ",x+1);
    #endif 
   }// end for x


  // free memory
  LifoTermine(LIFO);

  for (x = 0; x <l1.nb_cpts ; x++)
    delete [] Union[x];
  delete [] Union;
  for (x = 0; x <l1.nb_cpts ; x++)
    delete [] Inter[x];
  delete [] Inter;
  delete [] Marked;
  for (x = 0; x <l1.nb_cpts ; x++)
    delete [] DistanceCentroids[x];
  delete [] DistanceCentroids;
  
  return G;
}


/**********************************************************************/
universe * labeling_graph(graph G, int *labels)
{
  int label;
  // make a disjoint-set forest
  universe * u = new universe(G.nbV); 
  
  // for each edge, in non-decreasing weight order...
  for (int i = 0; i < G.nbE; i++) {
    // components connected by this edge
    int a = u->find(G.edges[i].a);
    int b = u->find(G.edges[i].b);
    if (a != b) 
      if (  !((labels[a]>=0) && (labels[b]>=0) && (labels[a]!=labels[b]))) {
	
	label = max(labels[a], labels[b]);
	if (label == labels[b]) u->join_root(b,a);
	else u->join_root(a,b);
	
	labels[a] =label; 
       	labels[b] =label; 
      }
  }
  return u;
}
  
/**********************************************************************/
// case 0 function creating a new marker for a region of image2 that is not present in image 1 

void label_comp_case0(int x, 
		      label_image l1,
		      label_image l2,
		      edge ** edges,
		      int * num_edges,
		      int * labels, 
		      rgb *colors, 
		      uchar * semantic, 
		      Lifo * LIFO, 
		      unsigned int ** Mrk, 
		      int * nb_calls 
		      ) {  
  int N = l1.width*l1.height;  
  int y,w,k, y1, y2;
  (*nb_calls) ++;
  (*Mrk)[x] = *nb_calls;
  LifoPush(LIFO, x);
  // small component appears with label of image 2
  
#ifdef VERBOSE
  printf("CASE 0 : NEW LABEL : labeling of region l2=%d by the same label \n", l2.I[x]+l1.nb_cpts);
#endif
  while (!LifoVide(LIFO)) {
    w = LifoPop(LIFO); 
    for (k = 0; k < 8; k += 1) { 
      y = voisin(w, k, l1.width, N);
      if ((y != -1) && ((*Mrk)[y] < *nb_calls) && (l2.I[y] == l2.I[w]))
	{
	  (*Mrk)[y] = *nb_calls;
	  LifoPush(LIFO, y);
	  y1 = y%l1.width; y2 = y/l1.width;
	  (*edges)[(*num_edges)].a = y;
	  (*edges)[(*num_edges)].b = l2.I[y]+N+l1.nb_cpts;
	  (*edges)[(*num_edges)].w = -1;
	  labels[l2.I[y]+l1.nb_cpts+N]= l2.I[y]+N+l1.nb_cpts;
	  colors[l2.I[y]+N+l1.nb_cpts]= imRef(l2.colorI, y1,y2 );
	  semantic[l2.I[y]+N+l1.nb_cpts]= 0;
	  labels[y]= l2.I[y]+l1.nb_cpts;
	  (*num_edges)++;
	} 
    } 
  } 
}


/**********************************************************************/
// case 1 function creating markers for the propagation of colors of image1 in image2 
// x: index of one pixel of the region of interest of image2 to label

void label_comp_case1(int x, 
		      int v2, 
		      label_image l1, 
		      label_image l2, 
		      edge ** edges,
		      int * num_edges, 
		      int * labels,
		      rgb *colors,
		      uchar * semantic,
		      Lifo * LIFO,
		      unsigned int ** Mrk,
		      int *nb_calls) {  
  
  int N = l1.width*l1.height;  
  int y,w,k,x1,x2;
  (*nb_calls) ++;
  
  #ifdef VERBOSE
  printf("CASE 1 : Match REGIONS : labeling of region l2=%d by l1 = %d \n", l2.I[x]+l1.nb_cpts, l1.I[v2]);
#endif
  (*Mrk)[x] = *nb_calls;
  LifoPush(LIFO, x);

  while (!LifoVide(LIFO)) {
    w = LifoPop(LIFO);
    for (k = 0; k < 8; k += 1) { /* parcourt les voisins en 8-connexite */
      y = voisin(w, k, l1.width, N);
      if ((y != -1) && ((*Mrk)[y] < *nb_calls) && (l2.I[y] == l2.I[w]))
	{
	  (*Mrk)[y] = *nb_calls;
	  LifoPush(LIFO,y);
	  (*edges)[(*num_edges)].a = y;
	  (*edges)[(*num_edges)].b = l1.I[v2]+N;
	  (*edges)[(*num_edges)].w = -1;
 	  labels[l1.I[v2]+N]= l1.I[v2]+N;
	  x1 = v2%l1.width; x2 = v2/l1.width;
	  colors[l1.I[v2]+N]= imRef(l1.colorI, x1, x2 );
	  semantic[l1.I[v2]+N]= imRef(l1.semanticI, x1, x2 );
	  labels[y]= l1.I[v2];
	  (*num_edges)++;
	} /* if y ... */
    } /* for k ... */
  } /* while (! LifoVide(LIFO)) */
}

/**********************************************************************/
// case 2 function creating several markers from image1 for labeling one region of image2 
void label_comp_case2(int x, int * list_neigh_regions, int nb_regions, label_image l1, label_image l2, edge ** edges, int * num_edges, int * labels, rgb *colors, uchar * semantic, int ** centroids1, int ** centroids2, Lifo * LIFO, unsigned int ** Mrk,
		      int *nb_calls) {  
  
  int N = l1.width*l1.height;  
  int i,j,y,w,k, y1, y2, ny;
  (*nb_calls) ++;

  #ifdef VERBOSE
  printf("CASE 2 : MARKERS : labeling of region l2=%d by ", l2.I[x]+l1.nb_cpts);
  for (i=0;i<nb_regions;i++)
    printf(" %d, ", list_neigh_regions[i]);
  printf(" \n");
  #endif
  (*Mrk)[x] = *nb_calls;
  LifoPush(LIFO, x);
  
  // compute a normalized centroid
  int cx=0;int cy=0; double size=0;
  for (i=0;i<nb_regions;i++)
    size +=l1.sizes[list_neigh_regions[i]];

  if (nb_regions>1) {
    for (i=0;i<nb_regions;i++) {  
      j=list_neigh_regions[i];
      cx += centroids1[j][0] * (1.0*l1.sizes[j]/size) ;
      cy += centroids1[j][1] * (1.0*l1.sizes[j]/size) ;
    }
  } else { // one region: no shift of centroids 
    cx = centroids2[list_neigh_regions[0]][0];
    cy = centroids2[list_neigh_regions[0]][1];
  }

  while (!LifoVide(LIFO)) {
    w = LifoPop(LIFO);
    for (k = 0; k < 8; k += 1) {
      y = voisin(w, k, l1.width, N);
      // for every pixel of the region of x in image2
      if ((y != -1) && ((*Mrk)[y] < *nb_calls) && (l2.I[y] == l2.I[w])) { 
	
	  (*Mrk)[y] = *nb_calls;
	  LifoPush(LIFO, y);
 
	  y1 = y/l1.width + cx -centroids2[l2.I[y]][0];
	  y2 = y%l1.width + cy -centroids2[l2.I[y]][1];
	 
	  if ((y1>=0) && (y2>=0) && (y1<l1.height) && (y2<l1.width)) {
	   
	    y = y1*l1.width+y2;
	   
	    // if the marker belong to the list
	    for (i=0;i<nb_regions;i++) 
	      {
		if (l1.I[y]==list_neigh_regions[i]) { 
		  (*edges)[(*num_edges)].a = y;
		  (*edges)[(*num_edges)].b = l1.I[y]+N; 
		  (*edges)[(*num_edges)].w = -1;
		  labels[l1.I[y]+N]= l1.I[y]+N;
		  colors[l1.I[y]+N]= imRef(l1.colorI, y%l1.width,y1 );
		  semantic[l1.I[y]+N]= imRef(l1.semanticI, y2,y1 );
		  labels[y]= l1.I[y];
		  (*num_edges)++;
		}
	      }
	  }	  
      } /* if y ... */
    } /* for k ... */
  } /* while (! LifoVide(LIFO)) */
}

/**********************************************************************/
// case 3 function creating markers for the propagation of colors of image1 in image2 
// x: index of one pixel of the region of interest of image2 to label

void label_comp_case3(int x, int point1, label_image l1, label_image l2, edge ** edges, int * num_edges, int * labels, rgb *colors, uchar * semantic, Lifo * LIFO, unsigned int ** Mrk, int * nb_calls ) {  
  
  int N = l1.width*l1.height;  
  int y,w,k,y1,y2;
  (*nb_calls) ++;
  (*Mrk)[x] = *nb_calls;
  LifoPush(LIFO, x);
 
  // small component appears with label of image 2
  
  if (l2.sizes[l2.I[x]] < SMALL_SIZE_NEW) {
  #ifdef VERBOSE
    printf("CASE 3 : NEW LABEL : labeling of region l2=%d by the same label \n", l2.I[x]+l1.nb_cpts);
  #endif
      while (!LifoVide(LIFO)) {
       w = LifoPop(LIFO); 
       for (k = 0; k < 8; k += 1) { 
	 y = voisin(w, k, l1.width, N);
	 if ((y != -1) && ((*Mrk)[y] < *nb_calls) && (l2.I[y] == l2.I[w]))
	   {
	     (*Mrk)[y] = *nb_calls;
	     LifoPush(LIFO, y);
	     y1 =  y%l1.width; y2 = y/l1.width;
	     (*edges)[(*num_edges)].a = y;
	     (*edges)[(*num_edges)].b = l2.I[y]+N+l1.nb_cpts;
	     (*edges)[(*num_edges)].w = 0;
	     labels[l2.I[y]+l1.nb_cpts+N]= l2.I[y]+N+l1.nb_cpts;
	     colors[l2.I[y]+N+l1.nb_cpts]= imRef(l2.colorI, y1, y2 );
	     semantic[l2.I[y]+N+l1.nb_cpts]= 0;
	     labels[y]= l2.I[y]+l1.nb_cpts;
	     (*num_edges)++;
	   } 
       } 
     } 
   }
   // too large component or largest component: label of image 1
   else {
        
#ifdef VERBOSE
     printf("CASE 3 : Propagation : labeling of a part of region l2=%d by the label %d\n", l2.I[x]+l1.nb_cpts, l1.I[point1]);
#endif
     while (!LifoVide(LIFO)) {
       w = LifoPop(LIFO);
       for (k = 0; k < 8; k += 1) { 
	 y = voisin(w, k, l1.width, N);
	 if ((y != -1) && ((*Mrk)[y] < *nb_calls) && (l2.I[y] == l2.I[w]))
	   {
	     (*Mrk)[y] = *nb_calls;
	     LifoPush(LIFO, y);
	     if( (l2.sizes[l2.I[x]] < SMALL_SIZE)||(l1.I[y]==l1.I[point1])) { // just added 
	     (*edges)[(*num_edges)].a = y;
	     (*edges)[(*num_edges)].b = l1.I[point1]+N;
	     (*edges)[(*num_edges)].w = -1;
	     labels[l1.I[point1]+N]= l1.I[point1]+N;
	     colors[l1.I[point1]+N]= imRef(l1.colorI, point1%l1.width,point1/l1.width );
	     semantic[l1.I[point1]+N]= imRef(l1.semanticI, point1%l1.width,point1/l1.width );
	     labels[y]= l1.I[point1];
	     (*num_edges)++;
	     }
	   } 
       } 
     }  
   }
}




/**********************************************************************/
// case 2 function creating several markers from image1 for labeling one region of image2 
void label_comp_case2bis(int x, int z, label_image l1, label_image l2, edge ** edges, int * num_edges, int * labels, rgb *colors, uchar * semantic, Lifo * LIFO, unsigned int ** Mrk, int * nb_calls) {  
  
  // x point of disapearing region in image 1 
  // z point of corresponding region in image 2

  int N = l1.width*l1.height;  
  int i,j,y,w,k, y1, y2, ny;
  (*nb_calls) ++;
  
  #ifdef VERBOSE
     printf("CASE 4 : LARGE DISAPPEARING LABEL : labeling region l2=%d by l1=%d \n", l2.I[z]+l1.nb_cpts, l1.I[x]);
  #endif

  (*Mrk)[x] = *nb_calls;
  LifoPush(LIFO, x);

  while (!LifoVide(LIFO)) {
    w = LifoPop(LIFO);
    for (k = 0; k < 8; k += 1) {
      y = voisin(w, k, l1.width, N);
      // for every pixel of the region of x in image1
      if ((y != -1) && ((*Mrk)[y] < *nb_calls) && (l1.I[y] == l1.I[w])) { 
	
	(*Mrk)[y] = *nb_calls;
	LifoPush(LIFO, y);
	if ((l2.I[y] == l2.I[z])) 
	  {
	    y1 =  y%l1.width; y2 = y/l1.width;
	    (*edges)[(*num_edges)].a = y;
	    (*edges)[(*num_edges)].b = l1.I[y]+N; 
	    (*edges)[(*num_edges)].w = -1;
	    labels[y]= l1.I[y];
	    labels[l1.I[y]+N]= l1.I[y]+N;
	    colors[l1.I[y]+N]= imRef(l1.colorI, y1,y2 );
	    semantic[l1.I[y]+N]= imRef(l1.semanticI, y1,y2 );
	    (*num_edges)++;
	  }	  
      } /* if y ... */
    } /* for k ... */
  } /* while (! LifoVide(LIFO)) */
}


/**********************************************************************/
// case 5 function creating 
void label_comp_case5(int x,
		      int z,
		      label_image l1,
		      edge ** edges,
		      int * num_edges,
		      int * labels,
		      rgb *colors,
		      uchar * semantic, 
		      Lifo * LIFO, 
		      unsigned int ** Mrk, 
		      int *nb_calls,
		      image<uchar> * eroded_l1,
		      int ** centroids1 ) {  
  
  // x point of disappearing region in image 1 
  // z label to make disappear

  int N = l1.width*l1.height;  
  int i,j,y,w,k, y1, y2;
  (*nb_calls) ++;
  
  #ifdef VERBOSE
  fprintf(stderr,"CASE 5 :  RE-LABELing : labeling region l1=%d and erasing label %d from this region\n", l1.I[x], z);
  #endif

  (*Mrk)[x] = *nb_calls;
  LifoPush(LIFO, x);

  while (!LifoVide(LIFO)) {
    w = LifoPop(LIFO);
    for (k = 0; k < 8; k += 1) {
      y = voisin(w, k, l1.width, N);
      // for every pixel of the region of x in image1
      if ((y != -1) && ((*Mrk)[y] <*nb_calls) && (l1.I[y] == l1.I[w])) { 
	
	(*Mrk)[y] = *nb_calls;
	LifoPush(LIFO, y);
	y1 =  y%l1.width; y2 = y/l1.width;
	if  ((labels[y] == z) && (imRef(eroded_l1, y1,y2 )==0)) {
	  (*edges)[(*num_edges)].a = y;
	  (*edges)[(*num_edges)].b = l1.I[y]+N; 
	  (*edges)[(*num_edges)].w = 0;
	  labels[y]= l1.I[y];
	  labels[l1.I[y]+N]= l1.I[y]+N;
	  colors[l1.I[y]+N]= imRef(l1.colorI, y1,y2);
	  semantic[l1.I[y]+N]= imRef(l1.semanticI,y1,y2);
	  (*num_edges)++;
	} 
      } /* if y ... */
    } /* for k ... */
  } /* while (! LifoVide(LIFO)) */
 
  
}

/**********************************************************************/
int argmax(int * T, // array
	   int n,  //nb of elements
	   int forbiden_index
	   ) {
  int argmax = 0;
  int maxi = T[0];
  for (int x = 1; x < n; x++) {
    if((T[x]>maxi)&&(x-1!= forbiden_index)) {
      maxi = T[x];
      argmax = x; 
    }
  }
  return argmax-1;
}

/**********************************************************************/
  // returns a seed image
void graph_matching(graph G, // in : Graph for the matching procedure
		    label_image l1, // in : segmentation S_t
		    label_image l2, // in : segmentation S'_{t+1}
		    int ** centroids1, 
		    int ** centroids2, 
		    edge ** edges, // in: parwise edges of the graph G_{t+1}. out: pairwise (unchanged) + unary edges 
		    int * num_edges, // out: total number of edges G_{t+1}
		    int * nedges_bef_corr, // out: total nb of edges G_{t+1} without counting the edges associated to "correcting markers" 
		    int *seeds, // out: markers values for the final segmentation step
		    rgb *colors, // out: markers colors for the final segmentation step
		    uchar * semantic, // out: semantic values for the final segmentation step
		    int static_mode) // in: 0 if the camera is moving, 1 otherwise
{

  // (0) some initializations
  int i,j,k,z,x,y, nb_regions;
  int nb_calls = 0;
  int N = l1.width*l1.height; 
  Lifo *LIFO= CreeLifoVide(N); 
  unsigned int * Mrk = new unsigned int[N]();
 
  for (i =0 ; i < G.nbV; i++) 
    seeds[i+N]= 0;
  
  SMALL_SIZE_NEW= (250.0*N)/(320*240); // threshold for new region appearance
  SMALL_SIZE= (500.0*N)/(320*240); 
  EROSION_PARAM= (14.0*N)/(320*240);
  int LARGE_SIZE= (1250.0*N)/(320*240); 
   
  if  (static_mode > 0) { 
    SMALL_SIZE_NEW = SMALL_SIZE_NEW/2.0;
    SMALL_SIZE = SMALL_SIZE/2.0;
    EROSION_PARAM = EROSION_PARAM/2.0;
    LARGE_SIZE = LARGE_SIZE/2.0;
 }

  // (1) compute an eroded map of regions of the image 1

  image<uchar>* eroded_l1 = new image<uchar>(l1.width, l1.height);
  boundary(l1.I, l1.width, l1.height, eroded_l1);
  ldilatdisc(eroded_l1, EROSION_PARAM*EROSION_PARAM, false);
  
  image<uchar>* eroded_l2 = new image<uchar>(l1.width, l1.height);
  boundary(l2.I, l1.width, l1.height, eroded_l2);
  ldilatdisc(eroded_l2, EROSION_PARAM*EROSION_PARAM, false);

  for (int y = 0; y < l1.height; y++) 
    for (int x = 0; x < l1.width; x++) {
      if ((imRef(eroded_l1,x,y) == 255)||(imRef(eroded_l2,x,y) == 255))
	imRef(eroded_l1,x,y) =255;
      else imRef(eroded_l1,x,y) =0;
    }
  
  // (2) compute intermediate result: the 2 labeled graphs
  int * labels = new int [G.nbV];  
  
  // sort edges by weight
  std::sort(G.edges, G.edges + G.nbE);
  
  // computing label of image1 
  for (x=0;x<l1.nb_cpts;x++)
    labels[x]=x;
  for (x=l1.nb_cpts;x<G.nbV;x++)
    labels[x]=-1;
  universe * u1 = labeling_graph(G, labels);
  
  // computing label of image2  
  for (x=0;x<l1.nb_cpts;x++)
    labels[x]=-1;
  for (x=l1.nb_cpts;x<G.nbV;x++)
    labels[x]=x;
  universe * u2 = labeling_graph(G, labels);
  
  #ifdef VERBOSE
    printf("--------components u2----------- \n");
    for (i = 0; i <l1.nb_cpts ; i++) 
      printf("(%d:) %02d [s %d]\n",i, u2->find(i), u2->size(u2->find(i)));
    printf("\n --------components u1----------- \n");
    for (i =l1.nb_cpts ; i<G.nbV; i++) 
      printf("(%d:) %02d [s %d]\n",i, u1->find(i), u1->size(u1->find(i))); 
  #endif
  
  // (3) pre-processing to create the unary weights, and seeds  
  
  // a) create a list of corresponding regions for every region of image1
  
  int * list_markers1 = new int[l1.nb_cpts]();
  for (i =0 ; i < l1.nb_cpts; i++) 
    list_markers1[i]= u2->find(i);
  
  // b) create a list of corresponding regions for every region of image2
  int *sl2 = new int [l2.nb_cpts]();
  int ** list_markers2 = new int* [l2.nb_cpts]; 
  for (x = 0; x <l2.nb_cpts ; x++)
    list_markers2[x] = new int [l1.nb_cpts]();
  
  for (i =0 ; i<l1.nb_cpts; i++) {
    if (u2->find(i) >= l1.nb_cpts){ // the region has at least one corresponding region
      list_markers2[u2->find(i)-l1.nb_cpts][sl2[u2->find(i)-l1.nb_cpts]]= i;
      // cout <<" list[  " << u2->find(i) << " ][ "<< sl2[u2->find(i)-l1.nb_cpts] <<" ]= "<< i << endl;
      sl2[u2->find(i)-l1.nb_cpts]++;
    }
  }
  
  *num_edges = (l1.width-1)*l1.height + l1.width*(l1.height-1) + 2*(l1.width-1)*(l1.height-1); 
 
  //--------------------------------------------------------------------------------------------
  // (4) Markers creation 

  // (a) Prevent large regions disappearing 

  // In cases where  s' has several corresponding regions s_1, ..., s_r (Case 2)
  // and one of the corresponding regions s is large : 
  // Instead of placing markers for each corresponding regions, just label the region with the label of s.  
  for (i = 0 ; i <l1.nb_cpts ; i++) {
    z = u2->find(i);
    j = u1->find(z);
    if (u1->size(j) != 2) {
      if (j!=i ) {
	//printf("The label %d (region of size %d )is disappearing \n", i,  l1.sizes[i]);
	if ( l1.sizes[i] > LARGE_SIZE ) {
	  label_comp_case2bis(centroids1[i][2], centroids2[z-l1.nb_cpts][2], l1, l2, edges,  num_edges, seeds, colors, semantic, LIFO, &Mrk, &nb_calls);  
	}
      }
    }
  }

  // (b) general loop
  for (i = l1.nb_cpts ; i <G.nbV ; i++) {
    j = u1->find(i);
    z = u2->find(j);
    x = centroids2[i-l1.nb_cpts][2]; 
    if (u1->size(j) == 2) { 
      // CASE 1: one to one mapping
      if (u2->size(z) == 2) 
	label_comp_case1(x, centroids1[j][2], l1, l2, edges, num_edges, seeds, colors, semantic, LIFO, &Mrk, &nb_calls);
      else {// CASE 2 : markers
	nb_regions = u2->size(l2.I[x]+l1.nb_cpts)-1;
	if (nb_regions==0) {list_markers2[l2.I[x]][0]=j ; nb_regions++; }
	label_comp_case2(x, list_markers2[l2.I[x]], nb_regions, l1, l2, edges, num_edges, seeds, colors, semantic, centroids1, centroids2, LIFO, &Mrk, &nb_calls);
      }
    }
    else {
      y = centroids1[j][2]; 
      if(j>=l1.nb_cpts) // CASE 0: no matching region: new label in image 2
	label_comp_case0(x, l1, l2, edges, num_edges, seeds, colors, semantic, LIFO, &Mrk, &nb_calls );
      else {// CASE 3: propagation of labels of image 1 or new label
	label_comp_case3(x, y, l1, l2, edges, num_edges, seeds, colors, semantic, LIFO, &Mrk, &nb_calls );
      }
    }
  }
  
  // (c) Check large differences between corresponding region sizes to correct markers 

  int label2;
  int *difference = new int [l1.nb_cpts](); 
  int **coverage = new int* [l1.nb_cpts];
  for (int y = 0; y < l1.nb_cpts; y++) 
     coverage[y] = new int [G.nbV+1](); // stores the nb of pixel of each label in region y  
  
  for (int y = 0; y < l1.height; y++) {
    for (int x = 0; x < l1.width; x++) {
      j=x+y*l1.width;
      if  (seeds[j] != l1.I[j]) difference[l1.I[j]]++;
      coverage[l1.I[j]][seeds[j]+1]++;
    }
  }
  
  *nedges_bef_corr = *num_edges;
  for (i = 0 ; i <l1.nb_cpts ; i++) {
    if ( ( ( (1.0*difference[i])/l1.sizes[i] > 0.5)&& (difference[i]> SMALL_SIZE))  || (difference[i]> LARGE_SIZE)) {
      // fprintf(stderr,"Noticeable difference[%d] = %d in a region of size %d \n",i, difference[i], l1.sizes[i]);
      y = centroids1[i][0]+centroids1[i][1]*l1.width;
      label2 = argmax(coverage[i], G.nbV+1, i);
      label_comp_case5(centroids1[i][2], label2, l1, edges,  num_edges, seeds, colors, semantic, LIFO, &Mrk, &nb_calls, eroded_l1, centroids1);   
    }
  }
 
  // (d) create an eroded map of the large non-seeded area and label then with l1.
  image<uchar>* eroded_marqueurs = new image<uchar>(l1.width, l1.height);  
  for (i =0 ; i < N; i++) 
    if (seeds[i]==-1) 
      imRef(eroded_marqueurs, i%l1.width, i/l1.width) = 255;
  
  ldilatdisc(eroded_marqueurs, EROSION_PARAM*EROSION_PARAM, true);
   
  for (y =0 ; y < N; y++) 
    if ((imRef(eroded_marqueurs, y%l1.width, y/l1.width) == 0)&& (seeds[y]==-1)){    
      (*edges)[(*num_edges)].a = y;
      (*edges)[(*num_edges)].b = l1.I[y]+N; 
      (*edges)[(*num_edges)].w = -1;
      seeds[y]= l1.I[y];
      seeds[l1.I[y]+N]= l1.I[y]+N;
      colors[l1.I[y]+N]= imRef(l1.colorI, y%l1.width,y/l1.width );
      semantic[l1.I[y]+N]= imRef(l1.semanticI, y%l1.width,y/l1.width );
      (*num_edges)++;
    }
 
  for (int y = 0; y < N; y++) 
    seeds[y]=-1;
  
   // (6) free memory
   eroded_marqueurs->image::~image();
   eroded_l1->image::~image();
   eroded_l2->image::~image();
   delete [] G.edges;
   LifoTermine(LIFO);
   delete [] labels; 
   delete u1,u2;
   delete [] list_markers1;
   for (x = 0; x <l2.nb_cpts ; x++)
     delete [] list_markers2[x];
   delete [] list_markers2;
   delete [] sl2;
   delete [] Mrk;
}






#endif


