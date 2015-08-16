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

#ifndef LABEL_comp
#define LABEL_comp


#include <iostream>
#include "list_utils.h"
#include <image.h>
#include <misc.h>
				 

/* i : index du point dans l'image 
 * k : direction du voisin 
 * rs : taille d'une rangee 
 * nb : taille de l'image 
 * retourne -1 si le voisin n'existe pas 
 */
int voisin(int i, int k, int rs, int nb) {
switch(k)
  {
 case 0: if (i%rs!=rs-1)              return i+1;    else return -1;break;
 case 1: if ((i%rs!=rs-1)&&(i>=rs))   return i+1-rs; else return -1;break;
 case 2: if (i>=rs)                   return i-rs;   else return -1;break;
 case 3: if ((i>=rs)&&(i%rs!=0))      return i-rs-1; else return -1;break;
 case 4: if (i%rs!=0)                 return i-1;    else return -1;break;
 case 5: if ((i%rs!=0)&&(i<nb-rs))    return i-1+rs; else return -1;break;
 case 6: if (i<nb-rs)                 return i+rs;   else return -1;break;
 case 7: if ((i<nb-rs)&&(i%rs!=rs-1)) return i+rs+1; else return -1;break; 
 default:  fprintf(stderr, "bad index value %d\n", k); break;
  }
} // voisin()



/*
 * Label components of an image with ordered labels from 0 to the total nb of different components-1
 * Returns an image where each component is labeled with an index
 */
void label_components(int * image, int width, int height, int **LABEL) { 
/* le LABEL initialement est mis a -1 (minimum) */
  
  int N = width*height;
  int x,y,w,k;
  for (x = 0; x < N; x++) (*LABEL)[x] = -1;
  int nlabels=-1;
  
  Lifo *LIFO= CreeLifoVide(N);//OK
  
  for (x = 0; x < N; x++) {
    if ((*LABEL)[x] == -1) {
      nlabels += 1;
      (*LABEL)[x] = nlabels;
      LifoPush(LIFO, x);
	while (!LifoVide(LIFO)) {
          w = LifoPop(LIFO);
	  for (k = 0; k < 8; k += 1) { /* parcourt les voisins en 8-connexite */
	    y = voisin(w, k, width, N);
	    if ((y != -1) && ((*LABEL)[y] == -1) && (image[y] == image[w]))
	      {
		(*LABEL)[y] = nlabels;
		LifoPush(LIFO, y);
	      } /* if y ... */
	  } /* for k ... */
	} /* while (! LifoVide(LIFO)) */
    } /* if ((*LABEL)[x] == 0) */
  } /* for (x = 0; x < N; x++) */
  LifoTermine(LIFO);  
  // return LABEL;
}

/*
 * Label components of an image
 * Returns an image where each component is labeled with an index
 * Removes regions of small size
 */
void remove_small_components(int * im, /* in : initial components segmentation */
			     int width, // in : image width 
			     int height, // in : image height
			     int minsize, // in :  minimum component size (enforced by post-processing stage) 
			     rgb*colors,/* initial color map*/
			     image<rgb> **output, /* output color map*/
			     int ** size, // in-out: array containing each region size 
			     int * num_ccs,  // out : nb of regions in the new segmentations
			     int ** LABEL,  /* out : resulting segmentation S_{t+1}*/
			     uchar *semantic, // in : previous semantic result
			     image<uchar> **semantic_result //out : current semantic result 
			     ) { 
  
  int N = width*height;  
  int x,y,w,k, cx,cy;
  for (x = 0; x < N; x++) {
    (*LABEL)[x] = -1;
     (*size)[x]=0;
  }
  int nlabels=-1;
  
  Lifo *LIFO= CreeLifoVide(N);
  
  // first pass to compute the size of each component
  for (x = 0; x < N; x++) {
    if ((*LABEL)[x] == -1) {
      nlabels += 1;
      (*LABEL)[x] = nlabels;
      LifoPush(LIFO, x);
	while (!LifoVide(LIFO)) {
          w = LifoPop(LIFO);
	  for (k = 0; k < 8; k += 2) { /* neighbors in 4-connectivity */
	    y = voisin(w, k, width, N);
	    if ((y != -1) && ((*LABEL)[y] == -1) && (im[y] == im[w]))
	      {
		(*LABEL)[y] = nlabels;
		LifoPush(LIFO,y);
                (*size)[nlabels]++;
	      } /* if y ... */
	  } /* for k ... */
	} /* while (! LifoVide(LIFO)) */
    } /* if ((*LABEL)[x] == 0) */
  } /* for (x = 0; x < N; x++) */ 


  // second pass to update size, label, and color

  for (x = 0; x < N; x++) {
     cx = x%width;
     cy=x/width;
     if ((*size)[(*LABEL)[x]]<minsize) {
       // identify one neighbor of an other label
       for (k = 0; k < 8; k += 2) {
           y = voisin(x, k, width, N);
           if ((y!=-1) && ((*size)[(*LABEL)[y]]>=minsize)) {
             // merge
             (*size)[(*LABEL)[x]] --;
             (*size)[(*LABEL)[y]]= (*size)[(*LABEL)[y]]+ 1;
             (*LABEL)[x]= (*LABEL)[y];
	     im[x]= im[y];
	     break;
           } 
       }
     }
     imRef((*output),cx,cy)= colors[im[x]];
     imRef((*semantic_result), cx, cy)= semantic[im[x]];
  }

  // third pass to range values of LABEL into the range of the exact number of different labels
   
 *num_ccs=-1;
   for (x = 0; x < N; x++) 
    im[x] = -1;
  for (x = 0; x <nlabels ; x++) 
    (*size)[x] = 1;
  
  for (x = 0; x < N; x++) {
    if (im[x] == -1) {
      (*num_ccs) += 1;
      im[x] = (*num_ccs);
      LifoPush(LIFO, x);
	while (!LifoVide(LIFO)) {
          w = LifoPop(LIFO);
	  for (k = 0; k < 8; k += 2) { 
	    y = voisin(w, k, width, N);
	    if ((y != -1) && (im[y] == -1) && ((*LABEL)[y] == (*LABEL)[w]))
	      {
		im[y] = (*num_ccs);
		LifoPush(LIFO, y);
                (*size)[(*num_ccs)]++;
	      } /* if y ... */
	  } /* for k ... */
	} /* while (! LifoVide(LIFO)) */
    } /* if ((*LABEL)[x] == 0) */
  } /* for (x = 0; x < N; x++) */ 
  (*num_ccs) += 1;
  delete [] *LABEL;
  *LABEL = im;

LifoTermine(LIFO); 

}


#endif


