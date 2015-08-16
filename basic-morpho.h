// Original code from Michel Couprie

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <imutil.h>
#include <image.h>
#include <misc.h>
#include "list_utils.h"

#define INFINI 1000000000

int32_t lsedt_meijster(image<uchar>*img, uint32_t **res);

/* ==================================== */
void inverse(image<uchar>* I)
/* ==================================== */
{
  for (int y = 0; y < I->height(); y++) {
    for (int x = 0; x < I->width(); x++) {
      
      if ( imRef(I,x,y) == 0) 
	imRef(I,x,y) = 255;
	else imRef(I,x,y) = 0;
    }
  }
}




void boundary(int * img, int rs, int cs, image<uchar>* out)
//return all contours (white on black) present in the image I 
{
  for (int i=2; i<rs-1;i++) {
    for (int j=2; j<cs-1;j++) {
      if ((img[i-1+j*rs]!=img[i+j*rs]) || (img[i+1+j*rs]!=img[i+j*rs])  || (img[i+(j-1)*rs]!=img[i+j*rs]) || (img[i+(j+1)*rs]!=img[i+j*rs]))
	imRef(out, i,j) = 255;
      else imRef(out, i,j) = 0;
    }
  }
}
    
      
/* ==================================== */
int ldilatdisc(image<uchar>* ob, int32_t r, bool erode)
/* ==================================== */
// dilation by a disc of radius r
{
  uint32_t rs = ob->width();
  uint32_t cs = ob->height();
   
  unsigned long i, N = rs*cs;
  uint32_t *dist = new uint32_t [N];

  if (erode==false) inverse(ob);
  if (!lsedt_meijster(ob, &dist)) 
    return 0;
  
  for (int y = 0; y < cs; y++) {
    for (int x = 0; x < rs; x++) {
      i = y*rs+x;
      if (dist[i] > r)
	imRef(ob, x, y) = 0;
      else
	imRef(ob, x, y) = 255;
    }
  }
  delete [] dist;    
  return 1;
} // ldilatdisc()



// Functions for the exact Squared Euclidean Distance Transform (Meijster et al., Linear algorithm)
#define F_2_2d(y,yp,f,i) (f[y*rs+i]+(yp-y)*(yp-y))
#define Sep_2_2d(v,u,f,i) (((u*u)-(v*v)+f[u*rs+i]-f[v*rs+i])/(2*(u-v)))

/* ======================================================== */
void SEDT_line(image<uchar> *f, uint32_t *g, uint32_t rs, uint32_t cs)
/* ======================================================== */
{
  int32_t i, j;
  for (j = 0; j < (int32_t)cs; j++)
  {
    if (imRef(f, 0 ,j) == 0) g[0 + rs*j] = 0; else g[0 + rs*j] = rs*cs; // infinity
    for (i = 1; i < (int32_t)rs; i++) {
      if (imRef(f, i, j) == 0) g[i + rs*j] = 0; 
      else                  g[i + rs*j] = 1 + g[i-1 + rs*j]; 
    }
    for (i = rs-2; i >= 0; i--)
      if (g[i+1 + rs*j] < g[i + rs*j]) g[i + rs*j] = 1 + g[i+1 + rs*j];
    for (i = 0; i < (int32_t)rs; i++) {
      if (g[i + rs*j] < rs*cs) // NECESSAIRE pour éviter un overflow
	g[i + rs*j] = g[i + rs*j] * g[i + rs*j];
    }
  }
} //  SEDT_line()

/* ======================================================== */
void SEDT_column(uint32_t *f, uint32_t **g, uint32_t rs, uint32_t cs)
/* ======================================================== */
{
  int32_t i, u, q;
  int32_t w;
  uint32_t *s, *t;
  s = (uint32_t *)calloc(1,cs * sizeof(uint32_t));
  t = (uint32_t *)calloc(1,cs * sizeof(uint32_t));

  for (i = 0; i < (int32_t)rs; i++) {
    q = 0; s[0] = 0; t[0] = 0;
    for (u = 1; u < (int32_t)cs; u++) {
      while ( (q >= 0) && (F_2_2d(s[q],t[q],f,i) > F_2_2d(u,t[q],f,i)) ) q--;
      if (q < 0) {
        q = 0;
        s[0] = u;
      }
      else {
        w = 1 + Sep_2_2d(s[q],u,f,i);
        if (w < (int32_t)cs)
        {
          q++; s[q] = u; t[q] = w;
        }
      }
    } 
    for (u = cs-1; u >= 0; u--) {
      (*g)[rs*u + i] = F_2_2d(s[q],u,f,i);
      if (u == (int32_t)(t[q])) q--;
    }
  }
  free(s); free(t);
} //  SEDT_column()



/* ==================================== */
int32_t lsedt_meijster(image<uchar>*img,   /* data: binary image */       
		       uint32_t **res    /* result: distances */
		       )
/* ==================================== */
//Calls the SEDT linear algorithm (Meijster & al.)
{ 
  uint32_t rs = img->width();
  uint32_t cs = img->height();
   
  uint32_t * tmp = new uint32_t [rs*cs];//allocimage(NULL, rs, cs, ds, VFF_TYP_4_BYTE);
    
  SEDT_line(img, tmp, rs, cs);
  SEDT_column(tmp, res, rs, cs);
    
  delete [] tmp;
  return(1);
} // lsedt_meijster()
