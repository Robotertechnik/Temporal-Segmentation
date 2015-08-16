// color_flow.h 
// color-code motion field
// normalizes based on specified value, or on maximum motion present otherwise
// code from Daniel Scharstein, 7/2/07

#include <stdlib.h>
#include <math.h>

int ncols = 0;
#define MAXCOLS 60
int colorwheel[MAXCOLS][3];
#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))

void setcols(int r, int g, int b, int k)
{
    colorwheel[k][0] = r;
    colorwheel[k][1] = g;
    colorwheel[k][2] = b;
}

void makecolorwheel()
{
    // relative lengths of color transitions:
    // these are chosen based on perceptual similarity
    // (e.g. one can distinguish more shades between red and yellow 
    //  than between yellow and green)
    int RY = 15;
    int YG = 6;
    int GC = 4;
    int CB = 11;
    int BM = 13;
    int MR = 6;
    ncols = RY + YG + GC + CB + BM + MR;
    //printf("ncols = %d\n", ncols);
    if (ncols > MAXCOLS)
	exit(1);
    int i;
    int k = 0;
    for (i = 0; i < RY; i++) setcols(255,	   255*i/RY,	 0,	       k++);
    for (i = 0; i < YG; i++) setcols(255-255*i/YG, 255,		 0,	       k++);
    for (i = 0; i < GC; i++) setcols(0,		   255,		 255*i/GC,     k++);
    for (i = 0; i < CB; i++) setcols(0,		   255-255*i/CB, 255,	       k++);
    for (i = 0; i < BM; i++) setcols(255*i/BM,	   0,		 255,	       k++);
    for (i = 0; i < MR; i++) setcols(255,	   0,		 255-255*i/MR, k++);
}

rgb computeColor(float fx, float fy)
{
  rgb pix;
    if (ncols == 0)
	makecolorwheel();

    float rad = sqrt(fx * fx + fy * fy);
    float a = atan2(-fy, -fx) / M_PI;
    float fk = (a + 1.0) / 2.0 * (ncols-1);
    int k0 = (int)fk;
    int k1 = (k0 + 1) % ncols;
    float f = fk - k0;
    //f = 0; // uncomment to see original color wheel
    for (int b = 0; b < 3; b++) {
	float col0 = colorwheel[k0][b] / 255.0;
	float col1 = colorwheel[k1][b] / 255.0;
	float col = (1 - f) * col0 + f * col1;
	if (rad <= 1)
	    col = 1 - rad * (1 - col); // increase saturation with radius
	else
	    col *= .75; // out of range
	if (b==0) pix.b = (int)(255.0 * col);
	if (b==1) pix.g = (int)(255.0 * col);
	if (b==2) pix.r = (int)(255.0 * col);
	//pix[2 - b] = (int)(255.0 * col);
    }
    return pix;
}


void MotionToColor(int *ftabx, int * ftaby, image<rgb> **output)
{

  int width = (*output)->width();
  int height = (*output)->height();
    int x, y;
    // determine motion range:
    float maxx = -999, maxy = -999;
    float minx =  999, miny =  999;
    float maxrad = -1;
    for (y = 0; y < height; y++) {
	for (x = 0; x < width; x++) {
	    float fx = ftabx[x+ y*width];
	    float fy = ftaby[x+ y*width];
	    
	    maxx = __max(maxx, fx);
	    maxy = __max(maxy, fy);
	    minx = __min(minx, fx);
	    miny = __min(miny, fy);
	    float rad = sqrt(fx * fx + fy * fy);
	    maxrad = __max(maxrad, rad);
	}
    }
    //printf("max motion: %.4f  motion range: u = %.3f .. %.3f;  v = %.3f .. %.3f\n",
    //	   maxrad, minx, maxx, miny, maxy);
    if (maxrad == 0) // if flow == 0 everywhere
	maxrad = 1;

    for (y = 0; y < height; y++) {
	for (x = 0; x < width; x++) {
	    float fx = ftabx[x+ y*width];
	    float fy = ftaby[x+ y*width];
            rgb pix;
	    pix = computeColor(fx/maxrad, fy/maxrad);
	    imRef((*output), x, y)= pix;
	}
    }
}

