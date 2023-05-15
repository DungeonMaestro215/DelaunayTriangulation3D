/* d3permute.h                Jack Snoeyink July 2003
Using the various readers (pdbreader, plyreader, etc), 
we read in coordinates and compute bounding boxes
and levels */

#include "delaunay3.h"

/* 
d3permute uses a  reader (pdbreader, plyreader, etc) to obtain coordinates from a file, 
then permutes them using a Hilbert curve into a good order for incremental Delauany.
It adds a 0th vertex at infinity. 
*/
void d3permute(FILE *fid, int (*reader)(FILE *, ppointType), ppointType *pv, int *nvert); 

#define MASKTAIL  0x3F // mask for bits in tail (for levels)
#define NUMTAIL   0x40 // how many bit patterns in tail
#define NUMTAIL2  0x80 // 2*NUMTAIL

// coordinate info to record while reading
typedef int bboxType[6];

// Bounding box & coord histogram functions
// WARNING: I'm being lazy here, and using globals
extern bboxType bb;  /* global: assume that we open one file at a time. */
#define BBMAXCOORD 0xFFFffff // 28 bits should be enough
#define bbInit() { bb[0] = bb[1] = bb[2] =  BBMAXCOORD; /*min x,y,z*/\
                   bb[3] = bb[4] = bb[5] = -BBMAXCOORD; /*max x,y,z*/}

// -DNOCOORDHIST will turn of the use of coordinate tail histograms for choosing levels.
#ifdef NOCOORDHIST
#define cHistInit() 
#define LEVELNO(pv) levelLUT[((((int)(pv->x)) | ((int)(pv->y)) \
 	            | ((int)(pv->z)) ) & MASKTAIL)]
#else
typedef short coordHistType[3][NUMTAIL2]; // histograms for last 7 coord bits
extern coordHistType cHist; /* global: assume that we open one file at a time. */
#define cHistInit() { bzero((void *)cHist, sizeof(cHist)); }
#define LEVELNO(pv) levelLUT[(((xorBits[XX]^((int)(pv->x))) | (xorBits[YY]^((int)(pv->y))) \
 	            | (xorBits[ZZ]^((int)(pv->z))) ) & MASKTAIL)]
#endif

// Called by readers to set coordinates. WARNING: uses global variables bb and cHist 
#define TWO26 0x2000000 // 2^26 added to coords to make them positive
inline void setVert(ppointType v, int indx, double xx, double yy, double zz, 
	     double rad, double mult) {  // this part of the code may be I/O bound anyway, 
  double mr = mult*rad;                  // so we count most common coordinate bit tails
  int ix, iy, iz, i;
  (v)->index = indx;
  (v)->x = ix = (int)(mult*(xx)+0.5); 
  (v)->y = iy = (int)(mult*(yy)+0.5); 
  (v)->z = iz = (int)(mult*(zz)+0.5); 
  (v)->sq = (int)(-mr*mr);

  if (bb[0] > ix) bb[0] = ix; // update bounding box
  if (bb[1] > iy) bb[1] = iy; 
  if (bb[2] > iz) bb[2] = iz; 
  if (bb[3] < ix) bb[3] = ix; 
  if (bb[4] < iy) bb[4] = iy; 
  if (bb[5] < iz) bb[5] = iz; 

#ifndef NOCOORDHIST
  i = (ix+TWO26)&MASKTAIL; cHist[XX][i]++; // add to coordinate histograms: 6 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; // 5 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; // 4 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; // 3 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; // 2 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; // 1 bit

  i = (iy+TWO26)&MASKTAIL; cHist[YY][i]++; // add to coordinate histograms: 6 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; // 5 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; // 4 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; // 3 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; // 2 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; // 1 bit

  i = (iz+TWO26)&MASKTAIL; cHist[ZZ][i]++; // add to coordinate histograms: 6 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; // 5 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; // 4 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; // 3 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; // 2 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; // 1 bit
#endif
}


// permutation by boxOrder
// box indices
#define BX(x,mask) (((x)>>bitshift)&mask) 
typedef int boxMatrixType[8][8][8]; // a boxMatrix has an int for each 3 bits of x,y,z
typedef int *pboxMatrixEntryType; // pointer to a box matrix entry

//increment boxMatrix entry (based on bitshift & mask)
#define boxv(bm, v, mask) (bm[BX((int)((v)->x),mask)][BX((int)((v)->y),mask)][BX((int)((v)->z),mask)])/**/
#define boxIncr(bm, x,y,z, mask) (bm[BX(x, mask)][BX(y, mask)][BX(z, mask)]++)
#define boxIncrv(bm, v, mask) boxIncr(bm, (int)((v)->x), (int)((v)->y), (int)((v)->z), mask)
