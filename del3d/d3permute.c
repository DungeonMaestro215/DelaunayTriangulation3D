/* d3permute.c                Jack Snoeyink July 2003
Using the various readers (pdbreader, plyreader, etc), 
we read in coordinates and compute bounding boxes
and levels.  */

#include <time.h>
#include "d3permute.h"

// boxOrder gives the order of boxMatrix elements for the Hilbert curve generators associated 
// with the 12 different edges.  hilbertCase gives the generator number for boxOrder[0].
#include "d3boxOrder.c"

//I'm being lazy here and using some globals for coordinate information
bboxType bb;  /* global: assume that we open one file at a time. */
#ifndef NOCOORDHIST
coordHistType cHist; /* global: " */
#endif

boxMatrixType bm, bm2;  /* global: two box levels max */


#define NLEVELS 7 
// NLEVELS-1 - # of trailing zeros is level number 
static const int levelLUT[0x40] = {0,6,5,6,4,6,5,6,3,6,5,6,4,6,5,6,
				   2,6,5,6,4,6,5,6,3,6,5,6,4,6,5,6,
				   1,6,5,6,4,6,5,6,3,6,5,6,4,6,5,6,
				   2,6,5,6,4,6,5,6,3,6,5,6,4,6,5,6};

static const int high[32] = {0, 1, 2,2, 3,3,3,3, 4,4,4,4,4,4,4,4, 
			     5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};  // for BBbits
		      
inline int BBbits(bboxType bb) { // returns the number of bits in bb range
  int range;
  int mask = 0xffffFFE0; // mask last 5 bits
  int shift = 0;
  range = bb[1]-bb[0]; // find max range in all three coords
  if (range < bb[3]-bb[2]) range = bb[3]-bb[2];
  if (range < bb[5]-bb[4]) range = bb[5]-bb[4];

  while (range & mask) { // some ones bits left above last 5.
    mask <<= 5;
    shift += 5;
  }
  return shift+high[range>>shift]; // shift needed to make all bits zero
}

// This function permutes points from v to newv by boxOrder
// From a coordinate x, it uses bits (x>>bitshift)&mask.  nboxes = (mask+1)^3;
// hcase is so we can handle hilbert cases for first level. 
// This function is here for documentation.  We don't actually call it, 
//  but weave this code into the reading, tuning it for efficiency.
//
/*
inline void permute(boxMatrixType bm, const ppointType v, const int nvert, ppointType newv, 
		    const short bitshift, const short mask, const short nboxes, const short hcase) {
  pboxMatrixEntryType pbm; 
  ppointType pv; 
  int i, n, tmp; 

  bzero((void *)(bm), sizeof(bm)); // init box counts 
  for (pv = v; pv < v+nvert; pv++) 
    boxIncrv(bm, pv, mask); // count for each box 

  n = 0; // prefix sum the boxes in boxOrder 
  for (i = 0; i < nboxes; i++) { 
    pbm = ((pboxMatrixEntryType) bm) + boxOrder[hcase][i]; 
    tmp = *pbm; 
    *pbm = n; 
    n += tmp; 
  }

  ASSERT(n == nvert, "total wrong after bm assignment");
  for (pv = v; pv < v+nvert; pv++) 
    newv[boxIncrv(bm, pv, mask)] = *pv; // move to place 
}
/**/
void bbPrint() {
  int i;
  printf("BB %d (%f %f %f; %f %f %f)\n", BBbits(bb),
	 bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
#ifndef NOCOORDHIST
  for (i = 0; i<NUMTAIL2; i++) 
    printf("%5d %5o  %4d %4d %4d\n", i, i, cHist[XX][i], cHist[YY][i], cHist[ZZ][i]);
#endif
}

// Do this once to a point, because it assumes point has been lifted by -rad^2
// move origin and lift (rad already there)
#define translateAndLift(v, xof, yof, zof) { \
     (v)->x -= xof; (v)->y -= yof; (v)->z -= zof; (v)->sq += DOT(v,v); }

#define NBOXES 64
// read vertices & return coordinfo
void d3permute(FILE *fid, int (*reader)(FILE *, ppointType), 
	    ppointType *vout, int *nvert) {
  int i, j, k, tmp;
  int bitshift, hcase;
  int n, nv, nvi;
  int nblock, sblock, eblock;
  int origx, origy, origz;
#ifndef NOCOORDHIST
  int xorBits[3];
#endif
  clock_t tic, toc; /**/
  ppointType vin, v; 
  ppointType pv, pvi;
  pboxMatrixEntryType pbm;
  int nlev[NLEVELS]; // Number/offset per level 

  bbInit(); // initialize bounding box
  cHistInit(); // initialize coordinate histogram
  vin = (ppointType) calloc(MAXVERT, sizeof(pointType)); // pdb has a five digit atom# field
  if (vin == NULL) { printf("ERROR: could not allocate memory for MAXVERT points\n"); exit(EXIT_FAILURE); }

  nvi = reader(fid, vin); // read in vertices & compute bbox and cHist

  tic = clock();/**/
  if (nvi > MAXVERT) { printf("FATAL ERROR: Too many points %d; increase MAXVERT\n", nvi); exit(EXIT_FAILURE); }
  vin = (ppointType) realloc(vin, (nvi+1)*sizeof(pointType));
  v = (ppointType) calloc(nvi+1, sizeof(pointType)); // one more for infinite point
  if (v == NULL) { printf("ERROR: could not allocate memory for %d points\n", nvi+1); exit(EXIT_FAILURE); }
  v[0].index = -1; v[0].x = v[0].y = v[0].z = 0; v[0].sq = 1; 
  *nvert = nvi+1; // output results
  *vout = v;
  
#ifndef NOCOORDHIST
  for (k = XX; k<= ZZ; k++) { // i,i+1 are bit patterns we decide between
    i = (MASKTAIL-1)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; // 1 bit 
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; // 2 bit 
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; // 3 bit 
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; // 4 bit 
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; // 5 bit 
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; // 6 bit 
    xorBits[k] = i; // frequent bit pattern for this coordinate
  }
#endif

  bzero((void *)(nlev), sizeof(nlev)); // zero, then accumulate level counts for vin
  for (pvi = vin; pvi < vin+nvi; pvi++) 
    nlev[LEVELNO(pvi)]++;

  n = 1; 
  for (i = 0; i<NLEVELS; i++) { // prefix sum level counts for level offsets 
    tmp = nlev[i]; 
    nlev[i] = n; 
    n += tmp; 
  }
  ASSERT(n==1+nvi, "bad count after level number prefix sum");

  for (pvi = vin; pvi < vin+nvi; pvi++) { // copy levels to v; adjust coords to 0--range
    pv = v + (nlev[LEVELNO(pvi)]++); 
    VSUBASSN(pvi, bb[0], bb[0], bb[0], pv);
  }    

  bitshift = BBbits(bb)-3; // find how many bits we need for entire range

  nv = 1; // index of next point to consider in v
   printf("Levels ");
  for (j = 0; j < NLEVELS; j++) { // order each level
    /*    printf("Level %d has %d [%d,%d): ", j, nlev[j]-nv, nv,nlev[j]);  /**/
    printf("%d:%d ", j, nlev[j]-nv);  /**/
    if (nlev[j]-nv > 64) { // 512 boxes on this level:  v[nv..nlev[j]) permuted using vin[nv..nlev[j])

      bzero((void *)(bm), sizeof(bm)); /* init box counts */
      for (pv = v+nv; pv < v+nlev[j]; pv++) 
	boxIncrv(bm, pv, 7); // count vertices of v[nv..nlev[j]) in each box 
      n = nv; 

      for (i = 0; i < 512; i++) { /* prefix sum the boxes in boxOrder */
	/*printf("<%d:%o ", n, i);/**/
	pbm = ((pboxMatrixEntryType) bm) + boxOrder[j&1][i]; // alternating directions on levels
	tmp = *pbm; *pbm = n; n += tmp; 
      }
      /*printf("<%d]\n", n);/**/

      ASSERT(n == nlev[j], "total wrong after bm assignment");
      /*      printf("\n L%d:"); /**/
      for (pv = v+nv; pv < v+nlev[j]; pv++) {
	/*	printf(" %d>%d,", pv->index, boxv(bm, pv, 7)); /**/
	vin[boxIncrv(bm, pv, 7)] = *pv; //  move to place in vin 
      } 
      // we now have this level in vin[nv..nlev[j]), and need it back in v[nv..nlev[j])
      // groups with few points can just be copied; those with many get radix sorted.
      // bm[] contains the prefix sum of counts INCLUDING current. I.e. old bm[1] is now bm[0]

      // We now consider radix sorting blocks in this level. 
      // invariants: sblock = start of curr block; nvi = start of block that needs to be copied to v
      // i is block number, eblock is end of current block, nblock is # in block
      eblock = nvi = nv;
      for (i = 0; i < 512; i++) { // undo prefix sum (shifted) to get number in block
	sblock = eblock;
	eblock = ((pboxMatrixEntryType) bm)[boxOrder[j&1][i]];
	nblock = eblock - sblock; 
	if (nblock > 32){ // shuffle this block 
	  if (nvi < sblock) // There are old blocks to copy first
	    memcpy((void *)(v+nvi), (const void *)(vin+nvi), (sblock-nvi)*sizeof(pointType));

	  /*      printf("\n   c(%d,%d)",nvi, sblock); /**/
	  bitshift -= 3; 
	  hcase = hilbertCase[j&1][i];
	  // Now we shuffle and copy vin[sblock..eblock) to v[sblock..eblock), with permutation

	  bzero((void *)(bm2), sizeof(bm2)); /* init box counts */
	  for (pvi = vin+sblock; pvi < vin+eblock; pvi++) 
	    boxIncrv(bm2, pvi, 7); // count vertices in each box 
	  n = sblock; 
	  for (k = 0; k < 512; k++) { /* prefix sum the boxes in boxOrder */
	    /*printf("{%d:%o|%o ", n, i, k);/**/
	    pbm = ((pboxMatrixEntryType) bm2) + boxOrder[hcase][k]; 
	    tmp = *pbm; 
	    *pbm = n; 
	    n += tmp; 
	  }
	  /*printf("{%d}\n", n);/**/
	  ASSERT(n == eblock, "total wrong after bm assignment");
	  for (pvi = vin+sblock; pvi < vin+eblock; pvi++) {
	    /*	    printf(" %d>%d,", pvi->index, boxv(bm2, pvi, 7)); /**/
	    k = boxIncrv(bm2, pvi, 7);
	    v[k] = *pvi; 
	  }

	  bitshift += 3; // undo bitshift change above
	  nvi = eblock; // Here's where we be after shuffle & copy 
	}
      }
      if (nvi < nlev[j]) // There are still old blocks to copy to finish this level
	memcpy((void *)(v+nvi), (const void *)(vin+nvi), (nlev[j]-nvi)*sizeof(pointType));
    }
    nv = nlev[j];
  }
  fflush(stdout);
  free(vin);

  origx = (bb[3]+bb[0])/2; // center the bbox
  origy = (bb[4]+bb[1])/2;
  origz = (bb[5]+bb[2])/2;
  for (pv = v+1; pv < v+(*nvert); pv++) 
    translateAndLift(pv, origx, origy, origz);

  toc = clock();
  printf("Hilbert time(secs) %f\n", ((double) (toc - tic)) / CLOCKS_PER_SEC);  (void)fflush(stdout); /**/

#ifndef NOASSERT
  {
    short *index = (short *) calloc(*nvert, sizeof(short));
    bzero(index, (*nvert)*sizeof(short));
    for (i = 1; i < (*nvert); i++) 
      if (index[v[i].index] != 0) 
	printf("POST: INDEX %d appears at %d and %d\n", v[i].index, i, index[v[i].index]);      
      else
	index[v[i].index] = i;
    free(index);
  }
#endif
}
