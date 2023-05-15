/* delaunay3 Delaunay/Power diagrams in 3d Jack Snoeyink Aug 2003
Implements Watson's incremental Delaunay, with sphere-based search and
simplical data structure storing vertices and opposite neighbors.
Points must be scaled to integers; guaranteed for differences of 10
bits, and uses leveling and hilbert curve to guarantee 16 bits if the
points are well distributed.  Works for pdb files, which are 20 bits.
Handles degeneracies by perturbing points by increasing infinitesimals
to guarantee that all simplices are full-dimensional.
*/

#ifndef DELAUNAY3_H
#define DELAUNAY3_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef float coord; // We'll just use floats for x,y,z coordinates

typedef struct { // Basic point data structure element
  int index; // index back to original line in file.
  coord x,y,z; double sq;  // coordinates: 3 spatial + lifted
} pointType, *ppointType;

//tetrahedra are groups of four consecutive corners
typedef struct cornerType { // data associated with each corner 
  ppointType v; // index of vertex
  int opp; // pointer to opposite corner in neighboring tetra
} cornerType, *pcornerType;

/* Tetrahedra are groups of four corners in order of increasing 
vertex index, except that the first two may be swapped to ensure 
that the orientation determinant is positive.  
I.e., take the lex smallest alternating permutation with positive sign. 
There must always be an odd number of swaps between two permutations; 
we swap the first two if necessary to achieve this.  */

/* d3.c Delaunay/Power diagram function 
d3batch takes input vertices, and returns a (compact) corner table for Delaunay.
REQUIRES that the first point is at infinity, 
and that the first 5 are in general position. 
(I should probably verify or relax this.)
Since all points have radii assigned already, it can handle power diagrams.
*/
void d3batch(ppointType vert, int nvert, // input vertices (pointType[]) & number of vertices
	pcornerType *result, int *ncorners); // output corner table


#define XX 0
#define YY 1
#define ZZ 2

#define NVERT     5000 // number of random vertices
#define MAXVERT 200000 // pdb has 5 digit atom# field
#define TETperV    10 // 
#define FORGETLIMIT  10000 // forget tetra after this many (for better locality of ref? Doesn't help.)

// Some useful definitions
#define EQUALV(i,j) (vert[i].x == vert[j].x && vert[i].y == vert[j].y && vert[i].z == vert[j].z)
#define EQUALPV(pv,v) (pv->x == v->x && pv->y == v->y && pv->z == v->z)

#define DET2(p,q,i,j)((i##p)*(j##q) - (j##p)*(i##q))
#define DOT(a,b) ((a)->x*(b)->x + (a)->y*(b)->y + (a)->z*(b)->z)

// Stack data structure operations
#define STACKMAX 2000
#define POP(stack) (stack##st[stack##sp--])
#define isEMPTY(stack) (stack##sp < 0)
#define stkINIT(stack) {stack##sp = -1; }

#ifndef STATS
#define PUSH(value, stack) { stack##st[++stack##sp] = value; }
#define stkDECLARE(stack,stn) int stack##sp, stack##st[STACKMAX];    
#else
#define PUSH(value, stack) { \
    stack##st[++stack##sp] = value; \
    if (stack##max < stack##sp) { stack##max = stack##sp; \
      if (stack##max >= STACKMAX) { \
        printf("ERROR: overflow stack %x pushing %d",  stack##st, value); exit(EXIT_FAILURE); } } /**/\
    }
#define stkDECLARE(stack,stn) int stack##sp, stack##st[STACKMAX]; int stack##max; //AUDIT /**/ 
#endif

typedef struct sphereType { // sphere equation
  double x, y, z, sq; // Invariant sq > 0 for all created tetra, unless they use pt at infty
 } sphereType, *psphereType;

typedef struct {
  ppointType vert; // vertices: 0th is point at infinity!!!!
  pcornerType s;  // corner table
  psphereType sph; // spheres
  int *active; // which spheres are active
  int freeTetra; // head for free list for tetrahedra kept in opp[CORNER(tetra,0)]
  int liveTetra; // latest tetra; known to be live. 
  int maxTetra; // AUDIT only
  int limitmaxTetra; // limit on # of created tetrahedra, spheres & corners/4.
  // stacks used in inserting pv
  stkDECLARE(dfs, "dfs");  // DFS stack to find dead tetras
  stkDECLARE(idfs, "idfs"); // DFS stack for tetras adj to infinite vertex (>4*30)
  stkDECLARE(nhbr, "nhbr"); // stack for dead corners with live neighbors
  stkDECLARE(kill, "kill"); // stack for base corners of tetras to recycle
} d3stateType, *pd3stateType;

// readers call setVert when they make a point
void setVert(ppointType v, int indx, double xx, double yy, double zz, 
	     double rad, double mult);

#define VSUBASSN(vin, xx,yy,zz, vout) {\
    vout->index = vin->index;  vout->sq = vin->sq;\
    vout->x = vin->x - xx; vout->y = vin->y - yy; vout->z = vin->z - zz;\
  }

//Some compiler/unix variants
#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef strcasecmp
#define strcasecmp(s1,s2)      strcmp(s1,s2)
#define strncasecmp(s1,s2,n)   strncmp(s1,s2,n)
#endif

#define HIGHCOORD 0x4000
#define COORDMASK 0x3fff

#ifdef bcc
#define RANDBIT (random(2))             // one random bit
#define RAND2BIT (random(4))            // two random bits
#define RANDPROB(mask) (random(mask+1)) // random bits masked
#define RANDCOORD (double)random(HIGHCOORD);
#define RANDOM(k) random(k)             // random number 0..k-1
#else
#define RANDBIT (random()&1)            // one random bit
#define RAND2BIT (random()&3)           // two random bits
#define RANDPROB(mask) (mask&random())  // random bits masked
#define RANDCOORD (double)(random()&COORDMASK)
#define RANDOM(k) random()%(k)          // random number 0..k-1
#endif

//  fprintf(stderr, "\nASSERT FAILED (line %d of %s ): %s\n", 
#ifndef NOASSERT
#define ASSERT(bool, string) if (!(bool)) {\
  printf("\nASSERT FAILED (line %d of %s ): %s\n", \
	  __LINE__, __FILE__, string); }
#else
#define ASSERT(bool, string) 
#endif

#endif
