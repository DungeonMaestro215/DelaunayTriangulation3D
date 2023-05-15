/* d3 Delaunay triangulator in 3d             Jack Snoeyink Aug 2003
Implements Watson's incremental Delaunay, with sphere-based search and
simplical data structure storing vertices and opposite neighbors.
Handles degeneracies by perturbing points by increasing infinitesimals.  
Guaranteed for integer coordinates of 10 bits. No flat simplices.
*/

#include <math.h>

/* d3.c Delaunay/Power diagram function 
d3batch takes input vertices, and returns a (compact) corner table for Delaunay.
REQUIRES that the first point is at infinity, 
and that the first 5 are in general position. (I should verify or relax this.)
Since all points have radii assigned already, it can produce power diagrams.
*/
void d3batch(ppointType vert, int nvert, // input vertices (pointType[]) & number of vertices
        pcornerType *result, int *ncorners); // output corner table

/* int d3initialize(const pd3stateType this, ppointType vertArray, int nvert) 
initialize d3state this from vertArray.  Allocates memory for corners and spheres, 
sets up the free list for tetrahedra, and creates spheres for the first five points.
REQUIRES: vertArray[0] contains the point at infinity and vertArray[1..4] contain 
four finite points; these first five points must be in general position. 
*/
int d3initialize(const pd3stateType this, ppointType vertArray, int nvert);

// LOCATION ROUTINES walk a mesh stored in s, active, & sph starting from corner start to find a point pv.
// result is the index of a simplex whose sphere strictly contains pv.  
// (The simplex itself need not contain pv.)
// The return code is positive if we succeed (# of location steps at present)
// Return code of 0 means that we failed, perhaps because points with small weight have no Voronoi cells
// Return code of -1 means that we found the duplicate of a vertex in the mesh. IN THIS CASE, result 
//   is the location of the corner where we found the duplicate!
int d3locSphere(const pd3stateType this, ppointType pv, int start,
                    int *result); // output

/* d3insert(const pd3stateType this, int vi, int p) 
   inserts the point this->vert[vi] that is contained in sphere p 
   into the delaunay triangulation stored in this.  
   (d3locSphere may be used to obtain p.)
 */
void d3insert(const pd3stateType this, int vi, int p);

/* d3compactCorners(const pd3stateType this, pcornerType *result, int *ncorners);
   Takes corner table this->s in which some corners/tetrahedra are unused, and 
   returns a compactified corner table result, and its length ncorners.
   DESTROYS the corner table s and active flags in the process.
*/
int d3compactCorners(const pd3stateType this, // arrays this->s & this->active are DESTROYED!
                    pcornerType *result, int *ncorners); // output 

// Functions to access corner tables
#define MOD4(a) (a & 3)
#define TETRA(corner) ((corner) >> 2)
#define INDEX(corner) (MOD4(corner))
#define CORNER(tetra,index) (((tetra)<<2)+(index))
#define BASECORNER(corner) ((corner) & 0xFFFFFFFC)
#define LASTCORNER(corner) ((corner)|3)

#define DEAD(p) (this->active[p] <= 0) // is this a dead or killed tetrahedron?
#define KILL(p) {this->active[p] = -1;} // kill tetrahedron
//#define DEAD(p) (sph[p].sq < 0) // is this a dead tetrahedron?
//#define KILL(p) {sph[p].sq = -1;} // kill tetrahedron

#define infiniteV(pv,vert) ((pv) == vert) // first point is at infinity
#define infiniteC(c) infiniteV(this->s[c].v, this->vert) // corner uses inf pt
#define infiniteP(p) (this->sph[p].sq == 0.0) // if tetra uses infinite point

/* Tetrahedra are groups of four corners in order of increasing 
vertex index, except that the first two may be swapped to ensure 
that the orientation determinant is positive.  
I.e., take the lex smallest alternating permutation with positive sign. 
There must always be an odd number of swaps between two permutations; 
we swap the first two if necessary to achieve this. 
*/

// set corner's vertex and opposite in tetrahedron structure
//void setCornerVC(int c, int vv, int op) {
#define setCornerVC(c, vv, op) { this->s[c].v = vv; this->s[c].opp = op;}
#define setCornerPairV(c, op, vv, ov) { setCornerVC(c, vv, op); setCornerVC(op, ov, c); }
#define setCornerVCN(c, vv, op) { setCornerVC(c, vv, op); this->s[op].opp = c; } // set corner & adjust nhbr opp
