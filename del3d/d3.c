/* d3 Delaunay triangulation in 3d             Jack Snoeyink Aug 2003
   Implements Watson's incremental Delaunay, with sphere-based location, 
   and simplical data structure storing vertices and opposite neighbors.
   Handles degeneracies by perturbing points by increasing infinitesimals.  
   Guaranteed for integer coordinates of 10 bits. No flat simplices.
*/

#include "delaunay3.h"
#include "d3.h"


#define spdot(sp,pv,sv) ((sp)->x*((pv)->x-(sv)->x)+(sp)->y*((pv)->y-(sv)->y)\
         +(sp)->z*((pv)->z-(sv)->z)+(sp)->sq*((pv)->sq-(sv)->sq))
#define spdotInf(sp,pv) (sp)->sq


/* Allocate or free space for a tetrahedron.  
   When allocating, this->liveTetra is the location fo the new tetrahedron. 
*/
#define STARTTETRA(this) \
{ startTetraCnt++; (this)->liveTetra = (this)->freeTetra; \
  (this)->freeTetra = (this)->s[CORNER((this)->liveTetra,0)].opp;  \
  ASSERT(DEAD((this)->liveTetra), "Reusing existing tetrahedron?"); (this)->active[(this)->liveTetra] = 1; \
  if ((this)->maxTetra <= (this)->liveTetra) { (this)->maxTetra = (this)->liveTetra+1; \
  if ((this)->maxTetra >= (this)->limitmaxTetra) { printf("AUDIT: %d > limitmaxTetra\n", (this)->liveTetra); exit(EXIT_FAILURE); }}/**/\
} // AUDIT 
 
#define FREETETRA(this, p) \
{ freeTetraCnt++;  /*TESTING*/ ASSERT((this)->active[p]<0, "Freeing already free tetrahedron?"); /**/\
  (this)->active[p]=0;  if (((this)->maxTetra-p)<FORGETLIMIT) {\
  (this)->s[CORNER(p,0)].opp = (this)->freeTetra; (this)->freeTetra = p; /**/  \
}}

//Tetrahedron manipulation tables
static const short offset[4][4] = { // c+offset[i][INDEX(c)] advances c to (c+i)mod5
  /*0*/{0,0,0,0}, /*1*/{1,1,1,-3}, /*2*/{2,2,-2,-2}, /*3*/{3,-1,-1,-1}};
#define INCREMENT(c) (c+offset[1][INDEX(c)])
// drop[i] contains new vertex order after vertex i is dropped and replaced by pv on same side.
// offdr[i] contains drop(i)-index(i)
// invdrop[i][k] = j whenever drop[i][j] = k.  4s signal i=k; bad because i is dropped. 
static const short drop[4][3] = {{2,1,3}, {0,2,3}, {1,0,3}, {0,1,2}};
static const short offdr[4][3] = {{2,1,3}, {-1,1,2}, {-1,-2,1}, {-3,-2,-1}};
static const short invdrop[4][4] = {{4,1,0,2}, {0,4,1,2}, {1,0,4,2}, {0,1,2,4}};

#include "d3audit.c"

const double TWO25 = 1024.0*1024.0*32.0;

inline double fm(double x) {
  double y = fmod(x, TWO25);
  return y<0 ? y+TWO25: y;
}

inline double InSpherev(psphereType sp, ppointType pv, ppointType sv) {
  double d = spdot(sp, pv, sv); // Return true if inside (==negative)
#ifndef NOASSERT
  if ((fabs(d/((sp)->x)*((pv)->x-(sv)->x)) < 1e-10) // massive cancellation
      || (fabs(d/((sp)->y)*((pv)->y-(sv)->y)) < 1e-10)
      || (fabs(d/((sp)->z)*((pv)->z-(sv)->z)) < 1e-10)
      || (fabs(d/((sp)->sq)*((pv)->sq-(sv)->sq)) < 1e-10)) {

  double dmod = fm(fm((sp)->x)*((pv)->x-(sv)->x)+fm((sp)->y)*((pv)->y-(sv)->y)+fm((sp)->z)*((pv)->z-(sv)->z)+fm((sp)->sq)*((pv)->sq-(sv)->sq));
  if (fm(d) != dmod) 
    printf("%% %10.0lf != %10.0lf: %18.0lf - %10.0lf = InSphere(<%8.0f %8.0f %8.0f %14.0f>*( %5.0f %5.0f %5.0f %10.0f))\n", 
	   fm(d), dmod, d, fm(d)-dmod, sp->x,sp->y,sp->z,sp->sq, pv->x-(sv)->x, pv->y-(sv)->y,pv->z-(sv)->z,pv->sq-(sv)->sq);
  }
#endif
  inSphereCnt++; /*TESTING*/
  return d; // perturb those on sphere to inside

}

inline void makeSphereV(psphereType sp, ppointType v0, ppointType v1, ppointType v2, ppointType pv, ppointType vert) {
  double x0, y0, z0, sq0, x1, y1, z1, sq1, x2, y2, z2, sq2;
  double xy, xz, xs, yz, ys, zs; // 2x2 minors
  // make sphere: only v0 or v1 may be infinte. 
  sphereCnt++; /*TESTING*/
  if(!infiniteV(v0, vert)) {
    x0 = v0->x - pv->x; y0 = v0->y - pv->y; 
    z0 = v0->z - pv->z; sq0 = v0->sq - pv->sq;
  } else { 
    x0 = v0->x; y0 = v0->y; 
    z0 = v0->z; sq0 = v0->sq; 
  }

  if(!infiniteV(v1, vert)) {
    x1 = v1->x - pv->x; y1 = v1->y - pv->y; 
    z1 = v1->z - pv->z; sq1 = v1->sq - pv->sq;
  } else { 
    x1 = v1->x; y1 = v1->y; 
    z1 = v1->z; sq1 = v1->sq; 
  }

  x2 = v2->x - pv->x; y2 = v2->y - pv->y; 
  z2 = v2->z - pv->z; sq2 = v2->sq - pv->sq;
  xy = DET2(0,1,x,y); 
  xz = DET2(0,1,x,z); 
  yz = DET2(0,1,y,z); 
  xs = DET2(0,1,x,sq); 
  ys = DET2(0,1,y,sq);
  zs = DET2(0,1,z,sq);
  sp->x  = -y2*zs +z2*ys -sq2*yz;
  sp->y  =  x2*zs -z2*xs +sq2*xz;
  sp->z  = -x2*ys +y2*xs -sq2*xy;
  sp->sq =  x2*yz -y2*xz +z2*xy;
  //  sp->w  = -p2->x*sp->x -p2->y*sp->y -p2->z*sp->z -p2->sq*sp->sq; 
  /*  printf("disp('Sphere equation: <%5.0f %5.0f %5.0f %5.0f %5.0f>')\n",sp->x,sp->y,sp->z,sp->sq); 
      (void)fflush(stdout);
      printf("%%Sp ck: %g %g %g %g\n", spdot(sp,v0), spdot(sp,v1), spdot(sp,v2), spdot(sp,pv)); /**/
  }


// when we initialize, this is what we fill in. 
//const int initialopp[] = {11,5,15,21,25, 1,10,16,20,26, 6,0,17,22,27, 2,7,12,23,28, 8,3,13,18,29, 4,9,14,19,24};
static const int initialopp[5][4] = {
  {CORNER(1,1), CORNER(2,0), CORNER(3,1), CORNER(4,0)}, 
  {CORNER(2,1), CORNER(0,0), CORNER(3,0), CORNER(4,1)},  
  {CORNER(0,1), CORNER(1,0), CORNER(3,2), CORNER(4,2)}, 
  {CORNER(1,2), CORNER(0,2), CORNER(2,2), CORNER(4,3)}, 
  {CORNER(0,3), CORNER(1,3), CORNER(2,3), CORNER(3,3)}};

static const int initialv[2][5][4] = {{{1,2,3,4}, {2,0,3,4}, {0,1,3,4}, {1,0,2,4}, {0,1,2,3}},
                                      {{0,2,3,4}, {2,1,3,4}, {1,0,3,4}, {0,1,2,4}, {1,0,2,3}}};


/* void d3initialize(const pd3stateType this, ppointType vertArray, int nvert) 
initialize d3state this from vertArray.  Allocates memory for corners and spheres, 
sets up the free list for tetrahedra, and creates spheres for the first five points.
REQUIRES: vertArray[0] contains the point at infinity and vertArray[1..4] contain 
four finite points; these first five points must be in general position. 
*/
int d3initialize(const pd3stateType this, ppointType vertArray, int nvert) 
{
  int j, p, last;
  double d;

  //  initBitTable(); // BITS
  this->vert = vertArray;
  this->limitmaxTetra = TETperV*nvert; // allocate space for spheres and corners
  this->s   = (pcornerType) calloc(4*this->limitmaxTetra, sizeof(cornerType)); // per corner: v, opp
  this->sph = (psphereType) calloc(this->limitmaxTetra, sizeof(sphereType)); // per tetra: sphere eqn
  this->active = (int *) calloc(this->limitmaxTetra, sizeof(int)); // flag -1 unused, 0 dead, 1 alive

  if ((this->s == NULL) || (this->sph == NULL) || (this->active == NULL)) {
    printf("ERROR: d3batch could not calloc memory for data structures\n");
    exit(EXIT_FAILURE);
  }

  // initialize tetrahedra
  last = -1; /* set up free list of tetrahedra */
  this->freeTetra = this->limitmaxTetra; 
  do {
    this->freeTetra--;
    this->active[this->freeTetra] = 0; // KILL(freeTetra);
    this->s[CORNER(this->freeTetra,0)].opp = last;
    last = this->freeTetra;
  } while (this->freeTetra > 5); 

  this->active[4] = 2; // create first sphere
  makeSphereV(this->sph+4, this->vert+0, this->vert+1, this->vert+2, this->vert+3, this->vert); 
  d = spdot(this->sph+4, this->vert+4, this->vert+3); // if d<0, then we need to swap

  if (d == 0.0) {
		return 0;
//    printf("ERROR: Need first five vertices to be in general position\n"); exit(EXIT_FAILURE);
  }

  if (d < 0) {
    this->sph[4].x = -this->sph[4].x; this->sph[4].y = -this->sph[4].y; this->sph[4].z = -this->sph[4].z; 
    this->sph[4].sq = -this->sph[4].sq; 
  }

  for (p=0; p<5; p++) {
    for (j=0; j<4; j++) {       // pay attention to orientation when assigning vertices
      setCornerVC(CORNER(p,j), this->vert+initialv[d<0][p][j], initialopp[p][j]); 
    }
    makeSphereV(this->sph+p, this->s[CORNER(p,0)].v, this->s[CORNER(p,1)].v, this->s[CORNER(p,2)].v, this->s[CORNER(p,3)].v, this->vert);
    this->active[p] = 1;
    // swap first two if d<0. 
    if (((d<0)&&(p==1)) || ((d>0) && (p==0))) 
      {ASSERT(this->sph[p].sq >0, "Somehow vertp at infinity is in or on sphere p in init");}
    else
      {ASSERT(spdot(this->sph+p, this->vert+p + (d<0)*(p<2)*(1-2*p), this->s[CORNER(p,3)].v) > 0, 
              "Somehow vertp is in or on sphere p in init.");}
  }
  this->liveTetra = 4;
  this->maxTetra = 5; 

	return 1;
}



#define AUDITDROP //printf("Dropping pv %d = (%5d; %5.0f %5.0f %5.0f; %10f)\n", pv-this->vert, pv->index, pv->x, pv->y, pv->z, pv->sq);  cornerPrint4(this, *result)

// LOCATION ROUTINES walk a mesh stored in s, active, & sph starting from corner start to find a point pv.
// result is the index of a simplex whose sphere strictly contains pv.  
// (The simplex itself need not contain pv.)
// The return code is positive if we succeed (# of location steps at present)
// Return code of 0 means that we failed, perhaps because points with small weight have no Voronoi cells
// Return code of -1 means that we found the duplicate of a vertex in the mesh. IN THIS CASE, result 
//   is the location of the corner where we found the duplicate!
int d3locSphere(const pd3stateType this, ppointType pv, int start,
                    int *result) { // output
  int j, guard; // loop variables
  int c1, c2; // corners
  psphereType s1, s2; // spheres 
  double I1, I2, d; // Insphere values
  
  //  printf("LOC:  %d = (%5d; %5.0f %5.0f %5.0f; %10f)\n", pv-this->vert, pv->index, pv->x, pv->y, pv->z, pv->sq);

  guard = 2*this->maxTetra+4; // prevent infinite loops  
  c1 = CORNER(start,0); // corner in start
  s1 = this->sph+start; // sphere at start
  I1 = InSpherev(s1, pv, this->s[c1+3].v); // Check if strictly inside start sphere.
  if (I1 < 0) {// found already
    *result = start; 
    return 1; // success on first try
  }

  while (--guard) { 
    c2 = this->s[c1].opp; s2 = this->sph + TETRA(c2); // nhbr corner, sphere, value
    I2 = InSpherev(s2, pv, this->s[LASTCORNER(c2)].v);
    if (I2 < 0) {  // found one!
      *result = s2 - this->sph;
      return 2*this->maxTetra+5-guard; // number of steps
    } 
    d = s2->sq * I1 - s1->sq * I2; // Warning: if s1 & s2 are same sphere, this is zero
    locateSideCnt++; // STATS
    if (d==0) { // We may be on two spheres---check for duplicate vertex
      if (EQUALPV(pv,this->s[c2].v)) { *result = c2; AUDITDROP; return -1; } 
      j = INDEX(c2);
      if (EQUALPV(pv,this->s[c2+offset[1][j]].v)) { *result = c2+offset[1][j]; AUDITDROP; return -1; }
      if (EQUALPV(pv,this->s[c2+offset[2][j]].v)) { *result = c2+offset[2][j]; AUDITDROP; return -1; }
      if (EQUALPV(pv,this->s[c2+offset[3][j]].v)) { *result = c2+offset[3][j]; AUDITDROP; return -1; }
      // otherwise no duplicate; we probably have s1 == s2. (Rare in protein data.)
      d = RANDBIT-0.5; // choose a random direction, as a hack. (Should do plane computation)
    }
    if (d < 0) // if on I1 side 
      c1 = INCREMENT(c1);
    else {// on I2 side
      c1 = INCREMENT(c2); s1 = s2; I1 = I2; 
    }
  }
  result = 0; // location failure
  return 0; 
}

/* d3compactCorners(const pd3stateType this, pcornerType *result, int *ncorners);
   Takes corner table this->s in which some corners/tetrahedra are unused, and 
   returns a compactified corner table result, and its length ncorners.
   DESTROYS the corner table s and active flags in the process.
   returns false if it is unable to allocate memory.
*/
int d3compactCorners(const pd3stateType this, // arrays this->s & this->active DESTROYED!
                    pcornerType *result, int *ncorners) { // output 
  int c, nc, i, j;
  pcornerType pc; 

  i = 0; // count actives
  for (j = 0; j < this->maxTetra; j++) 
    this->active[j] = (DEAD(j))? -1 : i++; // make old->new pointer dictionary for active tetra

  *ncorners = 4*i; // number of corners to return
  *result = pc = (pcornerType) calloc(*ncorners, sizeof(cornerType)); // return corner table: v, opp
  if (pc != NULL) {
    for (j = 0; j < this->maxTetra; j++) // compact the corners
      if (this->active[j] >= 0) { 
        c = CORNER(j,0);  // old corner c -->  new corner ptr pc = (*result)+0,1,2,...
        pc->v = this->s[c].v; nc = this->s[c].opp; 
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc)); // assign new tetra # w/ old index
        pc++; c++;
        pc->v = this->s[c].v; nc = this->s[c].opp; 
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc));
        pc++; c++;
        pc->v = this->s[c].v; nc = this->s[c].opp; 
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc));
        pc++; c++;
        pc->v = this->s[c].v; nc = this->s[c].opp; 
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc));
        pc++;
      }
  }
  free(this->s);  // free old corner list
  free(this->active);
  return (pc != NULL);  // true if we were successful
}



/* d3insert(const pd3stateType this, int vi, int p) 
   inserts the point this->vert[vi] that is contained in sphere p 
   into the delaunay triangulation stored in this.  
   (d3locSphere may be used to obtain p.)
 */
void d3insert(const pd3stateType this, int vi, int p) {
  int i, j, off;
  int b, c, newb; // corners
  int nc, ni, dead, jdead; // indices
  double d;
  ppointType v0, v1, v2; // pointers to vertices
  ppointType pv = this->vert+vi;

    // Tetrahedra containing pv are "dead", and are pushed onto kill stack.  
    // We use DFS with stack pst to find them and kill them
    // At live-dead boundary, we save dead tetras on stack nhbr,
    //  then make new tetras and hook in to live by setting the last opp pointer.
    //
    // Invariants/operations: Tetrahedron p is marked alive or dead on first visit.
    //    Corner c is pushed on stack when TETRA(this->s[c].opp) is marked dead. 
    //
    // On termination, stack nhbr contains dead corners with live neighbors 
    //    that have new tetras (so this->s[nhbr].opp != this->s[this->s[nhbr].opp].opp temporarily.)
    //    Stack kill contains old tetrahedra for final recycling.

    stkINIT(this->dfs);  // DFS stack holds corners opposite dead tetras
    stkINIT(this->idfs); // iDFS stack holds corners opposite infinite tetras with pv on bdry
                         //    (these are a special case: dead, but don't propagate)
    stkINIT(this->nhbr); // stack for dead corners with live nhbr tetras
    stkINIT(this->kill); // stack of dead tetras to recycle
    b = CORNER(p,0);
    PUSH(p, this->kill); KILL(p); // kill tetra initial p,
    PUSH(this->s[b++].opp, this->dfs); // stack neighbors
    PUSH(this->s[b++].opp, this->dfs);
    PUSH(this->s[b++].opp, this->dfs);
    PUSH(this->s[b  ].opp, this->dfs);

    while (!isEMPTY(this->dfs)) {
      c = POP(this->dfs); p = TETRA(c);
      /*        printf("::Popping %d with opp %d \n", c, this->s[c].opp);/**/
      ASSERT(DEAD(TETRA(this->s[c].opp)), "dfs stack element with non-dead neighbor");
      if (DEAD(p)) continue; // dead already
      d = InSpherev(this->sph + p, pv, this->s[LASTCORNER(c)].v); // Is pv in, out, or on?
      if (d < 0) { // kill and continue dfs if pv is strictly inside 
        KILL(p); PUSH(p, this->kill); // kill and stack tetra
        j = INDEX(c);
        PUSH(this->s[c+offset[1][j]].opp, this->dfs); // stack neighbors to check
        PUSH(this->s[c+offset[2][j]].opp, this->dfs);
        PUSH(this->s[c+offset[3][j]].opp, this->dfs);
      } 
      else if (d > 0  || this->sph[p].sq > 0) {  // pv is outside (or on with sp finite), so
	                                         // c is live neighbor of dead opp tetra this->s[c].opp
        PUSH(this->s[c].opp, this->nhbr); // remember old corner, so we can hook tetra into mesh later
        STARTTETRA(this); // make new tetrahedron liveTetra
        newb = CORNER(this->liveTetra,3); // last corner of new tetra
        setCornerVCN(newb, pv, c); // last corner is pv; also set opposite corner c. Do rest later.
      }
      else { // d==0 && sph[p] is infinite: handle two special cases
          if (this->sph[TETRA(this->s[c].opp)].sq == 0) { // if dead sphere is infinite, too
            PUSH(c, this->idfs); // then if c stays alive, we make tetra to it (flat, but infinite).
          } else { // dead sphere is finite; kill c and make tetras to neighbors, if they stay alive.
            KILL(p); PUSH(p, this->kill); // kill and stack tetra
            j = INDEX(c);
            PUSH(this->s[c+offset[1][j]].opp, this->idfs); // stack neighbors to check
            PUSH(this->s[c+offset[2][j]].opp, this->idfs);
            PUSH(this->s[c+offset[3][j]].opp, this->idfs);
          }
        }
    }
    
    while (!isEMPTY(this->idfs)) { // check the neighbors of infinite tetrahedra
      c = POP(this->idfs); p = TETRA(c);
      /*        printf("::Popping %d with opp %d \n", c, this->s[c].opp);/**/
      ASSERT(DEAD(TETRA(this->s[c].opp)), "dfs stack element with non-dead neighbor");
      if (DEAD(p)) continue; // dead already
      ASSERT(DEAD(TETRA(this->s[c].opp)), "Live corner c should have dead neighbor");
      PUSH(this->s[c].opp, this->nhbr); // remember old corner, so we can hook tetra into mesh later
      STARTTETRA(this); // make new tetrahedron liveTetra
      newb = CORNER(this->liveTetra,3); // last corner of new tetra
      setCornerVCN(newb, pv, c); // last corner is pv; also set opposite corner c. Do rest later.
    }
    
    // Now, we have stack of dead neighbors of live tetras, and we've hooked new tetras to them.  
    while (!isEMPTY(this->nhbr)) {
      dead = POP(this->nhbr); jdead = INDEX(dead); //  dead tetra and index of dropped corner.
      /*        printf("--Popped %d(%d)\n", dead, jdead); /**/
      ASSERT(DEAD(TETRA(dead)), "corner on nhbr stack is not dead!?");
      newb = this->s[this->s[dead].opp].opp-3; // base of new tetra. 

      dead -= jdead; // just use base of dead one.
      // new tetra has 0,1,2,3=pv; 
      // corresponding old indices before jdead is dropped: 
      //   drop[j][0],..,drop[j][3], (no corresp to pv)
      j = jdead; 
      i = drop[jdead][0]; // old index of new corner 0;
      c = dead+i; // note i = INDEX(C);
      v0 = this->s[c].v; // copy vertex v0
      nc = this->s[c].opp; // go to neighbor
      // In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
      // To avoid index calculations, maintain i = INDEX(c), nc = this->s[c].opp, ni = INDEX(nc).
      while (DEAD(TETRA(nc))) { 
        ni = INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
        j = ni; i = ni + off; c = nc + off; nc = this->s[c].opp; // fix new j, i, c, and try neighbor
      } 
      nc = this->s[nc].opp; // go to new tetra
      ASSERT(this->s[nc].v == pv, "Expected to find new tetra using pv after walking dead tetras. ");

      setCornerVC(newb, v0, nc-3+invdrop[i][j]); newb++;

      j = jdead; 
      i = drop[jdead][1]; // old index of new corner 1;
      c = dead+i; // note i = INDEX(C);
      v1 = this->s[c].v; // copy vertex v1
      nc = this->s[c].opp; // go to neighbor
      // In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
      // To avoid index calculations, maintain i = INDEX(c), nc = this->s[c].opp, ni = INDEX(nc).
      while (DEAD(TETRA(nc))) { 
        ni = INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
        j = ni; i = ni + off; c = nc + off; nc = this->s[c].opp; // fix new j, i, c, and try neighbor
      } 
      nc = this->s[nc].opp; // go to new tetra
      ASSERT(this->s[nc].v == pv, "Expected to find new tetra using pv after walking dead tetras. ");

      setCornerVC(newb, v1, nc-3+invdrop[i][j]); newb++;

      j = jdead; 
      i = drop[jdead][2]; // old index of new corner 2;
      c = dead+i; // note i = INDEX(C);
      v2 = this->s[c].v; // copy vertex v2
      nc = this->s[c].opp; // go to neighbor
      // In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
      // To avoid index calculations, maintain i = INDEX(c), nc = this->s[c].opp, ni = INDEX(nc).
      while (DEAD(TETRA(nc))) { 
        ni = INDEX(nc); off = indoff[i][ni][j]; // where j goes relative to i is our new i.
        j = ni; i = ni + off; c = nc + off; nc = this->s[c].opp; // fix new j, i, c, and try neighbor
      } 
      nc = this->s[nc].opp; // go to new tetra
      ASSERT(this->s[nc].v == pv, "Expected to find new tetra using pv after walking dead tetras. ");

      setCornerVC(newb, v2, nc-3+invdrop[i][j]); newb++;

      c = this->s[this->s[dead+jdead].opp].opp;
      ASSERT(v0==this->s[CORNERINOPP(0,c)].v, "v0 does not line up");
      ASSERT(v1==this->s[CORNERINOPP(1,c)].v, "v1 does not line up");
      ASSERT(v2==this->s[CORNERINOPP(2,c)].v, "v2 does not line up");

      makeSphereV(this->sph+TETRA(newb), v0, v1, v2, pv, this->vert); // use either this or makeSphereP above
    }
}



/* d3.c Delaunay/Power diagram function 
   d3batch takes lifted input vertices, and returns a (compact) corner table for Delaunay.
   REQUIRES that the first point is at infinity, 
   and that the first 5 are in general position. (I should verify or relax this.)
   Since all points have radii assigned already, it can compute power diagrams.
*/
void d3batch(ppointType vertArray, int nvert, // input vertices (pointType[]) & number of vertices
        pcornerType *result, int *ncorners)  // output corner table
{
  int k;
  int p; // tetra
  int vi; // vertex index in outermost loop
  ppointType pv; // vertex pointer in outermost loop
  pcornerType pc; // corner for copying to result
  d3stateType d3; // state data structure
  const pd3stateType this = &d3;

	int itry = 5;
	pointType vtemp;
	int j=0;
	while(!d3initialize(this, vertArray, nvert)) // JSS: does this leak memory?
	{
		 printf("Permute points to get the first five vertices to be in general position\n");  
		 j = (j+1)%5; 
		 vtemp = *(vertArray+j);
		 *(vertArray+j) = *(vertArray+itry);
		 *(vertArray+itry) = vtemp;
		 itry++;
		 
		 if (itry==nvert)
		 {
				printf("Can't find vertices to initialize. Quit\n"); 
				exit(EXIT_FAILURE);
		 }
	}
#ifdef STATS
  this->killmax = this->dfsmax = this->idfsmax = this->nhbrmax = -1; // init STATS
#endif
  maxLocate = locateSideCnt = sphereCnt = startTetraCnt = freeTetraCnt= inSphereCnt=0; 

  for (vi = 5; vi < nvert; vi++) { // incrementally insert vert[vi]
    //LOCATE: find some tetrahedron with sphere strictly containing vert[vi]
    if ((k = d3locSphere(this, this->vert+vi, this->liveTetra, &p)) < 1) 
      continue; // if we fail to locate (duplication or other reason) just skip vert[vi]
    if (k > maxLocate) maxLocate = k; // STATS
#ifndef NOASSERT
    if (vi%2000 == 0) { // print a few search results
      printf("%%Found %3d in %4d(%3d) after %d steps %f\n", 
             vi, CORNER(p,0), p, k, InSpherev(this->sph+p, this->vert+vi, this->s[CORNER(p,3)].v));
      auditCorners(this, 1); // audit, and check spheres, too
      (void)fflush(stdout); 
    } /**/
#endif
    d3insert(this, vi, p); //  insert vertex vi, which is in sphere of tetra p.
    while (!isEMPTY(this->kill)) { // recycle memory of dead tetrahedra
      p = POP(this->kill);
      FREETETRA(this, p);
    }
  }
#ifdef STATS
  printf(" %10d\t Points\n %10d\t Max Tetra\n", vi, this->maxTetra);
      printf("We performed:\n %10d\tinSphere tests\n %10d\tplane tests\n %10d\ttetrahedra created-\n %10d\tfreed =\n %10d\ttetrahedra\n %10d\tsphere equations computed.\n",  
             inSphereCnt, locateSideCnt, startTetraCnt, freeTetraCnt, startTetraCnt-freeTetraCnt, sphereCnt);
      printf("Max locate steps=%d, max killed=%d, max inf=%d, max created=%d, maxdfs=%d\n", 
             maxLocate, this->killmax, this->idfsmax, this->nhbrmax, this->dfsmax);
      fflush(stdout);
      auditCornersAux(this,1); /*TESTING*/
#endif
  auditCornersAux(this,1); /*TESTING*/

  free(this->sph);  // done with spheres; free them

  if (!d3compactCorners(this, result, ncorners)) { // disposes s and active!!
    printf("ERROR: d3.c could not allocate corner table to return\n");
    exit(EXIT_FAILURE);
  }
}
