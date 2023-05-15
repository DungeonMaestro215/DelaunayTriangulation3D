# README.txt for d3 delaunay computation for pdb   Jack Snoeyink Nov 03

Implements Watson's incremental Delaunay (and power diagram), 
with sphere-based search and simplical data structure storing 
vertices and opposite neighbors.

Points must be scaled to integers; guaranteed for differences of 10
bits, and uses leveling and hilbert curve to guarantee 16 bits if the
points are well distributed.  Works for pdb files, which are 20 bits.

Handles degeneracies by perturbing points by increasing infinitesimals
to guarantee that all simplices are full-dimensional.

Files: 
delaunay3.h	describes basic data structures
d3.c		incremental Delaunay: d3batch takes in a list of points and returns a corner table.
	assumptions: the 0th vertex is at infinity and next 4 are in general position.
		Also that the points are well distributed and inserted in an order so that spheres 
		can be constructed exactly in double precision.  
		perturbs points by order; guarantees all simplices are full dimensional. 
d3permute.c	permutes points by Hilbert curve to prepare them for d3batch
	assumption: points are provided by a "reader".
pdbreader.c	readers for pdb and txt files; also random generation of points
d3boxOrder.c	Lookup tables that d3permute uses to do sort along Hilbert curve. Produced by hilbOrder.
hilbOrder.c	Makes Hilbert curve orders from recursive rules.
test.c		example program that reads pdb files (named as command-line arguments)


d3batch takes in a list of points and returns a corner table:
A corner is a vertex-use in a tetrahedron.  
 Each has a pointer to its vertex and 
 to the corner in the neighboring tetra across the face opposite the vertex.
Tetrahedra are groups of four consecutive corners,  ordered by increasing 
vertex index, except that the first two may be swapped to ensure 
that the orientation determinant is positive.  
I.e., take the lex smallest alternating permutation with positive sign. 
There must always be an odd number of swaps between two permutations; 
we swap the first two if necessary to achieve this. 

Compiling without -DNOASSERT heavily audits data structures as it runs.
Compiling with -DSTATS keeps some stats.