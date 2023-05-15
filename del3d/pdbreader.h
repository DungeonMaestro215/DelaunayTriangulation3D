/* pdb reader for makehier. also makes random points and reads plain text. */
#include "delaunay3.h"

//PDB file offsets
//ATOM atomno  aname alt rname chain resno   x%8.3f y%8.3f z%8.3f
//1-6  7-11 12 13-16 17 18-20  22    23-26    31-38  39-46  47-54
#define SERIALOFFSET 7
#define NAMEOFFSET 13
#define XOFFSET 31
#define YOFFSET 39
#define ZOFFSET 47

FILE *pdbopener(char *fname, int argc, char **argv);
int pdbreader(FILE *fid, ppointType v);

FILE *rndopener(char *fname, int argc, char **argv);
int rndreader(FILE *fid, ppointType v);

FILE *txtopener(char *fname, int argc, char **argv); 
int txt3reader(FILE *fid, ppointType v); // x y z
int txt4reader(FILE *fid, ppointType v); // x y z rad

FILE *hieropener(char *fname, int argc, char **argv);
int hierreader(FILE *fid, ppointType v);



