#include "pdbreader.h"

int hetatm;
double mult;

FILE *pdbopener(char *fname, int argc, char **argv) {
  int i;
  FILE *fid = fopen(fname, "r");
  if (fid == NULL) {
    printf("Could not open input file %s.\n", fname);
    exit(EXIT_FAILURE);
  }
  hetatm = FALSE;   // look for -het option
  mult = 100.0; // default multiplier for pdb
  for (i = 1; i < argc; i++) {
    if ((0 == strncasecmp(argv[i],"-het", 4))) 
      hetatm = TRUE;
    if ((0 == strncasecmp(argv[i],"-m", 2))) 
      if (1 != sscanf(argv[i], "-m%f", &mult))
	printf("couldn't sscanf multiplier -m#\n");
  }
  return fid;
}

int AtomType(char* aname) {
  int t = 0;
  switch(aname[2]) {
  case 'C': t = 1; break;
  case 'N': t = 2; break;
  case 'O': t = 3; break;
  case 'H': t = 4; break;
  case 'S': t = 5; break;
  }
  return t;
}

int pdbreader(FILE *fid, ppointType v) {
  char line[132];
  double rad;
  int nvert = 0;
  static const double AtomRadius[] = {1.0, 0.8, 1.0, 0.6, 1.2, 1.0};
  while (fgets(line, 132, fid) != NULL) 
    if ((0 == strncasecmp(line,"ATOM  ", 6)) 
	|| (hetatm && (0 == strncasecmp(line,"HETATM", 6)))) {
      rad = AtomRadius[AtomType(line+NAMEOFFSET)]; // assign radius based on atom type
      setVert(v+nvert, nvert, atof(line+XOFFSET), atof(line+YOFFSET), atof(line+ZOFFSET), rad, mult);
      nvert++;
    }
  fclose(fid);
  return nvert;
}

int totcount; // count points read for index
FILE *txtopener(char *fname, int argc, char **argv) {
  // open input file
  int i;
  FILE *fid = fopen(fname, "r");

  if (fid == NULL) {
    printf("Could not open input file %s.\n", fname);
    exit(EXIT_FAILURE);
  }
  totcount = 0;

  mult = 1.0; // default multiplier for txt
  for (i = 1; i < argc; i++) {
    if ((0 == strncasecmp(argv[i],"-m", 2))) 
      if (1 != sscanf(argv[i], "-m%lf", &mult))
	printf("couldn't sscanf multiplier -m#\n");
      else
	printf("Multiplier %lf\n", mult);
  }
  return fid;
}

int txt3reader(FILE *fid, ppointType v) {
  coord rx, ry, rz; 
  int nvert = 0;
  while (fscanf(fid, "%f %f %f", &rx, &ry, &rz)==3) {
    setVert(v+nvert, nvert, rx, ry, rz, 0, mult);
    nvert++;
  }
  fclose(fid);
  return nvert;
}

int txt4reader(FILE *fid, ppointType v) {
  coord rx, ry, rz; 
  double rad;
  int nvert = 0;
  while (fscanf(fid, "%f %f %f %lf", &rx, &ry, &rz, &rad)==4) {
    setVert(v+nvert, nvert, rx, ry, rz, rad, mult);
    nvert++;
  }
  fclose(fid);
  return nvert;
}


int nrndpoints;
FILE *rndopener(char *fname, int argc, char **argv) {
  nrndpoints = NVERT;
  if (1 != sscanf(fname, "-%d", &nrndpoints))
    printf("couldn't sscanf\n");
  printf("Making %d points\n", nrndpoints);
  fflush(stdout);
  mult = 1.0; // no need to multiply random
  return 0;
}

int rndreader(FILE *dummy, ppointType v) {
  int i; 
  
  printf("Making %d points\n", nrndpoints);
  for (i = 0; i < nrndpoints; i++)
    setVert(v+i, i, RANDCOORD, RANDCOORD, RANDCOORD, 0.0, mult);
  return nrndpoints;
}

FILE *hieropener(char *fname, int argc, char **argv) {
  char line[132];
  FILE *fid = txtopener(fname, argc, argv);
  if (fid != NULL) {
    fgets(line, 132, fid);
    printf("%s\n", line);
  }
  return fid;
}

int hierreader(FILE *fid, ppointType v) {
  int nvert = fread(v, sizeof(pointType), MAXVERT, fid); // read all vertices
  fclose(fid);
  return nvert;
}











