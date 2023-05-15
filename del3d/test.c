/* test.c                   Jack Snoeyink Oct 2003
test program for d3.c, d3permute.c 
*/

#include <time.h>
#include "delaunay3.h"
#include "pdbreader.h"

int main(int argc, char** argv) {
  int i, j, ii, jj;
  int hashtb[HIGHCOORD];
  int nexttb[NVERT];
  ppointType vert,v;
  int nvert, nv;
  pcornerType result;
  int nresult;

  clock_t tic, toc;

  FILE *fid;

  /*  fid = rndopener("none", argc, argv);
  d3permute(fid, rndreader, &vert, &nvert);
  fclose(fid);/**/

  /* fid = txtopener("hilb6.txt", argc, argv);
  d3permute(fid, txt3reader, &vert, &nvert);
  fclose(fid); /**/
  
  ii = 1;
  for (ii = 1; ii < argc; ii++) 
    if (argv[ii][0] != '-') {
      fid = pdbopener(argv[ii], argc<4 ? argc : 4, argv);
      d3permute(fid, pdbreader, &v, &nv);
      fclose(fid);
      tic = clock(); // start timer
      jj = 1; 
#ifdef NOASSERT
      for (jj = 0; jj < 10; jj++) 
#endif
	{
	  d3batch(v, nv, &result, &nresult);
	  free(result); 
	}
      toc = clock();
      printf("%s %d %d Time(secs) %f\n\n", argv[ii], nv, nresult, 
	     ((double) (toc - tic)) / CLOCKS_PER_SEC / jj); (void)fflush(stdout); /**/
    }
  return EXIT_SUCCESS;
}
