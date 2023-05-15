#include <stdio.h>
#include <stdlib.h>

#define CORNERINOPP(i,c) (this->s[c].opp + indoff[INDEX(c)][INDEX(this->s[c].opp)][i])
static const short indoff[4][4][4] 
= {/*        0ABCD       1ABCD         2ACBD          3ABCD     */
  /* 0: */ {{5,5,5,5},  { 0,-1, 1, 2},  { 0,-1,-2, 1},  { 0,-3,-2,-1}},
  /* 1: */ {{1,0,2,3},  { 5, 5, 5, 5},  {-2, 0,-1, 1},  {-2, 0,-3,-1}}, 
  /* 2: */ {{2,1,0,3},  {-1, 1, 0, 2},  {-1,-2, 0, 1},  {-3,-2, 0,-1}}, 
  /* 3: */ {{1,2,3,0},  { 1,-1, 2, 0},  {-2,-1, 1, 0},  {-2,-3,-1, 0}}}; 

int main(int argc, char** argv) {
    int i, j, k;
    scanf("%d %d %d", &i, &j, &k);

    printf("%d\n", CORNERINOPP(0, 3));

    return EXIT_SUCCESS;
}