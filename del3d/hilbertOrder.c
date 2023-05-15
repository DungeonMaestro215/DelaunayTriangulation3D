/* hilbertOrder.c                        Jack Snoeyink July 2003
This program makes lookup tables for computing hilbert orders for 
boxes points in dimensions up to 6.  See usage string.
Some examples:  
planar hilbert curve with 3 levels of recursion: 
ho  xy 2 yYxXxyy 3
Expanded rule: yYxXxyy
 rule^2: xXyYyxxYyYxXxyyXyYxXxyyyXxYyYXX
 rule^3: yYxXxyyXxXyYyxxYxXyYyxxxYyXxXYYYxXyYyxxYyYxXxyyXyYxXxyyyXxYyYXXXxXyYyxxYyYxXxyyXyYxXxyyyXxYyYXXyYyXxXYYxXxYyYXXyXxYyYXXXyYxXxyy

3d snake in 8^3 boxes:
$ ./ho xyz 8 [[{[ZZ]7}/YYY~/Y]3.=/YYX~/X.=/Xyy~/.[y=/yyy~/]3X]4- 1 snake8
Expanded rule: ZZZZZZZZZZZZZZYYYzzzzzzzzzzzzzzYZZZZZZZZZZZZZZYYYzzzzzzzzzzzzzzYZ
...
zzzzzzzzzzzzzzyZZZZZZZZZZZZZZyyyzzzzzzzzzzzzzzyZZZZZZZZZZZZZZyyyzzzzzzzzzzzzzz

3D hilbert curve with three levels of recursion:
$ ./ho xyz 2 {ZZYYXzz}/.X.~/ 3 hilbxyz3.txt
Expanded rule: ZZYYXzzXZZXyyzz
 rule^2: XXYYZxxZXXZyyxxZZZXXYzzYZZYxxzzYZZYYXzzXZZXyyzzzxxYYzXXzxxzyyXXXXXYYZxx
ZXXZyyxxZZZYYXzzXZZXyyzzyZZxxyzzyZZyXXzzzxxYYzXXzxxzyyXX
 rule^3: ZZYYXzzXZZXyyzzXZZXXYzzYZZYxxzzYXXYYZxxZXXZyyxxxzzYYxZZxzzxyyZZZZZYYXzz
...
XZZXyyzzzzzYYxZZxzzxyyZZxxxYYzXXzxxzyyXXyZZxxyzzyZZyXXzzXZZYYXzzXZZXyyzz
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCOORDS 6 // maximum number of coordinates
#define TWICEMAXCOORDS 12
#define TMCPLUS1 13
#define LIMIT 300000

#define DOT(p,q) (p[0]*q[0] + p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4] + p[5]*q[5])


#define  checkDuplicates(s, n, msg) { \
  for (i = 0; i < n-1; i++)\
    for (j = i+1; j < n; j++)\
      if (s[i] == s[j]) {printf("%s has duplicates at %d and %d\n", msg, i,j); exit(EXIT_FAILURE); }\
}

/* expandRule processes the special character in a rule to leave only coordinate letters */
static char *expandUsage = 
"     The following special characters can be used when specifying rules:\n"
"       [...]digit : Repeat the enclosed the number of times specified by digit.\n"   
"       {   }char  : Remember the enclosed associated with digit.\n"   
"       =char      : Copy the remembered expression.\n"
"       ~char      : Reverse & complement as you copy remembered expression.\n"
"       . (period) : Ignored (for spacing)\n"
"       - (minus)  : back up one character\n"
"       ! (bang)   : start over (but remember saved strings)\n"
"     These are processed immediately, leaving only coordinate letters,\n"
"       so use ---.--- to repeat minus, not [-]6.\n";

char *expandRule(char * s) {
  int i,j,k,m, len, first;
  int stack[20], sp; // stack for parens, brackets: indices into rule
  int start[256], end[256]; // for remembering strings: indices into rule
  char *result;
  char *rule = (char *)calloc(LIMIT, sizeof(char *)); 

  len = strlen(s); 
  sp = 0; // uses stack for [] {}
  bzero((void *) start, sizeof(start)); // for saving with {}char
  bzero((void *) end, sizeof(end)); 

  first = j = 0; // copy s[i] to rule[j], expanding repetitions and macros
  for (i = 0; i < len; i++) {
    switch (s[i]) {
    case '.': break;
    case '!': first = j; break; // start over
    case '-': j--; rule[j] = '\0'; break; // back up one
    case '[':
    case '{': stack[sp++] = j; break;
    case '}': k = s[++i]; // get following ASCII char
      if (sp <= 0) { printf("} with stack pointer %d <= 0\n", sp); exit(EXIT_FAILURE);}
      end[k] = j;
      start[k] = stack[--sp];
      break;
    case ']': k = s[++i]-'0'; // get following digit
      if (k > 9) { printf("] not followed by digit at %d\n", i); exit(EXIT_FAILURE);}
      if (sp <= 0) { printf("] with stack pointer %d <= 0\n", sp); exit(EXIT_FAILURE);}
      m = j - stack[--sp]; // how many to copy
      while (k > 1) { 
	rule[j] = '\0';
	strncpy(rule+j, rule+stack[sp], m); // make copies
	j += m;
	k--;
      }
      break;
    case '=': k = s[++i]; // get following ASCII char
      if (start[k] >= end[k]) {printf(" Nothing saved in %c at %d\n", (char)k, i); exit(EXIT_FAILURE);}
      rule[j] = '\0';
      strncpy(rule+j, rule+start[k],  end[k]-start[k]); // copy
      j += end[k]-start[k];
      break;
    case '~': k = s[++i]; // get following ASCII char
      if (start[k] >= end[k]) {printf(" Nothing saved in %c at %d\n", (char)k, i); exit(EXIT_FAILURE);}
      for (m = end[k]-1; m >= start[k]; m--) 
	rule[j++] = rule[m] + ((rule[m] >= 'a')? 'A'-'a': 'a'-'A'); // copy in reverse & change case
      break;
    default: 
      rule[j++] = s[i];
    }
  }
  rule[j] = '\0';

   result = (char *)calloc(j-first+1, sizeof(char));
   strcpy(result, rule+first);
   free(rule);
   return result; 
}

/* checkRule takes in the rule, crange, coordinates (both lower & upper), and
ncoords (counting the lower only).  It parses the rule and, if successful, 
returns the total number of boxes, order, and hcase.  
Otherwise it exits the program with an error message.
*/
char *ruleUsage = 
"     In recursiveRule, the even coordinate letters are recursive cases;\n"
"       odd letters increment (upper) or decrement (lower) the coordinate\n"
"     E.g., hilbertOrder xy 1 yYxXxyy 3 uses 1 bit for each coordinate, xy,\n"
"     and the rule incrementing Y, X, then decrementing y is applied 3 times,\n"
"     making 4^3 entries per case and 8 cases.\n\n"
"     The recursive rule operates on ncoord integer and bit arrays,\n"
"     using two alternating operations:\n"
"       odd incr(upper) allowed only if bit=1 / decr(lower) only if bit=0;\n"
"       both even & odd letters toggle coordinate in bit array.\n"
"     The rule must go from integers 0000 to #000, where #=cRange-1,\n"
"       and from bits 0000 to 1000.  Thus, rule length is 2*cRange^ncoords-1.\n";

void checkRule(char *rule, short *inbits, short recurNo, int crange, char *coords, int ncoords, 
	       short *order, short *hcase) {
  int i,j, n;
  int nrule;
  short sign; // for increment or decrement
  short bits[MAXCOORDS];
  short ints[MAXCOORDS];  // stores the current integers
  int bmult[MAXCOORDS], hcasebit; // multipliers so bits.bmult = hcase index 
  int omult[MAXCOORDS]; // multipliers so ints.bmult = order of row major box storage
  short k;
  int nboxes; // number of boxes generated, at most LIMIT
  //  short order[LIMIT]; // to record the order that each pattern is generated 
  //  short hcase[LIMIT]; // to record the hilbert case for each pattern
  char *pc; 

  nrule = strlen(rule);
  nboxes = (int)(pow((double)crange,(double)ncoords)+0.5); // #boxes = (crange)^ncoords
  hcasebit = 1 << ncoords; // which bits are the hcase?
 
  if (nrule != 2*nboxes -1) {
    printf("Expanded rule has %d chars; expecting %d\n", nrule, 2*nboxes-1); exit(EXIT_FAILURE); 
  }

  for (i = 0; i<MAXCOORDS; i++) {
    bits[i] = inbits[i]; 
    ints[i] = inbits[i]*(crange-1); 
    bmult[i] = omult[i] = 0;
  }
  bmult[ncoords-1] = omult[ncoords-1] = 1; // here we set up multipliers so dot(ints, bmult) gives
  for (i = ncoords-1; i > 0; i--) {    // order/hcase array index consistent with row-major order
    bmult[i-1] = bmult[i]*2;
    omult[i-1] = omult[i]*crange;
  }

  for (i = 0; i < nboxes; i++) // set to -1 to indicate nothing there yet.
    order[i] = hcase[i] = -1;

  n = 0;
  for (i = 0; i < nrule; i++) {
    pc = strchr(coords, (int)(rule[i])); 
    if (pc == NULL) {
    printf("Expanded rule contains non-coordinate: %c at %d\n", rule[i], i); exit(EXIT_FAILURE); 
    }
    k = (pc-coords); // k = bit number for coordinate
    if (k >= ncoords) { k -= ncoords; sign = 1; }
    else sign = -1;

    if (0 == (i&1)) { // even cases: record hcase
      order[n] = DOT(ints,omult);  // dot product gives index
      hcase[n] = k*hcasebit + DOT(bits,bmult);
      n++;
      bits[k] ^= 1; // toggle kth bit
    } else {
      if ( (sign>0)!=(bits[k]>0) ) {
	printf("Expanded rule can't %c at %d because sign=%d and bits[%d]=%d\n", 
	       rule[i], i, sign, k, bits[k]); exit(EXIT_FAILURE); 
      }
      ints[k] += sign;
      bits[k] ^= 1;
      if ((ints[k] < 0)||(ints[k]>=crange)) {
	printf("Expanded rule can't %c at %d because makes ints[%d] =%d\n", 
	       rule[i], i, k, ints[k]); exit(EXIT_FAILURE); 
      }
    }
  }

  i = (bits[recurNo] != (1-inbits[recurNo])) || (ints[recurNo] != (1-inbits[recurNo])*(crange-1));
  for (j = 0; j<ncoords; j++)
    if (j != recurNo)
      i = i || (bits[j] != inbits[j]) || (ints[j] != inbits[j]*(crange-1));
  if (i) {
    printf("unexpected final from %d\ni: inbits  bits  ints\n", recurNo);
    for (j = 0; j<ncoords; j++) 
      printf("%1d:%5d %5d %5d\n", j, inbits[j], bits[j], ints[j]);
  }
  checkDuplicates(order, nboxes, "Duplicate box:  boxOrder");
}

void writeOutput(FILE *fid, char *varname, int id, short *array, int n) {
  int i;

  fprintf(fid, "// %s[%d] has %d\n", varname, id, n);
  fprintf(fid, " {%3d", array[0]);
   for (i = 1; i < n; i++) {
    if (i%16 == 0) fprintf(fid, "\n ");
    fprintf(fid, ",%3d", array[i]);
  }
  fprintf(fid, "}\n");
}


#define swap(a,b,tmp) { tmp = a; a = b; b = tmp; }

void ruleTransform(short *bits, short recurNo, char *coords, int ncoords, 
		   char *rule, char *xformed) {
  int i, j, nrule = strlen(rule); 
  char tmp, LUT[256];

  for (i = 0; i < ncoords; i++) { // set up LUT to swap case one bits
    LUT[coords[i]] = coords[i + bits[i]*ncoords];
    LUT[coords[i+ncoords]] = coords[i + ncoords - bits[i]*ncoords];
  }
  swap(LUT[coords[0]], LUT[coords[recurNo]], tmp); // swap based on 
  swap(LUT[coords[ncoords]], LUT[coords[recurNo+ncoords]], tmp); // swap based on 

  j = 0; // copy to transformed
  for (i = 0; i < nrule; i++) 
    xformed[j++] = LUT[rule[i]];
  xformed[j] = '\0';
}

char *ruleApply(char *rule, char *coords, int ncoords, 
	       char *oldstr) {
  
  int i,j, n;
  int nrule, nold;
  short sign; // for increment or decrement
  short ints[MAXCOORDS];  // stores the current integers
  short bits[MAXCOORDS];
  short k;
  char *pc, *newstr;

  for (i = 0; i<MAXCOORDS; i++)
    bits[i] = ints[i] = 0;

  nrule = strlen(rule);
  nold = strlen(oldstr);
  newstr = (char *)calloc((nold/2+1)*nrule+1, sizeof(char *)); 
  n = 0; // copy from oldstr[i] to newstr[n]
  for (i = 0; i < nold; i++) {
    pc = strchr(coords, (int)(oldstr[i])); 
    if (pc == NULL) {
    printf("Expanded oldstr contains non-coordinate: %c at %d\n", oldstr[i], i); exit(EXIT_FAILURE); 
    }
    k = (pc-coords); // k = bit number for coordinate
    if (k >= ncoords) { k -= ncoords; sign = 1; }
    else sign = -1;

    if (0 == (i&1)) { // even cases: process oldstr and expand rule into newstr
      ruleTransform(bits, k, coords, ncoords, rule, newstr+n);
      n += nrule;
      bits[k] ^= 1; // toggle kth bit according to oldstr
    } else { // odd cases: process oldstr and copy to newstr
      newstr[n++] = oldstr[i];
      ints[k] += sign;
      bits[k] ^= 1;
    }
  }
  newstr[n] = '\0';
  return newstr;
}

void pointOutput(FILE *outfile, char *coords, int ncoords, char *oldstr) {
  int i,j, n;
  int nold;
  short sign; // for increment or decrement
  short ints[MAXCOORDS];  // stores the current integers
  short bits[MAXCOORDS];
  short k;
  char *pc, *newstr;

  for (i = 0; i<MAXCOORDS; i++)
    bits[i] = ints[i] = 0;

  nold = strlen(oldstr);
  n = 0; 
  for (i = 0; i < nold; i++) {
    pc = strchr(coords, (int)(oldstr[i])); 
    if (pc == NULL) {
    printf("pointOutput contains non-coordinate: %c at %d\n", oldstr[i], i); exit(EXIT_FAILURE); 
    }
    k = (pc-coords); // k = bit number for coordinate
    if (k >= ncoords) { k -= ncoords; sign = 1; }
    else sign = -1;
    if (0 == (i&1)) { // even cases: process oldstr and expand rule into newstr
      fprintf(outfile,"%d %d %d\n", (ints[0]<<6), (ints[1]<<6), (ints[2]<<6));
      bits[k] ^= 1; // toggle kth bit according to oldstr
    } else { // odd cases: 
      ints[k] += sign;
      bits[k] ^= 1;
    }
  }
}


// This function is a quick hack to generate the specific tables 
// that I need for boxOrder.c in delaunay computations.
//  That is, I need order and hcase to contain the boxes in order and their hcases.
//  The hcases should be the first index into the order table such that 
//  0=000X, 1=100X, and the rest are consecutive integers in an arbitrary mapping, 
//   stored by HC2ind and ind2HC.
//  In other words, boxOrder[j][i] contains the ith box for the case ind2HC[j], 
//  and hcase[j][i] would be the corresponding hcase for recursion.  
//  I use hcases only for j=0 and j=1, so I don't output any others, 
//   which means that I need boxOrders only for a few js, saving more memory.
//
#define HCLIMIT (1<<MAXCOORDS)*MAXCOORDS
void doOutput(FILE *outfile, char *result, int crange, char *coords, int ncoords) {
  int i,j, n, hcaselimit;
  short k;
  int totalBoxes; // total number of boxes generated recursively. at most LIMIT
  short order[LIMIT]; // to record the order that each pattern is generated 
  short hcase[LIMIT]; // to record the hilbert case for each pattern
  short hcase1[LIMIT]; // to record the hilbert case for each pattern
  short myhcase[LIMIT]; // to record reindexed hilbert cases
  short HC2ind[HCLIMIT], ind2HC[HCLIMIT]; // which hcases are used
  short bits[MAXCOORDS];
  char *tmpstr = (char *)calloc(strlen(result)+1, sizeof(char)); // working storage

  totalBoxes = (int)pow((double)crange, (double)ncoords); // grand total # of boxes 
  if (totalBoxes > LIMIT) {
    printf("Too many total boxes: %d > LIMIT; increase LIMIT", totalBoxes); exit(EXIT_FAILURE); 
  } 
  hcaselimit = ncoords*(1 << ncoords); // above highest hcase

  for (i = 0; i < totalBoxes; i++) order[i] = hcase[i] = -1;
  for (i = 0; i < hcaselimit; i++) HC2ind[i] = ind2HC[i] = -1;
  
  // record hcases needed for the 0000X rule
  for (i = 0; i < ncoords; i++) bits[i] = 0;
  checkRule(result, bits, 0, crange, coords, ncoords, order, hcase); // get hcase0
  for (i = 0; i < totalBoxes; i++)
    HC2ind[hcase[i]]++; // determine which Hcases we need

  // record hcases needed for the 1000X rule
  bits[0] = 1;
  ruleTransform(bits, 0, coords, ncoords, result, tmpstr);
  checkRule(tmpstr, bits, 0, crange, coords, ncoords, order, hcase1); // get hcase1
  for (i = 0; i < totalBoxes; i++)
    HC2ind[hcase1[i]]++; // determine which Hcases we need
  HC2ind[0] = 0;            // cases 0=0000X 
  HC2ind[1<<(ncoords-1)] = 1; // and 1=1000X are special

  n = 2;
// reindex the Hcases by making HC2ind
  for (i = 1; i < hcaselimit; i++)  
    if (i != 1<<(ncoords-1)) // skip the 1000X case; already handled.
      if (HC2ind[i] >= 0) {
	/*	printf("%d->%d ", i, n);/**/
	HC2ind[i] = n++; // set the index for this hcase
      } 
  /*  printf("HC2ind\n");/**/

  for (i = 0; i < hcaselimit; i++) 
    if ((HC2ind[i] >= 0)) // construct inverse index ind2HC
      ind2HC[HC2ind[i]] = i; 

  /*  for (j = 0; j < n; j++) printf("%d ", ind2HC[j]); printf("ind2HC\n");/**/

  for (i = 0; i < totalBoxes; i++) myhcase[i] = HC2ind[hcase[i]];
  fprintf(outfile, "static const short hilbertCase[2][%d] = {\n", totalBoxes);
  //  writeOutput(outfile, "HCASE 0000X", ind2HC[0], hcase, totalBoxes);
  writeOutput(outfile, "hcase 0000X", ind2HC[0], myhcase, totalBoxes);
  for (i = 0; i < totalBoxes; i++) myhcase[i] = HC2ind[hcase1[i]];
  fprintf(outfile, ",", totalBoxes);
  writeOutput(outfile, "hcase 1000X", ind2HC[1], myhcase, totalBoxes);
  fprintf(outfile, "};\n");

  fprintf(outfile, "static const short boxOrder[%d][%d] = {\n", n, totalBoxes);
  for (j = 0; j < n; j++) { // write boxOrder for all Hcase indices
    for (i = 0; i < ncoords; i++) // recover bits & k from case index
      bits[ncoords-1-i] = (ind2HC[j]>>i)&1; 
    k = ind2HC[j]>>ncoords;
    if (i>0) fprintf(outfile, ",");
    ruleTransform(bits, k, coords, ncoords, result, tmpstr);
    checkRule(tmpstr, bits, k, crange, coords, ncoords, order, hcase); // get order
    writeOutput(outfile, "boxOrder", ind2HC[j], order, totalBoxes);  
  }
  fprintf(outfile, "};\n");
}

int main(int argc, char ** argv) {
  int i,j, n; 
  int ncoords, nrule, depth;
  char coords[TMCPLUS1];
  int crange; // coordinate range: [0..crange)
  char *rule, *tmpstr, *result; 
  int totalBoxes; // total number of boxes generated recursively. at most LIMIT
  short bits[MAXCOORDS] = {0,0,0, 0,0,0};
  short order[LIMIT]; // to record the order that each pattern is generated 
  short hcase[LIMIT]; // to record the hilbert case for each pattern
  FILE *outfile;

  if ((argc < 5)||(argc > 7)) {
    printf("hilbertOrder coordinates coordRange recursiveRule recursionDepth {outfile} {pointfile}\n"
	   "     produces boxOrder and hilbertCase tables for boxOrder.c\n\n"
	   "     The coordinates must be distinct letters, disregarding case.\n"
	   "     The box coordinate values come from [0..coordRange).\n"
	   "%s\n%s\n"
	   "     Output is written to outfile and a point in each box to pointfile.\n",
	   expandUsage, ruleUsage);
    exit(EXIT_SUCCESS); 
  }

  // check arguments: 1 == coords
  strcpy(coords, argv[1]);
  ncoords = strlen(argv[1]);
  checkDuplicates(coords, ncoords, "Coords");

  for (i = 0; i < ncoords; i++)
    coords[ncoords+i] = coords[i] + ('A'-'a'); // duplicate second part in upper case
  coords[2*ncoords] = '\0';

  if (strcasecmp(coords+ncoords, argv[1]) != 0) {
    printf("Coordinates %s must be all lower case\n", argv[1]); exit(EXIT_FAILURE); 
  }

  if ((1 != sscanf(argv[2], "%d", &crange)) || (crange < 1) || (crange > 32)) {
    printf("Bad coordinate range limit: 2 <= crange <= 32, not %d\n", crange); exit(EXIT_FAILURE); 
  }

  if ((1 != sscanf(argv[4], "%d", &depth)) || (depth < 0) || (depth > 6)) {
    printf("Couldn't get 1 <= depth <= 6\n"); exit(EXIT_FAILURE); 
  }

  totalBoxes = (int)pow((double)crange, (double)ncoords*depth); // grand total # of boxes 
  if (totalBoxes > LIMIT) {
    printf("Too many total boxes: %d > LIMIT; increase LIMIT", totalBoxes); exit(EXIT_FAILURE); 
  } 

  // Here is where the real work begins.
  // We expand special characters in the rule, then parse the rule to check it,
  // then apply it recursively and dump the output
  rule = expandRule(argv[3]);   // expand argv[3] into rule, processing special characters.  
  printf("Expanded rule: %s\n", rule);
  checkRule(rule, bits, 0, crange, coords, ncoords, order, hcase); // check rule before recursion

  result = rule; 
  // apply rule recursively (depth) times
  if (depth > 1) {  //
    result = ruleApply(rule, coords, ncoords, rule);
    printf(" rule^2: %s\n", result);
    for (i = 3; i <= depth; i++) {
      tmpstr = ruleApply(rule, coords, ncoords, result);
      free(result); 
      result = tmpstr;
      printf(" rule^%d: %s\n", i, result); /**/
    }
  }

  outfile = stdout;  // open output file if given
  if (argc >= 6) { outfile = fopen(argv[5], "w");
    if (outfile == NULL) outfile = stdout;
  }
  fprintf(outfile, "// Generated by hilbertOrder.c using line:\n// ");
  for (i = 0; i < argc; i++) fprintf(outfile, "%s ", argv[i]);
  fprintf(outfile, "\n");

  crange = (int)(pow((double) crange, (double) depth)+0.5);
  doOutput(outfile, result, crange, coords, ncoords);
  fclose(outfile);

  if (argc >= 7) {
    outfile = fopen(argv[6], "w");
    pointOutput(outfile, coords, ncoords, result); 
  }
  return EXIT_SUCCESS;
}
