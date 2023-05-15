// AUDIT routines for d3.c

/* Tetrahedra are groups of four corners in order of increasing
vertex index, except that the first two may be swapped to ensure
that the orientation determinant is positive.
I.e., take the lex smallest alternating permutation with positive sign.
There must always be an odd number of swaps between two permutations;
we swap the first two if necessary to achieve this.

The table indoff has the index offset for where
the vertex of index i will be found in the tetrahedron opposite c.
Vertex s[BASECORNER(c)+i].v will be at s[s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]].v
(except that s[c].v and s[s[c].opp]].v are different, and opposite sides of common pl).
Note that indoff[*][j] uses each offset -j:4-j exactly once,
that indoff[i][j][i] = 0 for all possible i and j (so s[c].v is opposite of s[s[c].opp].v),
and that the tetra will never change, TETRA(s[c].opp) == TETRA(CORNERINOPP(i,c)).
*/

// This is the table of where the indices go; used in ASSERTS only.
// Replacing  \ With V falling at position c.opp:
//   corner c  \0ABCD 1ABCD 2ABCD 3ABCD
//       A    0: xxxx  BVCD  CBVD  BCDV
//       B    1: VACD  xxxx  ACVD  CADV
//       C    2: vBAD  AVBD  BAVD  ABDV
//       D    3: Vabc  bvac  ABVC  BACV
//
// Offsets from c.opp  I=-1, Z=-2, B = -3, H = -4
// Replacing \0ABCD 1ABCD 2ACBD 3ABCD
//     A    0: xxxx  0I12  0IZ1  0BZI
//     B    1: 1023  xxxx  Z0I1  Z0BI
//     C    2: 1203  I102  IZ01  BZ0I
//     D    3: 1230  1I20  ZI10  ZBI0
// Get these by subtracting 01234 from columns in order
#define CORNERINOPP(i, c) (this->s[c].opp + indoff[INDEX(c)][INDEX(this->s[c].opp)][i])
static const short indoff[4][4][4] = {/*        0ABCD       1ABCD         2ACBD          3ABCD     */
									  /* 0: */ {{5, 5, 5, 5}, {0, -1, 1, 2}, {0, -1, -2, 1}, {0, -3, -2, -1}},
									  /* 1: */ {{1, 0, 2, 3}, {5, 5, 5, 5}, {-2, 0, -1, 1}, {-2, 0, -3, -1}},
									  /* 2: */ {{2, 1, 0, 3}, {-1, 1, 0, 2}, {-1, -2, 0, 1}, {-3, -2, 0, -1}},
									  /* 3: */ {{1, 2, 3, 0}, {1, -1, 2, 0}, {-2, -1, 1, 0}, {-2, -3, -1, 0}}};

void cornerPrint(const pd3stateType this, int c)
{
	printf("%%%3d(%2d,%1d)%2d=", c, TETRA(c), INDEX(c), this->s[c].v - this->vert);
	if ((this->s[c].v >= this->vert) && (this->s[c].v < this->vert + MAXVERT))
		printf("(%d %5.0f %5.0f %5.0f) opp:%4d(%3d,%2d) \n", (infiniteV(this->s[c].v, this->vert) ? 0 : 1),
			   this->s[c].v->x, this->s[c].v->y, this->s[c].v->z, this->s[c].opp,
			   TETRA(this->s[c].opp), INDEX(this->s[c].opp));
	(void)fflush(stdout);
}

void cornerPrint4(const pd3stateType this, int c)
{
	int b, j, k;
	psphereType sp;
	b = BASECORNER(c);
	if (this->sph != NULL)
	{
		sp = this->sph + TETRA(b);
		printf("disp('Sphere(%d) = <%5.0f %5.0f %5.0f %5.0f>')\n", TETRA(b),
			   sp->x, sp->y, sp->z, sp->sq);
	}
	printf("DetCheckH([");
	for (k = 0; k < 4; k++)
		printf(" %d %5.0f %5.0f %5.0f %5.0f; %% %d\n",
			   (infiniteV(this->s[b + k].v, this->vert) ? 0 : 1), this->s[b + k].v->x, this->s[b + k].v->y, this->s[b + k].v->z, this->s[b + k].v->sq, this->s[b + k].v->index);
	printf("]);\n");
	for (j = 0; j < 4; j++)
		cornerPrint(this, b + j);
}

// when keeping statistics...
int maxLocate = 0, locateSideCnt = 0, sphereCnt = 0, startTetraCnt = 0, freeTetraCnt = 0, inSphereCnt = 0, randbitCount = 0;

#ifndef NOASSERT
#define auditCorners(this, c) auditCornersAux(this, c);
#else /* NOASSERT */
#define auditCorners(this, c)
#endif /* NOASSERT */

void auditCornersAux(const pd3stateType this, int SphereCheck)
{
	int p, b, c, i, j, k, guard;
	double d;
	psphereType sp;
	ppointType vv, vp;

	for (p = 0; p < this->maxTetra; p++)
	{
		if (DEAD(p))
			continue; // don't audit tetras on free list
		b = CORNER(p, 0);
		for (c = b; c < CORNER(p, 4); c++)
		{						// per corner checks
			i = this->s[c].opp; // check opposite
			if (this->s[i].opp != c)
			{
				printf("%%AUDIT: wrong opp.opp \n");
				cornerPrint(this, c);
				cornerPrint(this, i);
			}
			if (this->s[c].v == this->s[i].v)
			{
				printf("%%AUDIT: Same vertex  s[%d(%d,%d)].v %d == opp[%d(%d,%d)].v %d \n",
					   c, TETRA(c), INDEX(c), this->s[c].v - this->vert, i, TETRA(i), INDEX(i), this->s[i].v - this->vert);
				cornerPrint4(this, c);
				cornerPrint4(this, i);
			}

			for (j = 0; j < 4; j++)
				if ((j != INDEX(c)) && (CORNERINOPP(j, c) < 3))
					if (this->s[BASECORNER(c) + j].v != this->s[CORNERINOPP(j, c)].v)
					{
						printf("%%AUDIT:Bad auditCornerInOpp(%d,%d) = %d  since vertex %d != %d\n",
							   j, c, CORNERINOPP(j, c), this->s[BASECORNER(c) + j].v - this->vert, this->s[CORNERINOPP(j, c)].v - this->vert);
						cornerPrint4(this, BASECORNER(c));
						cornerPrint4(this, BASECORNER(CORNERINOPP(j, c)));
						break;
					}

			for (j = 0; j < 4; j++)
			{
				k = CORNERINOPP(j, c);
				if ((TETRA(k) != TETRA(this->s[c].opp)) || (this->s[b + j].v != this->s[k].v))
				{
					if (TETRA(k) != TETRA(this->s[c].opp))
					{
						printf("%%AUDIT: CORNERINOPP(%d,%d) ==>%d: Accessing [%d:%d][%d:%d][%d]\n",
							   j, c, k - this->s[c].opp, c, INDEX(c), this->s[c].opp, INDEX(this->s[c].opp), j);
						k = this->s[c].opp;
					}
					else
					{
						if (b + j == c)
							if (k == this->s[c].opp)
								continue; // these vertices are supposed to differ; don't flag them
							else
								printf("%%AUDIT: CORNERINOPP(%d,%d) says %d(%d,%d) and %d(%d,%d) shouldn't happen\n",
									   j, c, b + j, TETRA(b), j, k, TETRA(k), INDEX(k));
						else
						{
							printf("%%AUDIT: CORNERINOPP(%d,%d) says %d(%d,%d) and %d(%d,%d) should agree\n",
								   j, c, b + j, TETRA(b), j, k, TETRA(k), INDEX(k));
							cornerPrint4(this, c);
							cornerPrint4(this, b + j);
							cornerPrint4(this, k);
							ASSERT(TETRA(this->s[c].opp) == TETRA(CORNERINOPP(i, c)), "CORNERINOP screws up tetras");
						}
					}
					cornerPrint4(this, b);
					cornerPrint4(this, BASECORNER(k));
					break;
				}
			}

			// check sphere opposite corner c
			if (SphereCheck > 0) // check sphere opposite corner
			{
				k = this->s[c].opp;
				sp = this->sph + TETRA(k);
				if (infiniteV(this->s[c].v, this->vert))
				{
					d = spdotInf(sp, this->s[c].v);
				}
				else
				{
					vp = this->s[CORNER(TETRA(k), 3)].v; // subtract from this
					d = spdot(sp, this->s[c].v, vp);
				}
				if (d < 0)
				{
					printf("disp('AUDIT: corner %d v%d in sphere %d(%d) =%5.0f');\n",
						   c, this->s[c].v - this->vert, TETRA(k), k, d);
					cornerPrint(this, c); // print vertex
					cornerPrint4(this, k);
				}
			}
		}

		if (SphereCheck > 0) // check sphere sqs (orient dets)
		{
			sp = this->sph + p; // only spheres using pt at infty have sq==0; none have sq < 0.
			b = CORNER(p, 0);
			if (sp->sq < 0 || (sp->sq == 0 && !infiniteV(this->s[b].v, this->vert) && !infiniteV(this->s[b + 1].v, this->vert)))
			{
				printf("disp('AUDIT: sq<=0 in tetra %d(%d) =%5.0f'); \n", CORNER(p, 0), p, sp->sq);
				cornerPrint4(this, b);
			} /**/
			  /*	if (pv-this->vert > 860 || ((pv-this->vert) % 100 == 0))
				  for (vv = this->vert; vv < pv; vv++) { // Check all vertices against all spheres
				  if (vv==this->s[b+0].v) continue;
				  if (vv==this->s[b+1].v) continue;
				  if (vv==this->s[b+2].v) continue;
				  if (vv==this->s[b+3].v) continue;
				  if (vv==this->s[b+4].v) continue;
				  d = spdot(sp, vv);//	  d = spdot(sp, this->s[c].v);
				  if (d < 0) {
				  printf("disp('AUDIT: vertex v%d in sphere %d(%d) =%5.0f');\n", vv-this->vert, p, b, d);
				  printf("   (%5d; %5.0f %5.0f %5.0f; %10f)\n",
				  vv->index, vv->x, vv->y, vv->z, vv->sq);
				  cornerPrint4(this, b);
				  }
				  }/**/
		}
	}
}
