# 3d delaunay    Jack Snoeyink Aug 03
# 

CXX = gcc
CFLAGS =  -O3 $-DNOASSERT # -DSTATS # -g
LFLAGS = -lm
HFILES = pdbreader.h d3permute.h d3.h delaunay3.h
SRCS = pdbreader.c d3permute.c d3.c 
OBJS = pdbreader.o d3permute.o d3.o 
EXECS = test.exe hilbOrder.exe

all: test hilbOrder

%.o: 	%.c %.h delaunay3.h
	$(CXX) $(CFLAGS) -c $*.c 

test:	$(OBJS) 
	$(CXX) $(OBJS) test.c $(LFLAGS) -o test 

hilbOrder:	hilbertOrder.c 
	$(CXX) $(CFLAGS) hilbertOrder.c -o hilbOrder $(LFLAGS)

nuke: 
	@rm -rf $(OBJS) ii_files $(EXECS) *~

clean: 
	@rm -rf $(OBJS) ii_files *~
