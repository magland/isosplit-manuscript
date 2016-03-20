# make MEX interface to C OpenMP dbscan Patwary code.
# Barnett 4/17/15

# To do: check if headers .h trigger re-compilation

MWRAP = ../mwrap-0.33.3/mwrap

# name of C file that MWrap will make, and MEX file
GATEWAY = db

# name of .mw file, and also .c file that MEX calls
HEAD = dbscanlib

CC = g++
CFLAGS= -c -fopenmp -O3
LDFLAGS=
LIBS= -fopenmp -O3
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)
APP = omp_dbscan

.c.o:
	gcc -fPIC -fopenmp -O3 -c $<
.cpp.o:
	g++ -fPIC -fopenmp -O3 -c $<

# for now make the cmdline version...
default: $(APP)

# cmd line version...
$(APP):	$(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(APP) $(LIBS)

# not finished...
mex:	$(HEAD).o $(HEAD).mw
# note -mb here makes separate .m file for each @ function in the .mw ...
	$(MWRAP) -list -mex $(GATEWAY) -mb $(HEAD).mw
	$(MWRAP) -mex $(GATEWAY) -c $(GATEWAY).c $(HEAD).mw
	mex $(GATEWAY).c $(HEAD).o -lm -lgomp

clean:
	rm -f *.o *.mex* $(APP)
