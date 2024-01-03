# list of executable files to produce

BFBD = getPost
DFBD = getAIC
GFBD = testAIC

# .o files necessary to build the executables

OBJ_BFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o GetModel.o 
OBJ_DFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o GetAIC.o
OBJ_GFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o TestAIC.o

########### MODIFY ONLY TO CHANGE OPTIONS ############

all = DIST DIBI TREE ESTI DRAW

# compiler and its options
CC = gcc
#CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT -g 
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -I/usr/include/gsl -I/usr/share/R/include -g -fopenmp

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(BFBD) $(DFBD) $(GFBD)

# build the executable
$(BFBD): $(OBJ_BFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g
# build the executable
$(DFBD): $(OBJ_DFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g
# build the executable
$(GFBD): $(OBJ_GFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lcairo -fopenmp  -g

############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release
