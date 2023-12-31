# list of executable files to produce
AFBD = afbd
BFBD = bfbd
CFBD = cfbd
DFBD = dfbd
EFBD = efbd
FFBD = ffbd
GFBD = gfbd

# .o files necessary to build the executables

OBJ_AFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o TestModel.o 
OBJ_BFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o GetModel.o 
OBJ_CFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o TestAIC.o 
OBJ_DFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o GetAIC.o
OBJ_EFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o GetDiffAIC.o
OBJ_FFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o GetMultiAIC.o
OBJ_GFBD =  MyR.o MyRandom.o SimulPiecewise.o PiecewiseModel.o Utils.o Tree.o Fossil.o FossilInt.o MCMCImportanceSamplingFossilInt.o FBDDensity.o TestAICBis.o

########### MODIFY ONLY TO CHANGE OPTIONS ############

all = DIST DIBI TREE ESTI DRAW

# compiler and its options
CC = gcc
#CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -D_REENTRANT -g 
CFLAGS = -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -I/usr/include/cairo -I/usr/include/gsl -I/usr/share/R/include -g -fopenmp

# linker and its options
LD = $(CC)

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(DIST) $(TREE) $(ESTI) $(DRAW)

# build the executable
$(AFBD): $(OBJ_AFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g
# build the executable
$(BFBD): $(OBJ_BFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g
# build the executable
$(CFBD): $(OBJ_CFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g
# build the executable
$(DFBD): $(OBJ_DFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g

# build the executable
$(EFBD): $(OBJ_EFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/  -lcairo -fopenmp  -g
# build the executable
$(FFBD): $(OBJ_FFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lcairo -fopenmp  -g
# build the executable
$(GFBD): $(OBJ_GFBD)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lcairo -fopenmp  -g
# build the executable
$(DIST): $(OBJ_DIST)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -g
	
# build the executable
$(DIBX): $(OBJ_DIBX)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -g

# build the executable
$(DIBT): $(OBJ_DIBT)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -g

# build the executable
$(DIAG): $(OBJ_DIAG)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -g

# build the executable
$(DIBI): $(OBJ_DIBI)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo  -L/usr/local/lib -L/usr/lib/  -lpthread  -g
# build the executable
$(DIBA): $(OBJ_DIBA)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo  -L/usr/local/lib -L/usr/lib/  -g
# build the executable
$(DIBB): $(OBJ_DIBB)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo  -L/usr/local/lib -L/usr/lib/  -g

# build the executable
$(TREE): $(OBJ_TREE)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  
# build the executable
$(TIME): $(OBJ_TIME)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm  -L/usr/local/lib -L/usr/lib/ 
# build the executable
$(TTTT): $(OBJ_TTTT)
	$(LD) -g  $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo -L/usr/local/lib -L/usr/lib/ 

# build the executable
$(BIS): $(OBJ_BIS)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo  -L/usr/local/lib -L/usr/lib/  -g

# build the executable
$(MIS): $(OBJ_MIS)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo  -L/usr/local/lib -L/usr/lib/  -g

# build the executable
$(ESTI): $(OBJ_ESTI)
	$(LD) $^ -o $@  -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lnlopt  -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3

# build the executable
$(ESTM): $(OBJ_ESTM)
	$(LD) $^ -o $@  -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lnlopt  -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  -O3

# build the executable
$(DRAW): $(OBJ_DRAW)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo -O3
	
# build the executable
$(TEST): $(OBJ_TEST)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm  -lgsl -lcairo -O3
# build the executable
$(DENS): $(OBJ_DENS)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm 
# build the executable
$(CONF): $(OBJ_CONF)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm 
# build the executable
$(SAMP): $(OBJ_SAMP)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lgsl -lgslcblas -lcairo
# build the executable
$(SINT): $(OBJ_SINT)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -lcairo
	
# build the executable
$(INTE): $(OBJ_INTE)
	$(LD) $^ -o $@ -L/usr/local/lib -lgsl -lgslcblas -L/usr/lib/ -lm
# build the executable
$(GEOS): $(OBJ_GEOS)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -g
# build the executable
$(REVG): $(OBJ_REVG)
	$(LD) $^ -o $@ -L/usr/local/lib -L/usr/lib/ -lm -g

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
