##############################################
# Makefile for Elmer/Ice LGGE library
#	---------------------------
# make -> Create LGGESolvers.$(EXT) LGGEUSF.$(EXT)	
# make install -> Create directory lib and copy libraries
# make clean -> delete *.o
# make cleanall ->delete *.o and lib directory
# #############################################

OS :=$(shell uname -s)
ifeq ($(OS), Linux)
	EXT:=so
else
	EXT:=dylib
endif

ELMERICE_HOME ?= $(ELMER_HOME)/share/elmersolver

PROJ_INC=/scratch/cnt0021/gge6035/SHARED/local/include
PROJ_LIB=/scratch/cnt0021/egi6035/SHARED/local/lib

MKDIR_P = mkdir -p

ifndef $(ElmerIceLGGE_HOME)
	ElmerIceLGGE_HOME=.
endif

ElmerIceLGGE_LIB=$(ElmerIceLGGE_HOME)/lib

EXECUTABLES := LGGESolvers.$(EXT) LGGEUSF.$(EXT)

include Solvers/Makefile
include UserFunctions/Makefile

solverobjects := $(filter %.o,$(patsubst %.F90,%.o,$(SOLVERSOURCE))) $(filter %.o,$(patsubst %.f,%.o,$(SOLVERSOURCE)))
usfobjects := $(patsubst %.F90,%.o,$(USFSOURCE))

all: compile

compile: $(EXECUTABLES)

LGGESolvers.$(EXT): $(solverobjects) 
	elmerf90 $(solverobjects) -o LGGESolvers.$(EXT)

LGGEUSF.$(EXT): $(usfobjects)
	elmerf90 $(usfobjects) -o LGGEUSF.$(EXT)

install: directories
	cp -f LGGESolvers.$(EXT) $(ElmerIceLGGE_LIB)/
	cp -f LGGEUSF.$(EXT) $(ElmerIceLGGE_LIB)/
	@echo '*-----------------------------------------*'
	@echo '  Add $(ElmerIceLGGE_LIB) to LD_LIBRARY_PATH  '
	@echo '*-----------------------------------------*'
#For root user:
#       cp -f LGGESolvers.$(EXT) $(ELMERICE_HOME)/lib/
#       cp -f LGGEUSF.$(EXT) $(ELMERICE_HOME)/lib/

directories:
	${MKDIR_P} $(ElmerIceLGGE_LIB)

clean:
	rm -f  Solvers/*.o UserFunctions/*.o

cleanall: clean
	rm -f $(EXECUTABLES) $(ElmerIceLGGE_LIB)/*.$(EXT)
	rm -r $(ElmerIceLGGE_LIB)

%.o: %.F90 
	echo elmerf90 -c $< -o $@
	elmerf90 -c $< -o $@

%.o: %.f
	elmerf90 -c $< -o $@

#Solvers/ToLonLat.o: Solvers/ToLonLat.F90
#	elmerf90 -I${PROJ_INC} $^ -o $@  -L${PROJ_LIB} -lproj -lfortranc -lfortrangis



