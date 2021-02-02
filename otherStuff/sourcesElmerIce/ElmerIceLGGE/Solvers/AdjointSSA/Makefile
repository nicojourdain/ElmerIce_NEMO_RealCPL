
EXECUTABLES := AdjointSSASolvers

OBJ=objects

_SOURCES := AdjointSSA_SSASolver.F90 \
	    AdjointSSA_GradientSolver.F90 \
	    AdjointSSA_AdjointSolver.F90 \
	    AdjointSSA_CostDiscSolver.F90 \
	    AdjointSSA_CostContSolver.F90 \
	    AdjointSSA_CostRegSolver.F90 \
	    AdjointSSA_CostFluxDivSolver.F90 \
	    USFs_VCHANGE.F90

SOURCES=$(patsubst $(SRC)/%, % ,$(_SOURCES))

_SOLVERS := $(patsubst %.F90, %.o ,$(SOURCES))
SOLVERS = $(patsubst %, $(OBJ)/% ,$(_SOLVERS))

compile: $(EXECUTABLES)

AdjointSSASolvers: $(SOLVERS)
	elmerf90 $^ -o $@

$(SOLVERS): $(OBJ)/%.o: $(SRC)/%.F90
	@if [ ! -d "./$(OBJ)" ]; then \
	     mkdir $(OBJ) ;\
	fi
	elmerf90 -c $< -o $@


clean: 
	rm -f $(EXECUTABLES)
	rm -f $(SOLVERS)
