
# =========================================================
# MAKEFILE
# =========================================================



#FCOMP   = ifort
FCOMP    = mpiifort
#FCFLAGS  = -traceback -O3 -heap-arrays -fast -qopenmp
FCFLAGS  = -O3 -heap-arrays -traceback -mkl
#FCFLAGS  = -g -heap-arrays -check all -traceback -mkl
#PROGRAM  = mainrun
PROGRAM  = ahhw_bgp

MOD      = MOD_PARAMETERS.o KIND_MODULE.o BOBYQA_MODULE.o
		   
SUBR     = DECISION_RULES.o DISTRIBUTION.o  DISTRIBUTION_Q.o GRID.o\
           INTERP.o MAIN_MPI.o MAIN_FUNCTION.o \
		   SPLINE.o STATIONARY_EQUIL.o STATIONARY_EQUIL_DEBT.o TOOLBOX.o \
		   toms659.o POLICY_FUNCTION_ITERATION.o \
		   PARTIAL_DERIVATIVES.o Q_TO_KAPPA.o RESIDUAL.o
  

OBJ      = $(MOD) $(SUBR)

# Set Phony targets
.PHONY : all clean

#@echo: Compiling the main program
all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(FCOMP) $(FCFLAGS) -o $@ $(OBJ)
	
	
BOBYQA_MODULE.mod: BOBYQA_MODULE.o KIND_MODULE.o BOBYQA_MODULE.f90
	$(FCOMP) -c BOBYQA_MODULE.f90
	
MOD_PARAMETERS.mod: MOD_PARAMETERS.o MOD_PARAMETERS.f90
	$(FCOMP) -c MOD_PARAMETERS.f90
	
KIND_MODULE.mod: KIND_MODULE.o KIND_MODULE.f90
	$(FCOMP) -c KIND_MODULE.f90	
	
%.o: %.f90
	$(FCOMP) -c $(FCFLAGS) $<

%.o: %.f
	$(FCOMP) -c $(FCFLAGS) $<
	
	
	
# Cleaning everything
#@echo: Cleaning files
clean:
	rm $(OBJ)
	
#@echo: Compilation COMPLETED

