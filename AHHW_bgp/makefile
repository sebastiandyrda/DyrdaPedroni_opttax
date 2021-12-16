
# =========================================================
# MAKEFILE
# =========================================================
FCOMP    = mpiifort
FCFLAGS  = -O3 -heap-arrays -traceback -mkl -g
PROGRAM  = ahhw_bgp

MOD      = MOD_PARAMETERS.o

SUBR     = DECISION_RULES.o DISTRIBUTION.o  DISTRIBUTION_Q.o GRID.o\
           INTERP.o MAIN_AHHW.o MAIN_FUNCTION.o \
		   SPLINE.o STATIONARY_EQUIL.o TOOLBOX.o \
		   POLICY_FUNCTION_ITERATION.o \
		   PARTIAL_DERIVATIVES.o Q_TO_KAPPA.o RESIDUAL.o
  

OBJ      = $(MOD) $(SUBR)

# Set Phony targets
.PHONY : all clean

#@echo: Compiling the main program
all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(FCOMP) $(FCFLAGS) -o $@ $(OBJ)
	
MOD_PARAMETERS.mod: MOD_PARAMETERS.o MOD_PARAMETERS.f90
	$(FCOMP) -c MOD_PARAMETERS.f90	
	
%.o: %.f90
	$(FCOMP) -c $(FCFLAGS) $<

%.o: %.f
	$(FCOMP) -c $(FCFLAGS) $<
	
	
# Cleaning everything
#@echo: Cleaning files
clean:
	rm $(OBJ)
	
@echo: Compilation COMPLETED

