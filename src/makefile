F90 = mpif90 -I/usr/local/software/spack/spack-views/.rhel8-icelake-202110272/uxqqj4xcjrltatqgtuoi2hp46uabtzom/intel-oneapi-mpi-2021.4.0/intel-2021.4.0/kypfgtnfzspxoby7tqy7yt6ykejpwk5n/mpi/2021.4.0/

FFLAGS = -O2
LDFLAGS = -lnetcdff

SRC = \
	tgu.f90		\
	params.f90	\
	vars.f90	\
	grow_up.f90

OBJ = $(SRC:.f90=.o)

tgu.exe : $(OBJ)
	$(F90) $(FFLAGS) -o tgu.exe $(OBJ) $(LDFLAGS)

# Main routine
tgu.o : params.o vars.o grow_up.o tgu.f90
	$(F90) $(FFLAGS) -c tgu.f90
	
# Subroutines
grow_up.o : params.o vars.o grow_up.f90
	$(F90) $(FFLAGS) -c grow_up.f90

# Modules
params.o : params.f90
	$(F90) $(FFLAGS) -c params.f90

vars.o : vars.f90
	$(F90) $(FFLAGS) -c vars.f90
