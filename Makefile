SHELL=/bin/zsh

.SUFFIXES:
FC=gfortran
CFLAGS=-Ofast -fopenmp
SRC=fortran_src
OBJ=Obj_folder

#all: simu main.o particle_management.o aux_mod.o stoermer_verlet.o
#	mkdir -p $(OBJ)
#
#.PHONY: all

dir_guard=mkdir -p $(OBJ)

simu: main.o aux_mod.o stoermer_verlet.o particle_management.o
	$(FC) $(CFLAGS) -o simu $(OBJ)/main.o $(OBJ)/aux_mod.o $(OBJ)/stver_mod.o $(OBJ)/part_manag.o

main.o: $(SRC)/program_particle_simulation.f08 aux_mod.o stoermer_verlet.o particle_management.o 
	$(FC) $(CFLAGS) -I$(OBJ) -c $(SRC)/program_particle_simulation_openmp.f08 -o $(OBJ)/main.o
  
particle_management.o: aux_mod.o stoermer_verlet.o
	$(FC) $(CFLAGS) -J$(OBJ) -I$(OBJ) -c $(SRC)/module_particle_management.f08 -o $(OBJ)/part_manag.o

aux_mod.o: $(SRC)/module_auxiliary_routines.f08
	$(dir_guard)
	$(FC) $(FLAGS) -J$(OBJ) -c $(SRC)/module_auxiliary_routines.f08 -o $(OBJ)/aux_mod.o

stoermer_verlet.o: $(SRC)/module_stoermer_verlet_method.f08
	$(dir_guard)
	$(FC) $(CFLAGS) -J$(OBJ) -c $(SRC)/module_stoermer_verlet_method.f08 -o $(OBJ)/stver_mod.o

clean:
	rm -r $(OBJ)
