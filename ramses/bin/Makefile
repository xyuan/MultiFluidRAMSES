#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@cea.fr
#############################################################################
# compilation options
#--------------------
# number of dimensions
ifndef NDIM
	NDIM = 3
endif
# number of variables
ifndef NVAR
	NVAR = 10
endif
# position-dependant gamma
VAR_G = 10
# floats precision
ifndef NPRE
	NPRE = 8
endif
# fluid transport solver : hydro or mhd
ifndef SOLVER
	SOLVER = hydro
endif
OPTIONS = -DNDIM=$(NDIM) -DNVAR=$(NVAR) -DNPRE=$(NPRE) -DSOLVER$(SOLVER) -DVAR_G=$(VAR_G)
#############################################################################
# compilation options
#--------------------
# user-defined code modifications
ifndef PATCH
	PATCH = ../patch/accelrsn/
endif
# executable name (number of dimensions will be added)
ifndef EXEC
	EXEC = ramses
endif
SAVE = source
# compiler options
# --- No MPI, pgf90 ----------------------------------
#FCOMP = pgf90
#FFLAGS = -Mpreprocess -DWITHOUTMPI
# --- No MPI, xlf ------------------------------------
#FCOMP = xlf
#FFLAGS = -WF -DWITHOUTMPI -qfree=f90 -qsuffix=f=f90 -qsuffix=cpp=f90
# --- No MPI, f90 ------------------------------------
#FCOMP = f90
#FFLAGS = -cpp -DWITHOUTMPI
# --- No MPI, ifort ----------------------------------
#FCOMP = ifort
#FFLAGS = -cpp -DWITHOUTMPI
# --- MPI, pgf90 syntax ------------------------------
#FCOMP = mpif90 -O3
#FFLAGS = -Mpreprocess
# --- MPI, ifort syntax ------------------------------
#FCOMP = mpif90 -O3
#FFLAGS = -cpp -DNOSYSTEM
# --- MPI, ifort syntax, additional checks -----------
#FCOMP = mpif90
#FFLAGS = -O3 -g -traceback -fpe0 -ftrapuv -cpp -DNOSYSTEM
# --- MPI with openMPI on my Mac ---------------------
FCOMP = mpif90
#FFLAGS = -g -cpp -fno-underscoring -fbounds-check #-O -Wall #-Wunused #-Wuninitialized
FFLAGS = -cpp -fno-underscoring -fbounds-check -O -Wall #-Wunused #-Wuninitialized
CCOMP = mpicc
CFLAGS = #-g
#############################################################################
# MPI librairies
LIBMPI = 
#LIBMPI = -lfmpi -lmpi -lelan
LIBS = $(LIBMPI)
#############################################################################
# Sources directories are searched in this exact order
VPATH = ../../em/non-thermal:../../em/thermal:../../Blasi:../../Chevalier:$(PATCH):../$(SOLVER):../hydro:../pm:../poisson:../amr
#############################################################################
MODOBJ = amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o
AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o amr_step.o update_time.o output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o refine_utils.o nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o units.o
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o sink_particle.o feedback.o
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o output_poisson.o
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o read_hydro_params.o synchro_hydro_fine.o
ifeq ($(PATCH),../patch/accelrsn/)
	NTOBJ = beschb.o bessik.o chebev.o antinue.o antinumu.o elec.o gamma.o nue.o numu.o posi.o sigma.o em.o
	THOBJ = TeTi.o abint.o spec_elt.o exponentiation.o math.o ionisation.o recombinaison.o
	SNROBJ1 = $(NTOBJ) $(THOBJ) thermal_V1.2_module.o Blasi_V2.1_module.o Chevalier_V2.1_module.o
	SNROBJ2 = accel.o emission.o
endif
# All objects
AMRLIB = $(SNROBJ1) $(MODOBJ) $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ) $(SNROBJ2)
#############################################################################
all: ramses save
#############################################################################
ramses:	$(AMRLIB) ramses.o
	$(FCOMP) $(FFLAGS) $(OPTIONS) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
#############################################################################
%.o:%.f90
	$(FCOMP) $(FFLAGS) $(OPTIONS) -c $^ -o $@
	cp -f $^ $(SAVE)
#############################################################################
%.o:%.f
	$(FCOMP) $(FFLAGS) -ffixed-line-length-none -c $^ -o $@
	cp -f $^ $(SAVE)
#############################################################################
%.o:%.c
	$(CCOMP) $(CFLAGS) -c $^ -o $@
	cp -f $^ $(SAVE)
#############################################################################
save:
	cp -f Makefile $(SAVE)
	tar -zcf $(SAVE).tgz $(SAVE)
#############################################################################
clean :
	rm -f *.o *.mod
	rm -f $(EXEC)*
	rm -f $(SAVE)/* $(SAVE).tgz
#############################################################################
