.SUFFIXES : .F .f .c .o .a .f90 .f95 .F90
########################################################################
# COSP Makefile
########################################################################

basedir=$(HOME)/github/Tools/COSPv2.0/src

obj  = $(basedir)/obj
lib = $(basedir)/lib
mod = $(basedir)/mod

# Compiler and flags
# Intel (ifort)
F90      = ifort
##F90FLAGS = -g -debug all -fpe0 -traceback -check bounds -check all -O0 
F90FLAGS = -O0 -module $(mod) -fpp
# F90_LIB  = /sw/rhel6-x64/intel/intel-17.0.2

# Portland Group (pgf90)
#F90      = pgf90
#F90FLAGS = -Mpreprocess -O0 -g
#F90_LIB  = /usr/local/pgf90

#F90      = gfortran
#F90FLAGS = -O0 -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -g
##F90FLAGS = -cpp -O0 -fbacktrace -fbounds-check


# Libraries
INC      = /usr/local/include
LIB      = /usr/local/lib

# Simulator directories
SIM_PATH     = simulator
RT_PATH      = $(SIM_PATH)/rttov
RS_PATH      = $(SIM_PATH)/quickbeam
CS_PATH      = $(SIM_PATH)/actsim
LLNL_PATH    = $(SIM_PATH)/llnl
ISCCP_PATH   = $(SIM_PATH)/icarus-scops-4.1-bsd
MISR_PATH    = $(SIM_PATH)/MISR_simulator
MODIS_PATH   = $(SIM_PATH)/MODIS_simulator
PARASOL_PATH = $(SIM_PATH)/parasol
RTTOVSTUB_PATH = $(SIM_PATH)/rttov

# RTTOV setup
RTTOV_PATH       = $(HOME)/work/Tools/RTTOV/rttov113
RTTOV_LIB_PATH   = $(RTTOV_PATH)/lib 
RTTOV_INC_PATH   = $(RTTOV_PATH)/include 
RTTOV_MOD_PATH   = $(RTTOV_PATH)/mod
RTTOV_LIBS       = -lrttov11_mw_scatt -lrttov11_brdf_atlas -lrttov11_emis_atlas \
                   -lrttov11_other -lrttov11_parallel -lrttov11_coef_io        \
                   -lrttov11_main

# List of library files that will be created
# ---------------------------------------------------------------------------------------------------------------------------------------
LIB_COSP = $(lib)/libcosp.a

FLAG_LIB_COSP = -L$(lib) -lcosp
# ---------------------------------------------------------------------------------------------------------------------------------------

# Final flag indicating the library and include paths
# ---------------------------------------------------------------------------------------------------------------------------------------
FLAGS_COSP = $(FLAG_LIB_COSP)
# ---------------------------------------------------------------------------------------------------------------------------------------

# List of all objects
# ---------------------------------------------------------------------------------------------------------------------------------------
LIST_OBJ_COSP = $(obj)/cosp_kinds.o \
		$(obj)/cosp_config.o \
		$(obj)/cosp_stats.o \
		$(obj)/cosp_constants.o \
		$(obj)/quickbeam.o \
		$(obj)/parasol.o \
		$(obj)/icarus.o \
		$(obj)/cosp_rttovSTUB.o \
		$(obj)/MISR_simulator.o \
		$(obj)/modis_simulator.o \
		$(obj)/lidar_simulator.o \
		$(obj)/cosp_calipso_interface.o \
		$(obj)/cosp_cloudsat_interface.o \
		$(obj)/cosp_isccp_interface.o \
		$(obj)/cosp_misr_interface.o \
		$(obj)/cosp_parasol_interface.o \
		$(obj)/cosp_modis_interface.o \
		$(obj)/cosp_rttov_interfaceSTUB.o \
		$(obj)/cosp.o
# ---------------------------------------------------------------------------------------------------------------------------------------

install : $(LIST_OBJ_COSP)
	ar r $(LIB_COSP) $(LIST_OBJ_COSP)

clean:	
	rm -f $(obj)/*.o $(mod)/*.mod $(lib)/*.a


# Creation of the object files
# ---------------------------------------------------------------------------------------------------------------------------------------
$(obj)/icarus.o : $(ISCCP_PATH)/icarus.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/quickbeam.o: $(RS_PATH)/quickbeam.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_calipso_interface.o : $(SIM_PATH)/cosp_calipso_interface.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_cloudsat_interface.o : $(SIM_PATH)/cosp_cloudsat_interface.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_isccp_interface.o : $(SIM_PATH)/cosp_isccp_interface.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_misr_interface.o : $(SIM_PATH)/cosp_misr_interface.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_parasol_interface.o : $(SIM_PATH)/cosp_parasol_interface.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_modis_interface.o : $(SIM_PATH)/cosp_modis_interface.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_rttov_interfaceSTUB.o : $(SIM_PATH)/cosp_rttov_interfaceSTUB.F90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -c $< -o $@

$(obj)/MISR_simulator.o : $(MISR_PATH)/MISR_simulator.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/modis_simulator.o : $(MODIS_PATH)/modis_simulator.F90 
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_rttovSTUB.o : $(RTTOVSTUB_PATH)/cosp_rttovSTUB.F90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -c $< -o $@

$(obj)/lidar_simulator.o : $(CS_PATH)/lidar_simulator.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/parasol.o : $(PARASOL_PATH)/parasol.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_constants.o :  hooks/cosp_constants.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_kinds.o : hooks/cosp_kinds.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/mo_rng.o : rng/mo_rng.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_config.o :  cosp_config.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_stats.o :  cosp_stats.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp.o :  cosp.F90
	$(F90) $(F90FLAGS) -c $< -o $@
# ---------------------------------------------------------------------------------------------------------------------------------------
