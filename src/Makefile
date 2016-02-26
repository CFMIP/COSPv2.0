.SUFFIXES : .F .f .c .o .a .f90 .f95 .F90
########################################################################
# COSP Makefile
########################################################################
# Compiler and flags
F90      = ifort
F90FLAGS = -fpp -O3

# Build cosp as a library to link against.
TARGET  = lib/libcosp.a
SRC_PATH = src/
OBJ_PATH = obj/

# Libraries
F90_LIB  = /usr/local/ifort
INC      = /usr/local/include
LIB      = /usr/local/lib
NCDF_INC    = $(F90_LIB)/include
NCDF_LIB    = $(F90_LIB)/lib
CMOR_INC    = $(F90_LIB)/include
CMOR_LIB    = $(F90_LIB)/lib
UDUNITS_LIB = $(F90_LIB)/lib
UDUNITS_INC = $(F90_LIB)/include
UUID_LIB    = $(F90_LIB)/lib
UUID_INC    = $(F90_LIB)/include

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

# RTTOV setup
RTTOV_PATH       = /Projects/Clouds/dswales/RTTOV/rttov_11.3
RTTOV_LIB_PATH   = $(RTTOV_PATH)/lib 
RTTOV_INC_PATH   = $(RTTOV_PATH)/include 
RTTOV_MOD_PATH   = $(RTTOV_PATH)/mod
RTTOV_LIBS       = -lrttov11_coef_io_789 -lrttov11_coef_io -lrttov11_main \
                   -lrttov11_mw_scatt -lrttov11_mw_scatt_coef -lrttov11_other

OBJS =  cosp_kinds.o cosp_constants.o cosp_cloudsat_interface.o cosp_config.o  \
        prec_scops.o cosp.o cosp_stats.o quickbeam.o cosp_optics.o\
	array_lib.o mo_rng.o quickbeam_optics.o scops.o math_lib.o optics_lib.o\
        parasol.o mrgrnk.o lidar_simulator.o icarus.o cosp_calipso_interface.o \
        cosp_errorHandling.o cosp_isccp_interface.o cosp_misr_interface.o      \
        MISR_simulator.o cosp_modis_interface.o modis_simulator.o              \
	cosp_rttov_interface.o cosp_rttov11.3.o cosp_parasol_interface.o

# To create cosp library, use this.
#$(TARGET): $(OBJ_PATH)$(OBJS)
#	rm -f $@
#	ar cr $@ $^



%.o: $(SIM_PATH)/%.f90
	@echo $(F90) $(F90FLAGS) -c -I$(NCDF_INC) -I$(CMOR_INC) $<
	$(F90) $(F90FLAGS) -c -I$(NCDF_INC) -I$(CMOR_INC) $<
	@echo "-----------------------------"

%.o: %.f90
	@echo $(F90) $(F90FLAGS) -c -I$(NCDF_INC) -I$(CMOR_INC) $<
	$(F90) $(F90FLAGS) -c -I$(NCDF_INC) -I$(CMOR_INC) $<
	@echo "-----------------------------"

# Dependencies
cosp.o                    : cosp_kinds.o cosp_modis_interface.o cosp_constants.o prec_scops.o      \
                            cosp_rttov_interface.o cosp_misr_interface.o cosp_isccp_interface.o    \
                            cosp_calipso_interface.o cosp_cloudsat_interface.o cosp_stats.o        \
			    cosp_parasol_interface.o cosp_rttov11.3.o cosp_rttov_interface.o       \
                            quickbeam.o MISR_simulator.o lidar_simulator.o parasol.o icarus.o      \
			    quickbeam_optics.o cosp_errorHandling.o
cosp_config.o             : cosp_kinds.o
cosp_stats.o              : cosp_kinds.o cosp_config.o
cosp_calipso_interface.o  : cosp_kinds.o lidar_simulator.o
cosp_cloudsat_interface.o : cosp_kinds.o cosp_config.o quickbeam.o
cosp_isccp_interface.o    : cosp_kinds.o icarus.o
cosp_misr_interface.o     : cosp_kinds.o 
cosp_modis_interface.o    : cosp_kinds.o cosp_config.o modis_simulator.o
cosp_rttov_interface.o    : cosp_kinds.o cosp_config.o cosp_rttov11.3.o
cosp_parasol_interface.o  : cosp_kinds.o 
cosp_rttov11.3.o          : cosp_kinds.o cosp_config.o
optics_lib.o	          : cosp_kinds.o
array_lib.o   	          : cosp_kinds.o
math_lib.o                : cosp_kinds.o array_lib.o mrgrnk.o
MISR_simulator.o          : cosp_kinds.o cosp_config.o cosp_stats.o
modis_simulator.o         : cosp_kinds.o cosp_config.o cosp_stats.o
lidar_simulator.o         : cosp_kinds.o cosp_config.o cosp_stats.o
icarus.o	          : cosp_kinds.o cosp_constants.o cosp_optics.o cosp_stats.o
quickbeam.o               : cosp_kinds.o cosp_config.o cosp_stats.o
quickbeam_optics.o        : cosp_kinds.o cosp_config.o array_lib.o math_lib.o optics_lib.o         \
                            cosp_constants.o quickbeam.o
mrgrnk.o                  : cosp_kinds.o
prec_scops.o	          : cosp_kinds.o
cosp_optics.o             : cosp_kinds.o cosp_constants.o 
parasol.o                 : cosp_kinds.o cosp_config.o cosp_constants.o
cosp_errorHandling.o      : cosp_kinds.o

install:
	make
	mv *.mod *.o $(OBJ_PATH)
clean:	
	rm -f $(OBJ_PATH)$(OBJS) *.mod *.o fort.*

icarus.o : $(ISCCP_PATH)/icarus.f90
	$(F90) $(F90FLAGS) -c $<

quickbeam.o: $(RS_PATH)/quickbeam.f90
	$(F90) $(F90FLAGS) -c $<

MISR_simulator.o : $(MISR_PATH)/MISR_simulator.f90
	$(F90) $(F90FLAGS) -c $<

modis_simulator.o : $(MODIS_PATH)/modis_simulator.F90 
	$(F90) $(F90FLAGS) -c $<

cosp_rttov11.3.o : $(RT_PATH)/cosp_rttov11.3.F90 
	$(F90) $(F90FLAGS) -c -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) $<

lidar_simulator.o : $(CS_PATH)/lidar_simulator.f90
	$(F90) $(F90FLAGS) -c $<

parasol.o : $(PARASOL_PATH)/parasol.f90
	$(F90) $(F90FLAGS) -c $<

cosp_optics.o :  optics/cosp_optics.f90
	$(F90) $(F90FLAGS) -c $<

quickbeam_optics.o :  optics/quickbeam_optics.f90
	$(F90) $(F90FLAGS) -c $<

cosp_constants.o :  hooks/cosp_constants.f90
	$(F90) $(F90FLAGS) -c $<

cosp_kinds.o : hooks/cosp_kinds.f90
	$(F90) $(F90FLAGS) -c $<

cosp_errorHandling.o : hooks/cosp_errorHandling.f90
	$(F90) $(F90FLAGS) -c $<

mo_rng.o : rng/mo_rng.f90
	$(F90) $(F90FLAGS) -c $<

prec_scops.o : subcol/prec_scops.f90
	$(F90) $(F90FLAGS) -c $<

scops.o : subcol/scops.f90
	$(F90) $(F90FLAGS) -c $<

cosp_config.o :  cosp_config.f90
	$(F90) $(F90FLAGS) -c $<

cosp.o :  cosp.f90
	$(F90) $(F90FLAGS) -c $<

cosp_stats.o :  cosp_stats.f90
	$(F90) $(F90FLAGS) -c $<

array_lib.o : optics/array_lib.f90
	$(F90) $(F90FLAGS) -c $<

math_lib.o : optics/math_lib.f90
	$(F90) $(F90FLAGS) -c $<

mrgrnk.o : optics/mrgrnk.f90
	$(F90) $(F90FLAGS) -c $<

optics_lib.o : optics/optics_lib.f90
	$(F90) $(F90FLAGS) -c $<
