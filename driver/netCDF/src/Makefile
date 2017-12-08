.SUFFIXES : .F .f .c .o .a .f90 .f95 .F90
########################################################################
# Makefile for cosp2_test offline driver
########################################################################
# Intel (ifort)
F90      = ifort
F90FLAGS = -O3 -traceback
F90_LIB  = /usr/local/ifort
#F90_LIB   = /scratch/bmcinnes/usr/local/ifort
# Portland Group (pgf90)
#F90      = pgf90
#F90FLAGS = -Mpreprocess -O0 -g
#F90_LIB  = /usr/local/pgf90
# Gfortran (gfortran)
#F90      = gfortran
#F90FLAGS = -cpp -O3 -fbacktrace -fbounds-check
#F90_LIB  = /usr/local/gfortran

INC = /usr/local/include
LIB = /usr/local/lib

NCDF_INC    = $(F90_LIB)/include
NCDF_LIB    = $(F90_LIB)/lib

# COSP source directory
HOME_PATH = /home/dswales/Projects/COSP/COSPv2.0
SRC_F90   = ${HOME_PATH}/src
OBJ_PATH  = ${HOME_PATH}/obj
PROG      = cosp2_test

# RTTOV setup
RTTOV_PATH       =
RTTOV_LIB_PATH   =
RTTOV_INC_PATH   =
RTTOV_MOD_PATH   = 
RTTOV_LIBS       =

# Object files
OBJS =  $(OBJ_PATH)/cosp_kinds.o  $(OBJ_PATH)/cosp_constants.o                  \
	$(OBJ_PATH)/cosp_cloudsat_interface.o $(OBJ_PATH)/cosp_config.o         \
        $(OBJ_PATH)/cosp.o $(OBJ_PATH)/cosp_stats.o cosp_optics.o               \
	$(OBJ_PATH)/quickbeam.o $(OBJ_PATH)/parasol.o                           \
	$(OBJ_PATH)/lidar_simulator.o $(OBJ_PATH)/icarus.o                      \
	$(OBJ_PATH)/cosp_calipso_interface.o $(OBJ_PATH)/cosp_isccp_interface.o \
	$(OBJ_PATH)/cosp_misr_interface.o $(OBJ_PATH)/MISR_simulator.o          \
	$(OBJ_PATH)/cosp_modis_interface.o $(OBJ_PATH)/modis_simulator.o        \
	$(OBJ_PATH)/cosp_rttov_interfaceSTUB.o $(OBJ_PATH)/cosp_rttovSTUB.o     \
	$(OBJ_PATH)/cosp_parasol_interface.o mo_rng.o scops.o cosp_utils.o      \
	cosp2_io.o cosp_errorHandling.o mrgrnk.o prec_scops.o array_lib.o     \
	quickbeam_optics.o math_lib.o optics_lib.o

all: $(PROG)
	mv *.mod *.o ../obj/

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(PROG).f90 $(OBJS) \
	-I$(NCDF_INC) -L${NCDF_LIB} -lnetcdff -lnetcdf -L${LIB} -lnetcdf \
	-I$(OBJ_PATH) \
	-o $(PROG)

$(PROG).o          : $(OBJ_PATH)/cosp_config.o $(OBJ_PATH)/cosp.o cosp2_io.o   \
                     $(OBJ_PATH)/quickbeam.o
cosp2_io.o          : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp.o                 \
                     $(OBJ_PATH)/cosp_config.o
cosp_utils.o       : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o
scops.o            : $(OBJ_PATH)/cosp_kinds.o mo_rng.o cosp_errorHandling.o
mo_rng.o           : $(OBJ_PATH)/cosp_kinds.o
optics_lib.o	   : $(OBJ_PATH)/cosp_kinds.o cosp_errorHandling.o
array_lib.o   	   : $(OBJ_PATH)/cosp_kinds.o cosp_errorHandling.o
math_lib.o         : $(OBJ_PATH)/cosp_kinds.o array_lib.o mrgrnk.o
quickbeam.o        : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o          \
                     $(OBJ_PATH)/cosp_stats.o
quickbeam_optics.o : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o          \
                     $(OBJ_PATH)/cosp_constants.o $(OBJ_PATH)/quickbeam.o        \
	             cosp_errorHandling.o array_lib.o math_lib.o optics_lib.o 
mrgrnk.o           : $(OBJ_PATH)/cosp_kinds.o
prec_scops.o	   : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o
cosp_optics.o      : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_constants.o

clean_objs:
	rm -f ../obj/*

clean_all:	
	rm -f $(PROG) $(OBJ_PATH)/*.o $(OBJ_PATH)/*.mod fort.* ../obj/*.o ../obj/*.mod

clean:	
	rm -f $(PROG) fort.* ../obj/*.o ../obj/*.mod *.mod *.o

cosp:
	make -C ../../../src/ install

cosp2_io.o: cosp2_io.f90
	$(F90) $(F90FLAGS) -c -I$(OBJ_PATH) -I$(NCDF_INC) -L${NCDF_LIB} -lnetcdff -lnetcdf -L${LIB} -lnetcdf  $<

mo_rng.o : $(HOME_PATH)/subsample_and_optics/subcol/rng/mo_rng.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

scops.o : $(HOME_PATH)/subsample_and_optics/subcol/scops.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

prec_scops.o : $(HOME_PATH)/subsample_and_optics/subcol/prec_scops.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

cosp_utils.o :  $(HOME_PATH)/subsample_and_optics/optics/cosp_utils.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

cosp_errorHandling.o : $(SRC_F90)/hooks/cosp_errorHandling.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

cosp_optics.o :  $(HOME_PATH)/subsample_and_optics/optics/cosp_optics.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

quickbeam_optics.o : $(HOME_PATH)/subsample_and_optics/optics/quickbeam_optics/quickbeam_optics.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

array_lib.o : $(HOME_PATH)/subsample_and_optics/optics/quickbeam_optics/array_lib.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

math_lib.o : $(HOME_PATH)/subsample_and_optics/optics/quickbeam_optics/math_lib.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

mrgrnk.o : $(HOME_PATH)/subsample_and_optics/optics/quickbeam_optics/mrgrnk.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

optics_lib.o : $(HOME_PATH)/subsample_and_optics/optics/quickbeam_optics/optics_lib.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<	
