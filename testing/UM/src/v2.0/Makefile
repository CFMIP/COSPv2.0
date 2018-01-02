.SUFFIXES : .F .f .c .o .a .f90 .f95 .F90
########################################################################
# Makefile for cosp_test_v1.5 driver 
########################################################################
# Intel (ifort)
F90      = ifort
F90FLAGS = -O3 -traceback
F90_LIB  = /usr/local/ifort
F90_LIB   = /scratch/bmcinnes/usr/local/ifort
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
CMOR_INC    = $(F90_LIB)/include
CMOR_LIB    = $(F90_LIB)/lib
UDUNITS_LIB = $(F90_LIB)/lib
UDUNITS_INC = $(F90_LIB)/include
UUID_LIB    = $(F90_LIB)/lib
UUID_INC    = $(F90_LIB)/include

# COSP source directory
HOME_PATH = /home/dswales/Projects/COSP/COSPv2.0
SRC_F90   = ${HOME_PATH}/src
OBJ_PATH  = ${HOME_PATH}/obj
PROG      = cosp_test

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
	cosp_io.o cosp_errorHandling.o mrgrnk.o prec_scops.o array_lib.o     \
	quickbeam_optics.o math_lib.o optics_lib.o

all: $(PROG)
	mkdir obj
	mv *.mod *.o obj/

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(PROG).f90 $(OBJS) \
	-L${CMOR_LIB} -L. -lcmor -I$(CMOR_INC) \
	-I$(NCDF_INC) -L${NCDF_LIB} -lnetcdff -lnetcdf -L${LIB} -lnetcdf \
	-L${RTTOV_LIB_PATH} $(RTTOV_LIBS) \
	-I$(OBJ_PATH) \
	-L${UDUNITS_LIB} -Wl,-rpath=${UDUNITS_LIB} -ludunits2 -lexpat -I${UDUNITS_INC} \
	-L${UUID_LIB} -Wl,-rpath=${UUID_LIB} -luuid -I$(UUID_INC) \
	-o $(PROG)

$(PROG).o          : $(OBJ_PATH)/cosp_config.o $(OBJ_PATH)/cosp.o cosp_io.o   \
                     $(OBJ_PATH)/quickbeam.o
cosp_io.o       : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp.o                 \
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

clean:	
	rm -rf fort.* obj/

clean_cosp:	
	rm -rf $(PROG) $(OBJ_PATH)/*.o $(OBJ_PATH)/*.mod fort.* obj/
	make -C ../../../../src/ clean

cosp:
	make -C ../../../../src/ install

mo_rng.o : subsample_and_optics/subcol/rng/mo_rng.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

scops.o : subsample_and_optics/subcol/scops.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

prec_scops.o : subsample_and_optics/subcol/prec_scops.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

cosp_utils.o :  subsample_and_optics/optics/cosp_utils.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

cosp_io.o: cosp_io.f90
	$(F90) $(F90FLAGS) -c -I$(CMOR_INC) -I$(OBJ_PATH) $<

cosp_errorHandling.o : $(SRC_F90)/hooks/cosp_errorHandling.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

cosp_optics.o :  subsample_and_optics/optics/cosp_optics.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

quickbeam_optics.o : subsample_and_optics/optics/quickbeam_optics/quickbeam_optics.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

array_lib.o : subsample_and_optics/optics/quickbeam_optics/array_lib.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

math_lib.o : subsample_and_optics/optics/quickbeam_optics/math_lib.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

mrgrnk.o : subsample_and_optics/optics/quickbeam_optics/mrgrnk.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<

optics_lib.o : subsample_and_optics/optics/quickbeam_optics/optics_lib.F90
	$(F90) $(F90FLAGS) -c -I $(OBJ_PATH) $<	
