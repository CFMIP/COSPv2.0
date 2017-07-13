.SUFFIXES : .F .f .c .o .a .f90 .f95 .F90
########################################################################
# Makefile for cosp_test_v1.5 driver 
########################################################################

basedir=$(HOME)/github/Tools/COSPv2.0/driver/src
cospdir=$(HOME)/github/Tools/COSPv2.0/src

obj  = $(basedir)/obj
lib = $(basedir)/lib
mod = $(basedir)/mod

obj_cosp  = $(basedir)/obj
lib_cosp = $(basedir)/lib
mod_cosp = $(basedir)/mod

# Intel (ifort)
F90      = ifort
F90FLAGS = -O3 -module $(mod) -fpp ##-g -traceback
#F90FLAGS = -g -debug all -fpe0 -traceback -check bounds -check all -O0 
##F90_LIB  = /sw/rhel6-x64/intel/intel-17.0.2
# Portland Group (pgf90)
#F90      = pgf90
#F90FLAGS = -Mpreprocess -O0 -g
#F90_LIB  = /usr/local/pgf90
# Gfortran (gfortran)
# F90      = gfortran
##F90FLAGS = -cpp -O0 -fbacktrace -fbounds-chec
# F90FLAGS = -O0 -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -g

INC = /usr/local/include

NCDF_C_LIB = /sw/rhel6-x64/netcdf/netcdf_c-4.3.2-intel14/lib
NCDF_INC    = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/include
NCDF_LIB    = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/lib
CMOR_INC    = /sw/rhel6-x64/cmor-2.9.2-intel14/include
CMOR_LIB    = /sw/rhel6-x64/cmor-2.9.2-intel14/lib
UDUNITS_LIB = /sw/rhel6-x64/util/udunits-2.2.17-intel15/lib
UDUNITS_INC = /sw/rhel6-x64/util/udunits-2.2.17-intel15/include
UUID_LIB    = /pf/b/b380333/usr/uuid-1.6.2/lib
UUID_INC    = /pf/b/b380333/usr/uuid-1.6.2/include

# NCDF_C_LIB = /sw/rhel6-x64/netcdf/netcdf_c-4.4.0-gcc48/lib
# NCDF_INC    = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc48/include
# NCDF_LIB    = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc48/lib
# CMOR_INC    = /sw/rhel6-x64/cmor-2.9.2-gcc48/include
# CMOR_LIB    = /sw/rhel6-x64/cmor-2.9.2-gcc48/lib
# UDUNITS_LIB = /sw/rhel6-x64/util/udunits-2.2.17-gcc48/lib
# UDUNITS_INC = /sw/rhel6-x64/util/udunits-2.2.17-gcc48/include
# UUID_LIB    = /pf/b/b380333/usr/uuid-1.6.2/lib
# UUID_INC    = /pf/b/b380333/usr/uuid-1.6.2/include

# Simulator directories
SIM_PATH     = $(cospdir)/simulator
RT_PATH      = $(SIM_PATH)/rttov
RS_PATH      = $(SIM_PATH)/quickbeam
CS_PATH      = $(SIM_PATH)/actsim
LLNL_PATH    = $(SIM_PATH)/llnl
ISCCP_PATH   = $(SIM_PATH)/icarus-scops-4.1-bsd
MISR_PATH    = $(SIM_PATH)/MISR_simulator
MODIS_PATH   = $(SIM_PATH)/MODIS_simulator
PARASOL_PATH = $(SIM_PATH)/parasol
RTTOVSTUB_PATH = $(SIM_PATH)/rttov

HOOKS_PATH = $(cospdir)/hooks
RNG_PATH = $(cospdir)/rng

# RTTOV setup
RTTOV_PATH       = $(HOME)/work/Tools/RTTOV/rttov113
RTTOV_LIB_PATH   = $(RTTOV_PATH)/lib 
RTTOV_INC_PATH   = $(RTTOV_PATH)/include 
RTTOV_MOD_PATH   = $(RTTOV_PATH)/mod
RTTOV_LIBS       = -lrttov11_mw_scatt -lrttov11_brdf_atlas -lrttov11_emis_atlas \
                   -lrttov11_other -lrttov11_parallel -lrttov11_coef_io        \
                   -lrttov11_main

## Other
PROG     = cosp_test_v2.0
COSP_INTERFACE_F90 = $(HOME)/github/Tools/COSPv2.0/cosp-interface


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
LIB_COSPTEST = $(lib)/librun.a

FLAG_LIB_COSP = -L$(lib_cosp) -lcosp -L$(lib) -lrun
# ---------------------------------------------------------------------------------------------------------------------------------------

# Final flag indicating the library and include paths
# ---------------------------------------------------------------------------------------------------------------------------------------
FLAG_CMOR= -L${CMOR_LIB} -lcmor -I$(CMOR_INC)
FLAG_NCDF= -I$(NCDF_INC) -L${NCDF_LIB} -lnetcdff -L${NCDF_C_LIB} -lnetcdf
FLAG_RTTOV= -L${RTTOV_LIB_PATH} $(RTTOV_LIBS) -Wl,-rpath,${NCDF_LIB} -Wl,-rpath,${NCDF_C_LIB}
FLAG_UDUNITS= -L${UDUNITS_LIB} -Wl,-rpath=${UDUNITS_LIB} -ludunits2 -lexpat -I${UDUNITS_INC}
FLAG_UUID= -L${UUID_LIB} -Wl,-rpath=${UUID_LIB} -luuid -I$(UUID_INC)

FLAGS_COSPTEST = $(FLAG_LIB_COSP) $(FLAG_CMOR) $(FLAG_NCDF) $(FLAG_RTTOV) $(FLAG_UDUNITS) $(FLAG_UUID)
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

LIST_OBJ_COSPTEST =     $(obj)/mo_rng.o \
			$(obj)/cosp_errorHandling.o \
			$(obj)/mrgrnk.o \
			$(obj)/scops.o \
			$(obj)/array_lib.o \
			$(obj)/math_lib.o \
			$(obj)/optics_lib.o \
			$(obj)/cosp_utils.o \
			$(obj)/cosp_optics.o \
			$(obj)/cosp_io.o \
			$(obj)/prec_scops.o \
			$(obj)/quickbeam_optics.o
# ---------------------------------------------------------------------------------------------------------------------------------------

install : $(LIST_OBJ_COSP) $(LIST_OBJ_COSPTEST)
	ar r $(LIB_COSP) $(LIST_OBJ_COSP)
	ar r $(LIB_COSPTEST) $(LIST_OBJ_COSPTEST)
	$(F90) $(F90FLAGS) $(PROG).f90 -o $(PROG) $(FLAGS_COSPTEST) 


# all: $(PROG)
# 	mv *.mod *.o obj/

# $(PROG): $(OBJS)
# 	$(F90) $(F90FLAGS) $(PROG).f90 $(OBJS) \
# 	-L${CMOR_LIB} -L. -lcmor -I$(CMOR_INC) \
# 	-I$(NCDF_INC) -L${NCDF_LIB} -lnetcdff -L${NCDF_C_LIB} -lnetcdf \
# 	-L${RTTOV_LIB_PATH} $(RTTOV_LIBS) \
# 	-I$(OBJ_PATH) \
# 	-L${UDUNITS_LIB} -Wl,-rpath=${UDUNITS_LIB} -ludunits2 -lexpat -I${UDUNITS_INC} \
# 	-L${UUID_LIB} -Wl,-rpath=${UUID_LIB} -luuid -I$(UUID_INC) \
# 	-Wl,-rpath,${NCDF_LIB} -Wl,-rpath,${NCDF_C_LIB} \
# 	-o $(PROG)

# $(PROG).o          : $(OBJ_PATH)/cosp_config.o $(OBJ_PATH)/cosp.o cosp_io.v2.o   \
#                      $(OBJ_PATH)/quickbeam.o
# cosp_io.v2.o       : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp.o                 \
#                      $(OBJ_PATH)/cosp_config.o
# cosp_utils.o       : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o
# scops.o            : $(OBJ_PATH)/cosp_kinds.o mo_rng.o cosp_errorHandling.o
# mo_rng.o           : $(OBJ_PATH)/cosp_kinds.o
# optics_lib.o	   : $(OBJ_PATH)/cosp_kinds.o cosp_errorHandling.o
# array_lib.o   	   : $(OBJ_PATH)/cosp_kinds.o cosp_errorHandling.o
# math_lib.o         : $(OBJ_PATH)/cosp_kinds.o array_lib.o mrgrnk.o
# quickbeam.o        : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o          \
#                      $(OBJ_PATH)/cosp_stats.o
# quickbeam_optics.o : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o          \
#                      $(OBJ_PATH)/cosp_constants.o $(OBJ_PATH)/quickbeam.o        \
# 	             cosp_errorHandling.o array_lib.o math_lib.o optics_lib.o 
# mrgrnk.o           : $(OBJ_PATH)/cosp_kinds.o
# prec_scops.o	   : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_config.o
# cosp_optics.o      : $(OBJ_PATH)/cosp_kinds.o $(OBJ_PATH)/cosp_constants.o

clean_objs:
	rm -f $(OBJ_PATH) *.mod *.o

cleanall:	
	rm -f $(PROG) $(OBJ_PATH)/*.o $(OBJ_PATH)/*.mod *.o *.mod fort.*

clean:	
	rm -f $(obj)/*.o $(mod)/*.mod $(lib)/*.a

cosp:
	make -C ../../src/ install


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

$(obj)/cosp_constants.o :  $(HOOKS_PATH)/cosp_constants.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_kinds.o : $(HOOKS_PATH)/cosp_kinds.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_config.o :  $(cospdir)/cosp_config.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_stats.o :  $(cospdir)/cosp_stats.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp.o :  $(cospdir)/cosp.F90
	$(F90) $(F90FLAGS) -c $< -o $@


$(obj)/mo_rng.o : $(RNG_PATH)/mo_rng.F90
	$(F90) $(F90FLAGS)  -c $<  -o $@

$(obj)/scops.o : $(cospdir)/subcol/scops.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/prec_scops.o : $(cospdir)/subcol/prec_scops.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/cosp_utils.o :  $(cospdir)/optics/cosp_utils.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/cosp_io.o: cosp_io.v2.f90
	$(F90) $(F90FLAGS) -I $(CMOR_INC) -I $(NCDF_INC)  -c $< -o $@

$(obj)/cosp_errorHandling.o : $(COSP_INTERFACE_F90)/cosp_errorHandling.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/quickbeam_optics.o :  $(cospdir)/optics/quickbeam_optics.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/cosp_optics.o :  $(cospdir)/optics/cosp_optics.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/array_lib.o : $(cospdir)/optics/array_lib.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/math_lib.o : $(cospdir)/optics/math_lib.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/mrgrnk.o : $(cospdir)/optics/mrgrnk.F90
	$(F90) $(F90FLAGS)  -c $< -o $@

$(obj)/optics_lib.o : $(cospdir)/optics/optics_lib.F90
	$(F90) $(F90FLAGS)  -c $< -o $@	
