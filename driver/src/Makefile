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
F90FLAGS = -module $(mod) -fpp -O0 -g -traceback

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
HDF5_LIB    = /sw/rhel6-x64/hdf5/hdf5-1.8.14-intel14/lib


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
RTTOVSIM_PATH = $(SIM_PATH)/rttov

HOOKS_PATH = $(cospdir)/hooks
RNG_PATH = $(cospdir)/rng

# RTTOV setup
RTTOV_PATH       = $(HOME)/work/Tools/RTTOV/rttov113
RTTOV_LIB_PATH   = $(RTTOV_PATH)/lib 
RTTOV_INC_PATH   = $(RTTOV_PATH)/include 
RTTOV_MOD_PATH   = $(RTTOV_PATH)/mod
RTTOV_LIBS       =  -lrttov11_wrapper -lrttov11_other -lrttov11_emis_atlas -lrttov11_brdf_atlas -lrttov11_mw_scatt -lrttov11_mw_scatt_coef -lrttov11_parallel -lrttov11_coef_io -lrttov11_hdf -lrttov11_main

## Other
PROG     = cosp_test_v2.0
COSP_INTERFACE_F90 = $(HOME)/github/Tools/COSPv2.0/cosp-interface


# List of library files that will be created
# ---------------------------------------------------------------------------------------------------------------------------------------
LIB_COSP = $(lib)/libcosp.a
LIB_COSPTEST = $(lib)/librun.a

FLAG_LIB_COSP = -L$(lib_cosp) -lcosp -L$(lib) -lrun
# ---------------------------------------------------------------------------------------------------------------------------------------

# Final flag indicating the library and include paths
# ---------------------------------------------------------------------------------------------------------------------------------------
FLAG_CMOR= -L${CMOR_LIB} -lcmor -I$(CMOR_INC)
FLAG_NCDF= -I$(NCDF_INC) -L${NCDF_LIB} -lnetcdff -L${NCDF_C_LIB} -lnetcdf -Wl,-rpath,${NCDF_LIB} -Wl,-rpath,${NCDF_C_LIB}
FLAG_RTTOV= -I${RTTOV_INC_PATH} -L${RTTOV_LIB_PATH} $(RTTOV_LIBS) 
FLAG_UDUNITS= -L${UDUNITS_LIB} -ludunits2 -lexpat -Wl,-rpath,${UDUNITS_LIB} -I${UDUNITS_INC}
FLAG_UUID= -L${UUID_LIB} -luuid -Wl,-rpath,${UUID_LIB} -I$(UUID_INC)
FLAG_HDF5= -L${HDF5_LIB} -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm -Wl,-rpath,${HDF5_LIB}

FLAGS_COSPTEST = $(FLAG_LIB_COSP) $(FLAG_CMOR) $(FLAG_RTTOV) $(FLAG_HDF5) $(FLAG_NCDF) $(FLAG_UDUNITS) $(FLAG_UUID)
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
		$(obj)/cosp_rttov11.o \
		$(obj)/MISR_simulator.o \
		$(obj)/modis_simulator.o \
		$(obj)/lidar_simulator.o \
		$(obj)/cosp_calipso_interface.o \
		$(obj)/cosp_cloudsat_interface.o \
		$(obj)/cosp_isccp_interface.o \
		$(obj)/cosp_misr_interface.o \
		$(obj)/cosp_parasol_interface.o \
		$(obj)/cosp_modis_interface.o \
		$(obj)/cosp_rttov_interface11.o \
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


clean_objs:
	rm -f $(OBJ_PATH) *.mod *.o

clean:	
	rm -f $(obj)/*.o $(mod)/*.mod $(lib)/*.a



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

$(obj)/cosp_rttov_interface11.o : $(SIM_PATH)/cosp_rttov_interface.F90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

$(obj)/MISR_simulator.o : $(MISR_PATH)/MISR_simulator.F90
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/modis_simulator.o : $(MODIS_PATH)/modis_simulator.F90 
	$(F90) $(F90FLAGS) -c $< -o $@

$(obj)/cosp_rttov11.o : $(RTTOVSIM_PATH)/cosp_rttov11.F90
	$(F90) $(F90FLAGS) -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) -L $(RTTOV_LIB_PATH) -c $< -o $@

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
