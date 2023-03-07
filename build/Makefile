# Paths to source codes in COSP library
HOOKS        = ../model-interface
SRC_PATH     = ../src
SIM_PATH     = $(SRC_PATH)/simulator
RT_PATH      = $(SIM_PATH)/rttov
RS_PATH      = $(SIM_PATH)/quickbeam
CS_PATH      = $(SIM_PATH)/actsim
ISCCP_PATH   = $(SIM_PATH)/icarus
MISR_PATH    = $(SIM_PATH)/MISR_simulator
MODIS_PATH   = $(SIM_PATH)/MODIS_simulator
PARASOL_PATH = $(SIM_PATH)/parasol

# RTTOV variables. You may need to change this. JKS
# RTTOV_PATH     = /glade/u/home/jonahshaw/w/RTTOV
# RTTOV_LIB_PATH = $(RTTOV_PATH)/lib 
# RTTOV_INC_PATH = $(RTTOV_PATH)/include 
# RTTOV_MOD_PATH = $(RTTOV_PATH)/mod

# Makefile.rttov needs to define: RTTOVDIR, FFLAGS, LIBDIR, LIBS, LLIBS
include Makefile.rttov
# RTTOV_LIB_PATH = $(RTTOVDIR)/lib  # JKS - this became LIBDIR in Makefile.rttov
# RTTOV_INC_PATH = $(RTTOVDIR)/include # JKS - this gets folded into FFLAGS
# RTTOV_MOD_PATH = $(RTTOVDIR)/mod     # JKS - this gets folded into FFLAGS



VPATH = $(SRC_PATH):$(HOOKS):$(SIM_PATH):$(RT_PATH):$(RS_PATH):$(CS_PATH):$(ISCCP_PATH):$(MISR_PATH):$(MODIS_PATH):$(PARASOL_PATH)

# Example subcolumn generaton and mapping to optical properties, following COSP 1.4
SUBCOLS = ../subsample_and_optics_example
VPATH += $(SUBCOLS)/optics:$(SUBCOLS)/optics/quickbeam_optics:$(SUBCOLS)/subcol:$(SUBCOLS)/subcol/rng

# Standalone driver
DRIVER = ../driver/src
VPATH += $(DRIVER)

# Makefile.conf needs to define $(F90), $(F90FLAGS), $(NC_LIB), $(NC_INC)
include Makefile.conf
F90FLAGS += $(NC_INC)
LDFLAGS  += $(NC_LIB) -lnetcdff

LDFLAGS += -L. -lcosp -lsubcol

ifdef RTTOV
# LDFLAGS += -L${RTTOV_LIB_PATH} -lrttov13
#LDFLAGS += -L${RTTOV_LIB_PATH} -lrttov_main
#LDFLAGS += -L${LIBDIR} -lrttov_main
#LDFLAGS += $(LIBS) $(LLIBS)
LDFLAGS += $(LIBS) $(LLIBS) $(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)
endif


# Rules
%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%: %.o
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)
#	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# Special rule for the RTTOV file. Had an ifdef line but it wasn't necessary.

# JKS is this needed because of the load_coefficient calls? Can I remove this by restructuring?
cosp_rttov_interface_v13.o : cosp_rttov_interface_v13.F90 
	$(F90) $(F90FLAGS) $(FFLAGS) -c $<
#	$(F90) $(F90FLAGS) -c -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) $<

cosp_rttov_v13.o : cosp_rttov_v13.F90 
	$(F90) $(F90FLAGS) $(FFLAGS) -c $<
#	$(F90) $(F90FLAGS) -c -I $(RTTOV_INC_PATH) -I $(RTTOV_MOD_PATH) $<

#
# The COSP library
#
COSP_OBJS = cosp.o cosp_config.o cosp_stats.o cosp_constants.o cosp_errorHandling.o cosp_kinds.o \
            cosp_isccp_interface.o icarus.o cosp_misr_interface.o MISR_simulator.o               \
            cosp_modis_interface.o  modis_simulator.o cosp_parasol_interface.o  parasol.o        \
            cosp_calipso_interface.o  lidar_simulator.o cosp_cloudsat_interface.o quickbeam.o    \
            cosp_grLidar532_interface.o cosp_atlid_interface.o

# Add RTTOV files appropriately. JKS: cosp should only need to depend on the interface, not the simulator.
ifdef RTTOV
COSP_OBJS += cosp_rttov_interface_v13.o cosp_rttov_v13.o
else
COSP_OBJS += cosp_rttov_interfaceSTUB.o cosp_rttovSTUB.o
endif


LIBRARY = libcosp.a

$(LIBRARY): $(COSP_OBJS)
	ar -rvs $(LIBRARY) $(COSP_OBJS)

#
# Alternative target for library
#
lib: $(LIBRARY)

clean:
	rm -f *.optrpt *.mod *.o *.a fort.* cosp2_test

driver:
	make cosp2_test
	mv cosp2_test ../driver/run

# I would like to set the env variable here, but it won't work.
driver_rttov:
	make driver


# Dependencies for normal COSP libary
cosp.o : cosp_kinds.o cosp_modis_interface.o cosp_constants.o                              \
         cosp_misr_interface.o cosp_isccp_interface.o cosp_calipso_interface.o             \
         cosp_cloudsat_interface.o  cosp_stats.o                                           \
         cosp_parasol_interface.o                                                          \
         quickbeam.o MISR_simulator.o lidar_simulator.o                                    \
         parasol.o icarus.o cosp_grLidar532_interface.o cosp_atlid_interface.o

# Conditionally add dependencies on the Stub or actual RTTOV simulator.
ifdef RTTOV
cosp.o: cosp_rttov_interface_v13.o cosp_rttov_v13.o
else
cosp.o: cosp_rttov_interfaceSTUB.o cosp_rttovSTUB.o
endif


cosp_constants.o            : cosp_kinds.o
cosp_config.o               : cosp_kinds.o
cosp_errorHandling.o        : cosp_kinds.o
cosp_stats.o                : cosp_kinds.o cosp_config.o cosp_constants.o
cosp_isccp_interface.o      : cosp_kinds.o icarus.o
icarus.o                    : cosp_kinds.o cosp_constants.o  cosp_stats.o
cosp_misr_interface.o       : cosp_kinds.o
MISR_simulator.o            : cosp_kinds.o cosp_config.o cosp_stats.o
cosp_modis_interface.o      : cosp_kinds.o cosp_config.o modis_simulator.o
modis_simulator.o           : cosp_kinds.o cosp_config.o cosp_stats.o
cosp_parasol_interface.o    : cosp_kinds.o
parasol.o                   : cosp_kinds.o cosp_config.o cosp_constants.o
cosp_calipso_interface.o    : cosp_kinds.o lidar_simulator.o
lidar_simulator.o           : cosp_kinds.o cosp_config.o cosp_stats.o
cosp_grLidar532_interface.o : cosp_kinds.o
cosp_atlid_interface.o      : cosp_kinds.o
cosp_cloudsat_interface.o   : cosp_kinds.o cosp_config.o quickbeam.o
quickbeam.o                 : cosp_kinds.o cosp_stats.o

# RTTOV Simulator code dependencies (only one interface and one RTTOV simulator are used at a time)
cosp_rttov_interface_v13.o  : cosp_kinds.o cosp_config.o cosp_rttov_v13.o
cosp_rttov_v13.o            : cosp_kinds.o cosp_config.o cosp_constants.o
cosp_rttov_interfaceSTUB.o  : cosp_kinds.o cosp_config.o cosp_rttovSTUB.o
cosp_rttovSTUB.o            : cosp_kinds.o cosp_config.o cosp_constants.o


# Example subcolumn generaton and mapping to optical properties, following COSP 1.4
SUBCOL_OBJS = mo_rng.o scops.o prec_scops.o cosp_utils.o cosp_optics.o quickbeam_optics.o array_lib.o math_lib.o mrgrnk.o optics_lib.o cosp_errorHandling.o
libsubcol.a: $(SUBCOL_OBJS) libcosp.a
	ar -rvs libsubcol.a $(SUBCOL_OBJS)


array_lib.o       : cosp_kinds.o
mrgrnk.o          : cosp_kinds.o
math_lib.o        : cosp_kinds.o cosp_errorHandling.o
optics_lib.o      : cosp_kinds.o cosp_errorHandling.o
quickbeam_optics.o: cosp_kinds.o cosp_errorHandling.o cosp_constants.o cosp_config.o mrgrnk.o array_lib.o optics_lib.o math_lib.o quickbeam.o
scops.o      : cosp_kinds.o cosp_errorHandling.o mo_rng.o
prec_scops.o : cosp_kinds.o cosp_config.o
cosp_utils.o : cosp_kinds.o cosp_config.o
cosp_optics.o: cosp_kinds.o cosp_config.o cosp_constants.o modis_simulator.o
mo_rng.o     : cosp_kinds.o

# Offline driver
cosp2_io.o    : cosp_kinds.o cosp_config.o cosp.o
cosp2_test.o  : cosp2_io.o libcosp.a libsubcol.a
cosp2_test    : cosp2_test.o cosp2_io.o


###############
### OLD
###############

# Offline driver with RTTOV is just adding the COSP-RTTOV file and the RTTOV library.
#cosp2_rttov_test : cosp2_test.o cosp2_io.o cosp_rttov_v13.o


# OLD Offline driver with RTTOV
#cosp2_rttov_io.o    : cosp_kinds.o cosp_config.o cosp_rttov.o
#cosp2_rttov_test.o  : cosp2_rttov_io.o libcosp_rttov.a libsubcol_rttov.a
#cosp2_rttov_test    : cosp2_rttov_test.o cosp2_rttov_io.o


# Dependencies for COSP-RTTOV libary
#cosp_rttov.o : cosp_kinds.o cosp_modis_interface.o cosp_constants.o \
#	cosp_misr_interface.o cosp_isccp_interface.o cosp_calipso_interface.o             \
#	cosp_cloudsat_interface.o  cosp_stats.o                                           \
#	cosp_parasol_interface.o cosp_rttov_v13.o \
#        cosp_rttov_interface_v13.o quickbeam.o MISR_simulator.o lidar_simulator.o         \
#        parasol.o icarus.o cosp_grLidar532_interface.o cosp_atlid_interface.o


### Library for COSP-RTTOV build
#COSP_RTTOV_OBJS = cosp.o cosp_config.o cosp_stats.o cosp_constants.o cosp_errorHandling.o cosp_kinds.o \
#                  cosp_isccp_interface.o icarus.o cosp_misr_interface.o MISR_simulator.o               \
#                  cosp_modis_interface.o  modis_simulator.o cosp_parasol_interface.o  parasol.o        \
#                  cosp_calipso_interface.o  lidar_simulator.o cosp_cloudsat_interface.o quickbeam.o    \
#                  cosp_rttov_interface_v13.o cosp_rttov_v13.o cosp_grLidar532_interface.o cosp_atlid_interface.o

#LIBRARY_RTTOV = libcosp_rttov.a

#$(LIBRARY_RTTOV): $(COSP_RTTOV_OBJS)
#	ar -rvs $(LIBRARY_RTTOV) $(COSP_RTTOV_OBJS)


# Is this right?
#libsubcol_rttov.a: $(SUBCOL_OBJS) libcosp_rttov.a
#	ar -rvs libsubcol_rttov.a $(SUBCOL_OBJS)

# Will need to be updated. Conditionally removed rttov_interface dependence on cosp_rttov_v13
#cosp_rttov_interface_v13.o  : cosp_kinds.o cosp_config.o

#ifdef RTTOV
#cosp_rttov_interface_v13.o  : cosp_kinds.o cosp_config.o cosp_rttov_v13.o
#else
#cosp_rttov_interface_v13.o  : cosp_kinds.o cosp_config.o cosp_rttovSTUB.o
#endif

#cosp_rttov_interface_v13.o  : cosp_kinds.o cosp_config.o cosp_rttov_v13.o
#cosp_rttov_v13.o            : cosp_kinds.o cosp_config.o cosp_constants.o

#ifdef RTTOV
#COSP_OBJS = cosp.o cosp_config.o cosp_stats.o cosp_constants.o cosp_errorHandling.o cosp_kinds.o \
#            cosp_isccp_interface.o icarus.o cosp_misr_interface.o MISR_simulator.o               \
#            cosp_modis_interface.o  modis_simulator.o cosp_parasol_interface.o  parasol.o        \
#            cosp_calipso_interface.o  lidar_simulator.o cosp_cloudsat_interface.o quickbeam.o    \
#            cosp_grLidar532_interface.o cosp_atlid_interface.o                                   \
#            cosp_rttov_interface_v13.o cosp_rttov_v13.o
#else
#COSP_OBJS = cosp.o cosp_config.o cosp_stats.o cosp_constants.o cosp_errorHandling.o cosp_kinds.o \
#            cosp_isccp_interface.o icarus.o cosp_misr_interface.o MISR_simulator.o               \
#            cosp_modis_interface.o  modis_simulator.o cosp_parasol_interface.o  parasol.o        \
#            cosp_calipso_interface.o  lidar_simulator.o cosp_cloudsat_interface.o quickbeam.o    \
#            cosp_grLidar532_interface.o cosp_atlid_interface.o                                   \
#            cosp_rttov_interfaceSTUB.o cosp_rttovSTUB.o
#endif