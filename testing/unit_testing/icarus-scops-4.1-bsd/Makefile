.SUFFIXES : .f90
# *****************************COPYRIGHT****************************
# (c) British Crown Copyright 2009, the Met Office.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the
# following conditions are met:
# 
#     * Redistributions of source code must retain the above 
#       copyright  notice, this list of conditions and the following 
#       disclaimer.
#     * Redistributions in binary form must reproduce the above 
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials 
#       provided with the distribution.
#     * Neither the name of the Met Office nor the names of its 
#       contributors may be used to endorse or promote products
#       derived from this software without specific prior written 
#       permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
# 
# *****************************COPYRIGHT*******************************
# *****************************COPYRIGHT*******************************
# *****************************COPYRIGHT*******************************

# Compiler and flags
F90        = ifort
DEBUG      = #-g #-C
F90FLAGS   = -O2 -fpp

# Source file location
ROOT       = /home/dswales/Projects/COSP/COSPv2.0
SRC        = $(ROOT)/src/simulator/icarus-scops-4.1-bsd

######################################################################## 
PROG1 = test_isccp_cloud_types
OBJS1 = cosp_kinds.o cosp_constants.o isccp_cloud_types.o scops.o \
        icarus.o cosp_optics.o cosp_stats.o cosp_config.o mo_rng.o\
	cosp_errorHandling.o

all: $(PROG1)

$(PROG1): $(OBJS1) 
	$(F90) $(F90FLAGS) $(PROG1).f90 $(OBJS1) \
	-o $(PROG1)

%.o : %.f90
	${F90} ${F90FLAGS} -c $<

$(PROG1).o               : cosp_kinds.o
cosp_errorHandling.o     : cosp_kinds.o
isccp_cloud_types.o      : cosp_kinds.o icarus.o mo_rng.o
scops.o                  : cosp_kinds.o mo_rng.o cosp_errorHandling.o
icarus.o	         : cosp_kinds.o cosp_constants.o cosp_optics.o \
			   cosp_stats.o scops.o cosp_config.o
mo_rng.o                 : cosp_kinds.o
cosp_stats.o             : cosp_kinds.o cosp_config.o

clean:
	rm -f *.o test_isccp_cloud_types ftn* test_congvec *.mod
test:
	./test_isccp_cloud_types.ksh

cosp_kinds.o : $(ROOT)/src/hooks/cosp_kinds.f90
	$(F90) $(F90FLAGS) -c $<

cosp_errorHandling.o : $(ROOT)/src/hooks/cosp_errorHandling.f90
	$(F90) $(F90FLAGS) -c $<

cosp_constants.o :  $(ROOT)/src/hooks/cosp_constants.f90
	$(F90) $(F90FLAGS) -c $<

scops.o :  $(ROOT)/src/subcol/scops.f90
	$(F90) $(F90FLAGS) -c $<

icarus.o :  $(SRC)/icarus.f90
	$(F90) $(F90FLAGS) -c $<

cosp_optics.o :  $(ROOT)/src/optics/cosp_optics.f90
	$(F90) $(F90FLAGS) -c $<

cosp_config.o :  $(ROOT)/src/cosp_config.f90
	$(F90) $(F90FLAGS) -c $<

cosp_stats.o :  $(ROOT)/src/cosp_stats.f90
	$(F90) $(F90FLAGS) -c $<

mo_rng.o : $(ROOT)/src/rng/mo_rng.f90
	$(F90) $(F90FLAGS) -c $<






