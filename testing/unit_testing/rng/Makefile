.SUFFIXES : .f90
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Copyright (c) 2016, Regents of the University of Colorado
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are 
# permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of 
#    conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list
#    of conditions and the following disclaimer in the documentation and/or other 
#    materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be 
#    used to endorse or promote products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
# THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# History
# Jan 2016 - D. Swales - Original version
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compiler and flags
F90      = ifort
F90FLAGS = # -O2 -g -C

# Source file location
SRC  = /home/dswales/Projects/COSP/COSPv2.0

######################################################################## 
PROG = test_rng
OBJS = cosp_kinds.o mo_rng.o

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(PROG).f90 $(OBJS) \
	-o $(PROG)

%.o : %.f90
	${F90} ${F90FLAGS} -c $<

$(PROG).o : cosp_kinds.o mo_rng.o
mo_rng.o  : cosp_kinds.o

clean:
	rm -f *.o test_rng *.mod

cosp_kinds.o : $(SRC)/src/hooks/cosp_kinds.f90
	$(F90) $(F90FLAGS) -c $<

mo_rng.o : $(SRC)/src/rng/mo_rng.f90
	$(F90) $(F90FLAGS) -c $<










