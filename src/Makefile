# File name: Makefile          
# 
# Copyright (C) 2009 Thomas Reerink & Michael Kliphuis. This program
# is distributed under the terms of the GNU General Public License.
#
# This file is part of OBLIMAP.
#
# OBLIMAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OBLIMAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License    
# along with OBLIMAP.  If not, see <http://www.gnu.org/licenses/>.
# 
#
# OBLIMAP is maintained by:
#
# Thomas Reerink
# Institute for Marine and Atmospheric Research Utrecht (IMAU)
# Utrecht University
# Princetonplein 5
# 3584 CC Utrecht
# The Netherlands
#
# email: <t.reerink@uu.nl>
#
# OBLIMAP is hosted on NeSCForge: NOT YET
#
# http://forge.nesc.ac.uk/projects/oblimap/
#

#    Makefile for compiling OBLIMAP
#
#    To execute the make file (in the src/ directory) e.g.:
#     make oblimap_gcm_to_im_program                        
#     make oblimap_im_to_gcm_program                        

include Makefile.gfortran
#include Makefile.pgf90
#include Makefile.g95
#include Makefile.mpfort

# Since on some platforms the compiler seems not to be able to put the object 
# files and modules in any directory other than the current one. If the OBJ_PATH 
# and MOD_PATH are not specified in the file included above, then we can set it 
# to any value we like at this point:
OBJ_PATH ?= obj
MOD_PATH ?= mod
        
# Pattern rules for creating object files from fortran 90 source files
$(OBJ_PATH)/%.o : %.f90
	mkdir -p $(OBJ_PATH) $(MOD_PATH) && $(F90) $(F90FLAGS) -c -o $@  $<

# This rule is used to build an executable from object files: 
% : 
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

    

# =============================================================
#	BUILDING THE OBLIMAP_GCM_TO_IM_PROGRAM
# =============================================================

# Source files used for the oblimap_gcm_to_im_program (this program only uses fortran 90 source files):
OBLIMAP_GCM_TO_IM_PROGRAM_SOURCES = \
        oblimap_configuration_module.f90 \
        oblimap_mapping_module.f90 \
        oblimap_read_and_write_module.f90 \
        oblimap_projection_module.f90 \
        oblimap_scan_contributions_module.f90 \
        oblimap_gcm_to_im_program.f90
        
# All object files used for the oblimap_gcm_to_im_program (this program only uses fortran 90 source files):
OBLIMAP_GCM_TO_IM_PROGRAM_OBJECTS = $(OBLIMAP_GCM_TO_IM_PROGRAM_SOURCES:%.f90=$(OBJ_PATH)/%.o)

# Update lists of all sources. This is used later on to create a dependency file automatically. See the "deps" target.
ALL_SOURCES += $(OBLIMAP_GCM_TO_IM_PROGRAM_SOURCES)

# Update a list of all programs, which is used by the "all" target
ALL_PROGRAMS += oblimap_gcm_to_im_program

oblimap_gcm_to_im_program : $(OBLIMAP_GCM_TO_IM_PROGRAM_OBJECTS) 



# =============================================================
#	BUILDING THE OBLIMAP_IM_TO_GCM_PROGRAM
# =============================================================

# Source files used for the oblimap_im_to_gcm_program (this program only uses fortran 90 source files):
OBLIMAP_IM_TO_GCM_PROGRAM_SOURCES = \
        oblimap_configuration_module.f90 \
        oblimap_mapping_module.f90 \
        oblimap_read_and_write_module.f90 \
        oblimap_projection_module.f90 \
        oblimap_scan_contributions_module.f90 \
        oblimap_post_processing_module.f90 \
        oblimap_im_to_gcm_program.f90
        
# All object files used for the oblimap_im_to_gcm_program (this program only uses fortran 90 source files):
OBLIMAP_IM_TO_GCM_PROGRAM_OBJECTS = $(OBLIMAP_IM_TO_GCM_PROGRAM_SOURCES:%.f90=$(OBJ_PATH)/%.o)

# Update lists of all sources. This is used later on to create a dependency file automatically. See the "deps" target.
ALL_SOURCES += $(OBLIMAP_IM_TO_GCM_PROGRAM_SOURCES)

# Update a list of all programs, which is used by the "all" target
ALL_PROGRAMS += oblimap_im_to_gcm_program

oblimap_im_to_gcm_program : $(OBLIMAP_IM_TO_GCM_PROGRAM_OBJECTS) 



# =============================================================
#	AUTOMATIC DEPENDENCIES
# =============================================================

# The dependencies are listed in the file 'deps'. The following rule updates the dependencies
# automatically with the makedepf90 program. If makedepf90 is not available comment the next
# two lines, then the already existing deps will be used. (Take care of a real tab).
##deps : $(ALL_SOURCES)
##	makedepf90 -b $(OBJ_PATH) $^ > deps

include deps



# =======================================
#       CLEANING UP EVERYTHING
# =======================================

.PHONY : clean

clean : 
	rm -fv $(ALL_PROGRAMS) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod $(MOD_PATH)/*.MOD *.mod *.lst

        

# ==================================
# 	BUILD EVERYTHING
# ==================================

.PHONY : all

all : $(ALL_PROGRAMS)
