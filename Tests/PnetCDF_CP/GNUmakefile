AMREX_HOME ?= ../../

DEBUG	= TRUE
#DEBUG	= FALSE

DIM	= 3

COMP    = gnu

PRECISION = DOUBLE

PNETCDF_HOME   = /home/khl7265/.local
USE_PNETCDF = TRUE
PNC_MF_LAYOUT=2

USE_MPI   = TRUE
#USE_OMP   = TRUE
USE_OMP   = FALSE

TINY_PROFILE = TRUE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules