# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2019 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for GNU Fortran compiler on Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# NF_CONFIG      NetCDF Fortran configuration script
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF library directory
# NETCDF_LIBS    NetCDF library switches
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := gfortran
           FFLAGS := -frepack-arrays
       FIXEDFLAGS := -132
        FREEFLAGS := -free
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional -w              # -w turns of warnings
           INCDIR := /usr/include /usr/local/bin
            SLIBS := -L/usr/local/lib -L/usr/lib
            ULIBS :=
             LIBS :=
       MOD_SUFFIX := mod
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := -r
            MKDIR := mkdir -p
               CP := cp -p -v
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

#--------------------------------------------------------------------------
# Compiling flags for ROMS Applications.
#--------------------------------------------------------------------------

ifdef ROMS_APPLICATION
 ifdef USE_DEBUG
           FFLAGS += -g
           FFLAGS += -fbounds-check
           FFLAGS += -fbacktrace
           FFLAGS += -fcheck=all
           FFLAGS += -fsanitize=address -fsanitize=undefined
           FFLAGS += -finit-real=nan -ffpe-trap=invalid,zero,overflow
 else
           FFLAGS += -O3
           FFLAGS += -ffast-math
 endif
        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)
endif

#--------------------------------------------------------------------------
# Compiling flags for CICE Applications.
#--------------------------------------------------------------------------

ifdef CICE_APPLICATION
          CPPDEFS := -DLINUS $(MY_CPP_FLAGS)
 ifdef USE_DEBUG
           FFLAGS += -g
#          FFLAGS += -O3
           FFLAGS += -fbounds-check
#          FFLAGS += -fcheck=all
           FFLAGS += -fsanitize=address -fsanitize=undefined
 else
           FFLAGS := -O3 -w
 endif
endif

#--------------------------------------------------------------------------
# Library locations, can be overridden by environment variables.
#--------------------------------------------------------------------------

          LDFLAGS := $(FFLAGS)

ifdef USE_NETCDF4
        NF_CONFIG ?= nf-config
    NETCDF_INCDIR ?= $(shell $(NF_CONFIG) --includedir)
             LIBS := $(shell $(NF_CONFIG) --flibs)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
else
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
      NETCDF_LIBS ?= -lnetcdf
             LIBS := -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
           INCDIR += $(NETCDF_INCDIR) $(INCDIR)
endif

ifdef USE_HDF5
      HDF5_INCDIR ?= /opt/gfortransoft/serial/hdf5/include
      HDF5_LIBDIR ?= /opt/gfortransoft/serial/hdf5/lib
             LIBS += -L$(HDF5_LIBDIR) -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lz
           INCDIR += $(HDF5_INCDIR)
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/gfortransoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/gfortransoft/ARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -lfmpi -lmpi
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -fopenmp -static-libgcc
#            LIBS += -lgomp
endif

ifdef USE_MCT
       MCT_INCDIR ?= /usr/local/mct/include
       MCT_LIBDIR ?= /usr/local/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
           INCDIR += $(MCT_INCDIR) $(INCDIR)
endif

ifdef USE_ESMF
          ESMF_OS ?= $(OS)
      ESMF_SUBDIR := $(ESMF_OS).$(ESMF_COMPILER).$(ESMF_ABI).$(ESMF_COMM).$(ESMF_SITE)
      ESMF_MK_DIR ?= $(ESMF_DIR)/lib/lib$(ESMF_BOPT)/$(ESMF_SUBDIR)
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

ifdef USE_COAMPS
             LIBS += $(COAMPS_LIB_DIR)/libashare.a
             LIBS += $(COAMPS_LIB_DIR)/libcoamps.a
             LIBS += $(COAMPS_LIB_DIR)/libaa.a
             LIBS += $(COAMPS_LIB_DIR)/libam.a
             LIBS += $(COAMPS_LIB_DIR)/libfishpak.a
             LIBS += $(COAMPS_LIB_DIR)/libfnoc.a
             LIBS += $(COAMPS_LIB_DIR)/libtracer.a
#            LIBS += $(COAMPS_LIB_DIR)/libnl_beq.a
endif

ifdef CICE_APPLICATION
            SLIBS += $(SLIBS) $(LIBS)
endif

ifdef USE_WRF
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
             LIBS += $(WRF_LIB_DIR)/main/module_wrf_top.o
             LIBS += $(WRF_LIB_DIR)/main/libwrflib.a
             LIBS += $(WRF_LIB_DIR)/external/fftpack/fftpack5/libfftpack.a
             LIBS += $(WRF_LIB_DIR)/external/io_grib1/libio_grib1.a
             LIBS += $(WRF_LIB_DIR)/external/io_grib_share/libio_grib_share.a
             LIBS += $(WRF_LIB_DIR)/external/io_int/libwrfio_int.a
             LIBS += $(WRF_LIB_DIR)/external/esmf_time_f90/libmyesmf_time.a
             LIBS += $(WRF_LIB_DIR)/external/RSL_LITE/librsl_lite.a
             LIBS += $(WRF_LIB_DIR)/frame/module_internal_header_util.o
             LIBS += $(WRF_LIB_DIR)/frame/pack_utils.o
             LIBS += $(WRF_LIB_DIR)/external/io_netcdf/libwrfio_nf.a
     WRF_MOD_DIRS  = main frame phys share external/esmf_time_f90
 else
             LIBS += $(WRF_LIB_DIR)/module_wrf_top.o
             LIBS += $(WRF_LIB_DIR)/libwrflib.a
             LIBS += $(WRF_LIB_DIR)/libfftpack.a
             LIBS += $(WRF_LIB_DIR)/libio_grib1.a
             LIBS += $(WRF_LIB_DIR)/libio_grib_share.a
             LIBS += $(WRF_LIB_DIR)/libwrfio_int.a
             LIBS += $(WRF_LIB_DIR)/libmyesmf_time.a
             LIBS += $(WRF_LIB_DIR)/librsl_lite.a
             LIBS += $(WRF_LIB_DIR)/module_internal_header_util.o
             LIBS += $(WRF_LIB_DIR)/pack_utils.o
             LIBS += $(WRF_LIB_DIR)/libwrfio_nf.a
 endif
endif

# Use full path of compiler.

               FC := $(shell which ${FC})
               LD := $(FC)

#--------------------------------------------------------------------------
# ROMS specific rules.
#--------------------------------------------------------------------------

# Turn off bounds checking for function def_var, as "dimension(*)"
# declarations confuse Gnu Fortran 95 bounds-checking code.

ifdef ROMS_APPLICATION
 $(SCRATCH_DIR)/def_var.o: FFLAGS += -fno-bounds-check
endif

# Allow integer overflow in ran_state.F.  This is not allowed
# during -O3 optimization. This option should be applied only for
# Gfortran versions >= 4.2.

ifdef ROMS_APPLICATION
 FC_TEST := $(findstring $(shell ${FC} --version | head -1 | \
                              awk '{ sub("Fortran 95", "Fortran"); print }' | \
                              cut -d " " -f 4 | \
                              cut -d "." -f 1-2), \
             4.0 4.1)

 ifeq "${FC_TEST}" ""
  $(SCRATCH_DIR)/ran_state.o: FFLAGS += -fno-strict-overflow
 endif
endif

# Set free form format in some ROMS source files to allow long string for
# local directory and compilation flags inside the code.

ifdef ROMS_APPLICATION
 $(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -ffree-form -ffree-line-length-none
 $(SCRATCH_DIR)/mod_strings.o: FFLAGS += -ffree-form -ffree-line-length-none
 $(SCRATCH_DIR)/analytical.o: FFLAGS += -ffree-form -ffree-line-length-none
 $(SCRATCH_DIR)/biology.o: FFLAGS += -ffree-form -ffree-line-length-none

 ifdef USE_ADJOINT
  $(SCRATCH_DIR)/ad_biology.o: FFLAGS += -ffree-form -ffree-line-length-none
 endif
 ifdef USE_REPRESENTER
  $(SCRATCH_DIR)/rp_biology.o: FFLAGS += -ffree-form -ffree-line-length-none
 endif
 ifdef USE_TANGENT
  $(SCRATCH_DIR)/tl_biology.o: FFLAGS += -ffree-form -ffree-line-length-none
 endif
endif

#--------------------------------------------------------------------------
# Model coupling specific rules.
#--------------------------------------------------------------------------

# Add COAMPS library directory to include path of ESMF coupling files.

ifdef USE_COAMPS
 $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
 $(SCRATCH_DIR)/esmf_esm.o: FFLAGS += -I$(COAMPS_LIB_DIR)
endif

# Add WRF library directory to include path of ESMF coupling files.

ifdef USE_WRF
 ifeq "$(strip $(WRF_LIB_DIR))" "$(WRF_SRC_DIR)"
  $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += $(addprefix -I$(WRF_LIB_DIR)/,$(WRF_MOD_DIRS))
 else
  $(SCRATCH_DIR)/esmf_atm.o: FFLAGS += -I$(WRF_LIB_DIR)
 endif
endif

# Supress free format in SWAN source files since there are comments
# beyond column 72.

ifdef USE_SWAN
 $(SCRATCH_DIR)/ocpcre.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/ocpids.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/ocpmix.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swancom1.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swancom2.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swancom3.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swancom4.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swancom5.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanmain.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanout1.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanout2.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanparll.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanpre1.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanpre2.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swanser.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swmod1.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/swmod2.o: FFLAGS += -ffixed-form
 $(SCRATCH_DIR)/m_constants.o: FFLAGS += -ffree-form -ffree-line-length-none
 $(SCRATCH_DIR)/m_fileio.o: FFLAGS += -ffree-form -ffree-line-length-none
 $(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += -ffree-form -ffree-line-length-none
 $(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -ffree-form -ffree-line-length-none
endif
