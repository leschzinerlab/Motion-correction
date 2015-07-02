# To USE
# 1) Modify paths for libraries
# 2) type "make -f Makefile alignparts_lmbfgs
#         "make -f Makefile alignframes_lmbfgs
FC = gfortran

FFLAGS = -ofast -Wall -fbounds-check -g -Wno-uninitialized 

LBFGSB      = Lbfgsb.3.0/lbfgsb.f
LINPACK     = Lbfgsb.3.0/linpack.f
BLAS        = Lbfgsb.3.0/blas.f
TIMER       = Lbfgsb.3.0/timer.f
IMLIB2010   = ${IMAGE_PUBLIC}/lib/imlib2010.a
GENLIB2010  = ${IMAGE_PUBLIC}/lib/genlib.a

ALIGNFRAMES_LMBFGS = alignframes_lmbfgs.f90
ALIGNPARTS_LMBFGS  = alignparts_lmbfgs.f90
ALIGNPARTS_LMBFGS_RELION  = alignparts_lmbfgs_relion.f90
ALIGNGRID_LMBFGS   = aligngrid_lmbfgs.f90
ALIGNPARTS_LMBFGS_AMEDEE = alignparts_lmbfgs_amedee.f90
SHIFTFRAMES_LIST = shiftframes_list.f90
ALIGNPARTS_LMBFGS_RADANDMOTION = alignparts_lmbfgs_radandmotion.f90
ALIGNPARTS_STARFILEHANDLER = alignparts_starfilehandler.f90

all :  alignframes_lmbfgs alignparts_lmbfgs aligngrid_lmbfgs shiftframes_list alignparts_lmbfgs_amedee alignparts_lmbfgs_radandmotion alignparts_starfilehandler

alignframes_lmbfgs : $(ALIGNFRAMES_LMBFGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNFRAMES_LMBFGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010) -o alignframes_lmbfgs.exe -lfftw3 -lm

alignparts_lmbfgs : $(ALIGNPARTS_LMBFGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNPARTS_LMBFGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010) -o alignparts_lmbfgs.exe -lfftw3 -lm

alignparts_lmbfgs_relion : $(ALIGNPARTS_LMBFGS_RELION) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNPARTS_LMBFGS_RELION) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010) -o alignparts_lmbfgs_relion.exe -lfftw3 -lm

aligngrid_lmbfgs : $(ALIGNGRID_LMBFGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNGRID_LMBFGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010) -o aligngrid_lmbfgs.exe -lfftw3 -lm

alignparts_lmbfgs_amedee : $(ALIGNPARTS_LMBFGS_AMEDEE) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNPARTS_LMBFGS_AMEDEE) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010) -o alignparts_lmbfgs_amedee.exe -lfftw3 -lm

alignparts_lmbfgs_radandmotion : $(ALIGNPARTS_LMBFGS_RADANDMOTION) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNPARTS_LMBFGS_RADANDMOTION) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(IMLIB2010) $(GENLIB2010) -o alignparts_lmbfgs_radandmotion.exe -lfftw3 -lm

alignparts_starfilehandler :  $(ALIGNPARTS_STARFILEHANDLER) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(ALIGNPARTS_STARFILEHANDLER) $(IMLIB2010) $(GENLIB2010) -o alignparts_starfilehandler.exe -lfftw3 -lm

shiftframes_list :  $(SHIFTFRAMES_LIST) $(IMLIB2010) $(GENLIB2010)
	$(FC) $(FFLAGS) $(SHIFTFRAMES_LIST) $(IMLIB2010) $(GENLIB2010) -o shiftframes_list.exe -lfftw3 -lm
