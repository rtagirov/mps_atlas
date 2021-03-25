#FC=ifort 
FC= mpif90    
PROG=mps-atlas-run.x


LFLAGS =
#ifort preprocessor flags
FPPFLAGS = "-DMPI -DORIGL "# -DDFSY" # -DDYNBIN" 

FFLAGS = "-c"# -traceback  -heap-arrays -check bounds" 
#  -stand f90  -assume realloc_lhs  -check all  -traceback   -fstack-protector  -assume protect_parens"  #-O2 


# NETCDF library routines

INCLUDE="-I../lib/include"
NETCDFLIB="-L../NetCDF -lnet -L../lib/lib -lnetcdf"

####################################################################

SUBDIRS = NetCDF  src
 
  

all: mps-atlas-run

mps-atlas-run:
	 @for i in $(SUBDIRS) ; do \
                cd $$i ; \
                $(MAKE)         \
                        FC=$(FC) \
                        EXEC=$(PROG) \
                        FFLAGS=$(FFLAGS) \
                        FPPFLAGS=$(FPPFLAGS) \
                        FFTLIB=$(FFTLIB) \
                        INCLUDE=$(INCLUDE) \
                        NETCDFLIB=$(NETCDFLIB) \
                        STATIC=$(STATIC) \
                        all ; cd .. ;\
        done



clean:
	@for i in $(SUBDIRS); do \
                (cd $$i; \
                $(MAKE) clean ); \
        done
	rm -f *~

