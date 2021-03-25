
INCLUDE="-I/scratch/witzke/Libraries/netcdf-3.6.1/include"
NETCDFLIB="-L../NetCDF -lnet -L/scratch/witzke/Libraries/netcdf-3.6.1/lib -lnetcdf"


#FC= ifort 
OBJ = CreateNetCDF.o OpenNetCDF.o ReadNetCDF.o CloseNetCDF.o \
	WriteNetCDF.o handle_error.o  

all: libnet.a

libnet.a: $(OBJ)
	rm -f  libnet.a
	ar crv libnet.a  $(OBJ)

CreateNetCDF.o: CreateNetCDF.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) CreateNetCDF.f90 $(INCLUDE) $(STATIC)

OpenNetCDF.o: OpenNetCDF.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) OpenNetCDF.F90 $(INCLUDE) $(STATIC)

ReadNetCDF.o: ReadNetCDF.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) ReadNetCDF.F90 $(INCLUDE) $(STATIC)

CloseNetCDF.o: CloseNetCDF.f90
	$(FC) $(FFLAGS) $(FPPFLAGS) CloseNetCDF.f90 $(INCLUDE) $(STATIC)

WriteNetCDF.o: WriteNetCDF.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  WriteNetCDF.F90 $(INCLUDE) $(STATIC)

handle_error.o: handle_error.f90
	$(FC) $(FFLAGS) $(FPPFLAGS)  handle_error.f90 $(INCLUDE) $(STATIC)


clean:
	rm -f *.o core *.a *.d *~ work.pc work.pcl *.mod

