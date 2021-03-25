# General information 
This the working version of the new dfsynthe code together with the ATLAS9 code:
- All three possible calculations procedures are included and using mps.control it is possible to chose between
  a) ODF calculations, b) atmosphere model calculations, and c) emergent flux calculations
- the current version might have come bugs when several modules are executed one after another. 
- execpt the mps.control file all other input files are in the INPUT folder
- ODF calculations are parallelised in the loop over the # of T if the pre-compilation flag DMPI is set
- The input file has been adapted: if calculating the ODFs the amount of velocities and the velocities can be specified
- if user defined binning is off, then KURUZC bin gird is used, if it is on then a file bin_grid_sizes.dat is expected

# Requirements
- Intel Fortran compiler (tested with versions 2017/2018)
- openmpi Intel compiler (for use with the pre-compilation flag -DMPI) 
- NetCDF library is needed to compile the routines in ./NetCDF: netcdf-3.6.1 (can be downloaded online)  -> change the path in the ./Makefile! 
- linelists have to be downloaded and put in ./linelists (https://owncloud.gwdg.de/index.php/s/lnC5Sjs9K0U7VEt)

# Instructions
- load appropriate modules (note the openMPI module already loads the Intel compiler module)

    `module load openmpi_intel`

- download linelists

    `bash download_linelists.sh`

- download and compile the NetCDF library
    
    `bash install_netcdf`

- compile SSWOP
    
    `make`

