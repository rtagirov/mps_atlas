#module load openmpi_intel/1.10.1_2017.1

COMMONS = common.abross common.abtotb \
          common.constb common.contbl common.convbl \
          common.depart common.edenbl common.elembl \
          common.fluxbl common.freqbl common.height \
          common.ifblkk common.ifequa common.ifopbl common.ionsbl \
          common.iterbl common.junkbl common.meudbl common.musblk \
          common.opsblk common.optotb common.parblk \
          common.putblk common.rhoxbl \
          common.stateb common.steplg common.freebl \
          common.taushj common.tcorrb common.odfnlte \
          common.teffbl common.tempbl common.turbpr common.waveyb \
          common.xabund common.xnfblk common.xnfpbl common.xnmolb 


#CAA filter part
AOBJECTS = mpsa.main.o   mpsa.read0.o mpsa.kapp.o \
           mpsa.pops.o   mpsa.josh.o mpsa.conv.o \
           mpsa.odfcalc.o mpsa.continuum.o mpsa.binsize.o \
           mpsa.calc.sub.o mpsa.calc.ross.o mpsa.comms.o\
           df.synthe.o df.vttab.o df.tabvoigt.o df.intervals.o \
	   df.bisort.o df.calc.o df.selectlines.o  \
           df.xlinop.o df.hprof.o df.write_odf.o \
           df.auxiliary.o df.intervals.filter.o df.calc.filter.o df.filter.o df.bisort.filter.o

OBJECTS = types.o dfcomm.o comm.atlas.o

MODULES = mpsa.kapp.o mpsa.pops.o mpsa.comms.o  

AOBJECTS_B = mpsa.main.o   mpsa.read0.o  mpsa.kapp.o \
             mpsa.pops.o   mpsa.josh.o mpsa.conv.o



# target definition
all: atlas 
# -------------------------- ATLAS ------------------------------------

atlas: $(AOBJECTS) $(OBJECTS)
	$(FC) $(FPPFLAGS) -o run-mpsa.x $(AOBJECTS) $(OBJECTS) $(INCLUDE) $(NETCDFLIB)
	mv run-mpsa.x ../bin 

#########################################################################
# ----------------------- ATLAS9 DEPENDENCIES ---------------------------
dfcomm.o : types.o dfcomm.f90
	$(FC) $(FFLAGS) dfcomm.f90

types.o : types.f90
	$(FC) $(FFLAGS) types.f90

comm.atlas.o : types.o comm.atlas.f90
	$(FC) $(FFLAGS) comm.atlas.f90

mpsa.read0.o: $(OBJECTS)  mpsa.read0.F
	$(FC) $(FFLAGS)  $(FPPFLAGS) mpsa.read0.F

mpsa.binsize.o: $(OBJECTS)  mpsa.binsize.f90
	$(FC) $(FFLAGS)  $(FPPFLAGS) mpsa.binsize.f90

mpsa.kapp.o: $(OBJECTS) mpsa.binsize.o  mpsa.kapp.F
	$(FC) $(FFLAGS)  $(FPPFLAGS) mpsa.kapp.F

mpsa.pops.o: $(OBJECTS)  mpsa.pops.F
	$(FC)  $(FFLAGS)  $(FPPFLAGS)  mpsa.pops.F

mpsa.josh.o: $(OBJECTS)  mpsa.josh.f
	$(FC)  $(FFLAGS)  $(FPPFLAGS) mpsa.josh.f

mpsa.conv.o: $(OBJECTS)  mpsa.conv.f
	$(FC)  $(FFLAGS)  $(FPPFLAGS) mpsa.conv.f

mpsa.odfcalc.o: $(OBJECTS)  mpsa.odfcalc.f90
	$(FC)  $(FFLAGS)  $(FPPFLAGS) mpsa.odfcalc.f90

mpsa.continuum.o: $(OBJECTS) mpsa.continuum.f90
	$(FC)  $(FFLAGS)  $(FPPFLAGS) mpsa.continuum.f90

mpsa.calc.sub.o: $(OBJECTS) $(MODULES) mpsa.calc.sub.f90 
	$(FC) $(FFLAGS)  $(FPPFLAGS) mpsa.calc.sub.f90

mpsa.calc.ross.o: $(OBJECTS) $(MODULES) mpsa.calc.ross.f90
	$(FC) $(FFLAGS)  $(FPPFLAGS) mpsa.calc.ross.f90 

mpsa.comms.o: $(OBJECTS) mpsa.comms.F90 
	$(FC) $(FFLAGS) $(FPPFLAGS) mpsa.comms.F90
	
mpsa.main.o: $(OBJECTS) $(MODULES) mpsa.main.F
	$(FC) $(FFLAGS)  $(FPPFLAGS) mpsa.main.F


# ----------------------- DFsynth DEPENDENCIES ---------------------------
#CAA
df.filter.o : dfcomm.o types.o df.filter.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.filter.F90

df.write_odf.o: df.write_odf.F90
	$(FC) $(FFLAGS)  $(FPPFLAGS)  df.write_odf.F90

df.vttab.o : dfcomm.o types.o df.vttab.f90
	$(FC) $(FFLAGS)  $(FPPFLAGS)  df.vttab.f90

df.tabvoigt.o : dfcomm.o types.o  df.tabvoigt.f90
	$(FC) $(FFLAGS)  $(FPPFLAGS) df.tabvoigt.f90

df.intervals.o : dfcomm.o types.o df.intervals.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.intervals.F90

df.bisort.o : dfcomm.o types.o df.bisort.f90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.bisort.f90

#CAA
df.intervals.filter.o : dfcomm.o types.o df.intervals.filter.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.intervals.filter.F90

df.bisort.filter.o : dfcomm.o types.o df.bisort.filter.f90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.bisort.filter.f90

df.calc.o : dfcomm.o types.o df.bisort.o df.calc.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.calc.F90
#CAA
df.calc.filter.o : dfcomm.o types.o df.bisort.filter.o df.calc.filter.F90
	 $(FC) $(FFLAGS) $(FPPFLAGS)  df.calc.filter.F90

df.selectlines.o : common.teffbl common.odfnlte \
 common.tempbl types.o dfcomm.o df.selectlines.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.selectlines.F90

df.synthe.o : df.synthe.F90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.synthe.F90

df.auxiliary.o : df.auxiliary.f90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.auxiliary.f90

df.xlinop.o : dfcomm.o types.o df.hprof.o df.xlinop.f90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.xlinop.f90

df.hprof.o : dfcomm.o types.o df.hprof.f90
	$(FC) $(FFLAGS) $(FPPFLAGS)  df.hprof.f90


mpsab.josh.o: $(COMMONS) mpsab.josh.f
	$(FC)  $(FFLAGS) mpsab.josh.f


# ========================================================================
clean: 
	rm -f *~ *.o *.mod runatlas* 
