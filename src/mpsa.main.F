      program mpsatlas

      use types
      use atlcomm
      use kappa_cal
      use population
#ifdef MPI
      use comms
#endif  
      implicit none
!-------------------------------------------------------------------------------!
!.... THE FOLLOWING FILES ARE USED:
!        FILE 11 = ROSSELAND OPACITY TABLES TO BE READ IF ABUNDANCE NOT 1X SOLAR
!                  THIS FILE IS USED IN SUBROUTINE READIN
!        FILE 12 = MOLECULAR EQUILIBRIUM CONSTANTS
!                  THIS FILE IS USED IN SUBROUTINE READMOL
!        FILE 13 = A FILE CONTAINING A PUNCHED MODEL.  TO BE READ UP TO
!                  BUT NOT INCLUDING THE "BEGIN".  THIS FILE IS USED IN
!                  THE SUBROUTINE READIN
!        FILE 15 = INPUT (for general variables)
!        FILE 16 = OUTPUT
!        FILE 17 = MODEL AND OR FLUX OUTPUT
!        FILE 18 = OUTPUT OF <JNU> FOR USE IN NON-LTE LINE TRANSFER
!
!
!        FILE 19 = LINE DISTRIBUTION FUNCTION INPUT WITH A SINGLE VALUE
!                  OF THE MICROTURBULENCE VELOCITY.  IF MICROTURBULENCE
!                  IS TO BE ALLOWED TO VARY WITH DEPTH, THEN A SEPARATE
!                  ODF IS NEEDED FOR EACH VELOCITY.  THIS LEADS TO 
!        FILE 20 = ODF AT 0 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]0
!        FILE 21 = ODF AT 1 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]1
!        FILE 22 = ODF AT 2 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]2
!        FILE 24 = ODF AT 4 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]4
!        FILE 28 = ODF AT 8 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]8
!        FILE 33 = actual atm model file with several models in one `row'
!
!----------------------------- COMMONS ---------------------------------
!


      include 'common.odfnlte'
      include 'common.constb'
      include 'common.abross'
      include 'common.abtotb'
      include 'common.contbl'
      include 'common.convbl'
      include 'common.depart'
      include 'common.edenbl'
      include 'common.elembl'
      include 'common.fluxbl'
      include 'common.height'
      include 'common.ifblkk'
      include 'common.ifopbl'
      include 'common.ionsbl'
      include 'common.iterbl'
      include 'common.junkbl'
      include 'common.musblk'
      include 'common.opsblk'
      include 'common.optotb'
      include 'common.putblk'
      include 'common.stateb'
      include 'common.steplg'
      include 'common.taushj'
      include 'common.tcorrb'
      include 'common.teffbl'
      include 'common.tempbl'
      include 'common.turbpr'
      include 'common.waveyb'
      include 'common.xabund'
      include 'common.xnfblk'
      include 'common.xnfpbl'
      include 'common.xnmolb'
      include 'common.freebl'

#ifdef MPI

      integer :: ncount, sendlog(3) 
      character(len=3) :: rnum
      double precision :: sendbuf(26),  time1, time2

!--- model parallelisation
      integer :: oldnulo,  oldnuhi, pprank, points, nfrcount 

#else 
      integer :: oldnulo 
      integer :: rank, comm
#endif
      

!
!--------------------------- LOCAL VARIABLES --------------------------
!
      double precision  contin, freq15, part(maxd, 6), rco,  
     &                  rcowt, stepwt, sumwt, wave, x, rosstab, asixth

! for model calculations
      double precision  taustd(maxd),  tmin, trat, trat4
      double precision dum1(maxd), dum2(maxd), dum3(maxd), dum4(maxd),
     &                 dum5(maxd), dum6(maxd), dum7(maxd), dum8(maxd)
      integer  i, j, mm,  mode1, n, nsteps, nu, nmod, imod, finum 
      integer  idum, map1  
      logical  finish, more, readi0, stopfl, cast_version 

      double precision introssk


      real(kind=8), allocatable :: stau(:), st(:), sp(:) 
      real(kind=8) :: maxTch1, maxTch2
      integer :: indexT

      integer :: structnum
!
!-------------------------------  EXTERNALS ----------------------------
!
      external convec, high,  josh, 
     &         output, outpt2, outpt3, outpt4, outpt5, outpt6,
     &         radiap, radiap2, radiap3, readi0, ross, ross2, 
     &         ross3, stateq, stateq2, stateq3, tcorr, tcorr2, tcorr3, 
     &         turb, rosstab, introssk 


      character lit*1

!------------------------------- CONTROL ------------------------------
      integer flagODF, flagmodel, flagflux
      

!--------------------------- INITIALIZATION ----------------------------

      itemp = 0 
      cast_version = .true. 
      recalxne = .false. 

! --- if parallel calculations of emergent flux, along frequencies: 

      freqsplit = .false.
      oldnulo = 1 

!------------------------------- EXECUTION -----------------------------
!
!.... INPUT SECTION
!     PREFIX P  PRESSURE
!     PREFIX T  TEMPERATURE
!     PREFIX X  ABUNDANCE FRACTION
!     PREFIX F  IONIZATION FRACTION
!     PREFIX R  FREQUENCY INTEGRAL OR INTEGRATION COEFFICIENT
!     PREFIX A OR AB  MASS ABSORPTION COEFFICIENT
!     PREFIX XNFP  NUMBER DENSITY OVER PARTITION FUNCTION
*
!     ABUND     = THE NORMALLY ASSUMED ABUNDANCES
!     ALPHA     = THE FRACTION OF OPACITY CAUSED BY SCATTERING
!     BHYD      = STATISTICAL EQUILIBRIUM FACTOR FOR HYDROGEN
!     BMIN      = STATISTICAL EQUILIBRIUM FACTOR FOR HMINUS
!     ELEM      = THE LETTER CODES FOR ELEMENTS
!     FREQID    = A LABEL FOR THE FREQUENCY SET
!     IFCORR    = TEMPERATURE CORRECTION .TRUE. = ON OR .FALSE. = OFF
!     IFEMPIR   = FALSE = DEFAULT - A THEORETICAL MODEL
!               = TRUE TO USE AN EMPIRCAL MODEL
!     IFMOL     = TRUE - SET UP EQUILIBRIUM EQUATIONS FOR NUMBER DENSITIES
!               = FALSE -  ASSUME NO MOLECULES AND ITERATE FOR NUMBER DENSITIES
!     IFPRES    = PRESSURE INTEGRATION .TRUE. = ON OR .FALSE. = OFF
!     IFPRNT    = 0 DO NOT PRINT ANYTHING 
!               = 1 PRINT MINIMAL SUMMARY TABLE AND <JUN(J)>
!               = 2 PRINT ALL FREQUENCY INDEPENDENT DATA
!               = 3 PRINT SNU, TAUNU, JNU, ETC.
!               = 4 PRINT OPACITIES
!               = 5 PRINT <JNU>
!     IFPNCH    = 0 DO NOT PUNCH  (no longer makes sense)
!               = 1 PUNCH STRUCTURE 
!               = 2 PUNCH STRUCTURE AND SURFACE FLUX OR INTENSITY
!               = 5 PUNCH 2 AND MOLECULAR NUMBER DENSITIES/PART FNS
!     IFSURF    = 0  CALCULATE FLUX FOR EVERY DEPTH.  THIS IS THE DEFAULT.
!               = 1  CALCULATE FLUX AT SURFACE ONLY
!               = 2  CALCULATE INTENSITY AT SURFACE, SUPPLY NMU ANGLES
!     IFWAVE    = TRUE- STEP NUMNU WAVELENGTHS STARTING AT WBEGIN BY WSTEP
!     LODF      = THE SWITCH THAT TELLS HOW THE LINE OPACITY DISTRIBUTION 
!                 FUNCTION WILL BE USED.  THE OPTIONS ARE: 1) "CONSTANT"
!                 FOR A CONSTANT VTURB WITH DEPTH, 2) "MEMORY  " FOR A 
!                 CONSTANT VTURB ALL IN MEMORY AT ONE TIME, AND 3) "VARIABLE"
!                 FOR A VARIABLE VTURB WITH DEPTH.  THE DEFAULT IS "CONSTANT".
!     NLTEON    = .FALSE. FOR LTE
!               = .TRUE. FOR NLTE
!     NUHI      = NUMBER OF THE FREQUENCY AT WHICH INTEGRATION STOPS
!     NULO      = NUMBER OF THE FREQUENCY AT WHICH INTEGRATION STARTS
!     NUMITS    = NUMBER OF ITERATIONS
!     NUMNU     = NUMBER OF FREQUENCIES IN THE FREQUENCY SET
!     RCOSET    = INTEGRATION COEFFICIENTS FOR THE FREQUENCIES IN FRESET
!     XABUND    = THE ABUNDANCES USED IN THE MODEL
!     XSCALE    = A SCALING FACTOR FOR METAL ABUNDANCES
!     NMOD      = number of models to be calculated
*
!---------------------------  -------------------------------------------
! --------- FIRST RUN all initialisation processes : 
!------------------------------------------------------------------------!

      comm = 0
#ifdef MPI
!      Initialise MPI communications 

       call initialise_comm
       time1 = MPI_WTIME()

#endif

!---- initialise what modules need to be calculated:

      call initialise_param(comm, flagODF, flagmodel, flagflux)

      if (flagmodel .eq. 1) mode1 = 1
      if (flagmodel .eq. 0) mode1 = 2


!---------------------------------------------------------------------
! check if we want to calculate ODFs, if so call the corresponding 
! routines
!----------------------------------------------------------------------
      
      if (flagodf .eq. 1) then 

      iterat = 1 
#ifdef MPI
      time1 = MPI_WTIME()
      if (rank .eq. 0) then 
! first calculate the condinuum and xnfdopmax 
      call calc_odf
!       call calc_continuum
!       call mpi_abort(comm, ierr)

! destribute the information of turb velocities that needs
! to be calculated to all processors 

      sendlog(1) = 1
      if (.not. ifkbin) sendlog(1) =0

      sendlog(2) = 0
      if (ifilter .eq. 1) sendlog(2) = 1

      sendlog(3) = neediv 


      do j = 1, numvt+1
       if (j .eq. 1) then 
          sendbuf(j) =dble(numvt) 
       else
          sendbuf(j) = ivt(j-1)
       end if
      end do 

      end if 
!Bcast
      call MPI_BCAST(sendlog, 3, MPI_INTEGER, 0, comm, ierr)

      if (rank .ne. 0) then 
        ifkbin = .false.
        if (sendlog(1) .eq. 1) ifkbin = .true.

        ifilter = 0 
        if (sendlog(2) .eq. 1) ifilter = 1

        neediv = sendlog(3)

      end if

      ncount = 26  
      call MPI_BCAST(sendbuf,ncount, MPI_DOUBLE_PRECISION,0,comm,ierr)
    

      if (rank .ne. 0) then 
       numvt = int(sendbuf(1))
       do j = 1, numvt
        ivt(j) = sendbuf(j+1)
       end do 
      end if
 

      time2 = MPI_WTIME()

      time2 = time2-time1
      if(rank .eq.0) print*,'time for pre-calculations  t1 = ',time2

      call MPI_Barrier(comm, ierr) 
! start the actual ODF calculations ...
!      call mpi_abort(comm, ierr)  
      call dfsynthe(rank, sizee, comm, ODFfile) 

      call MPI_Barrier(comm,ierr)

#else
! first get continuum and xnfdopmax 
      call calc_odf
! start the actual ODF calculations
      call dfsynthe(0, 1, 0, ODFfile) 
#endif

      
      end if 
!     ODF calculations finished!
!-------------------------------------------------------------------

!--- Here either an atmosphere model or emergent flux is calculated,
! or first a model, then it returns to 10000 and restarts with the 
! calculated model to get the flux 
! ---------------------------------

! initialiaze:

10000 itemp = 0  

! --- note, iter is used, when flux only is calculated
! --- iterat is used if there is more than 1 iteration, for example when a model is calculated

      iter = 1
      iterat = 1


      if ((flagmodel .eq. 1)  .or. (flagflux .eq. 1) ) then 

        asixth = 1.26d-15

          
        if (mode1 .eq. 1) then 
          INPUTfile = filemodelinput 
          modelfile = filestartmodel 
        else
          INPUTfile = filefluxinput 
          modelfile = filefluxmodel  
        endif
 


 
#ifdef MPI


!--------------------------------------------------------!
!.... OPEN THE I/O FILES
!
      open (unit = 15, file = INPUTfile, form = 'formatted',
     &      status = 'old', access= 'stream',  err = 7100)
!
      open (unit = 33, file = modelfile, form = 'formatted',
     &      status = 'old',  access = 'stream',   err=7200)

      write (rnum, '(I0)') rank

      open (unit = 16, file = 'mpsa.'//trim(rnum)//'print',
     & form = 'formatted', status = 'new', access = 'sequential')
      open (unit = 17, file = 'mpsa.'//trim(rnum)//'punch',
     & form = 'formatted', status = 'new', access = 'sequential')
      open (unit = 18, file = 'mpsa.'//trim(rnum)//'jnu',
     & form = 'formatted', status = 'new', access = 'sequential')






#else

      rank = 0
!--------------------------------------------------------!
!.... OPEN THE I/O FILES
!
      open (unit = 15, file = INPUTfile, form = 'formatted',
     &      status = 'old', access = 'sequential', err = 7100)

      call cut_ray_out_of_slice(modelfile)

      open (unit = 33, file = modelfile, form = 'formatted',
     &      status = 'old', access = 'sequential', err = 7200)

      open (unit = 16, file = 'mpsa.print', form = 'formatted',
     &      status = 'new', access = 'sequential')
      open (unit = 17, file = 'mpsa.punch', form = 'formatted',
     &      status = 'new', access = 'sequential')
      open (unit = 18, file = 'mpsa.jnu', form = 'formatted',
     &      status = 'new', access = 'sequential')
#endif

*.... THE ODF FILES ARE OPENED IN READIN WHEN IT HAS IFOP TO CHECK
*.... THEY ARE ALSO CLOSED IN READIN WHEN IT SEES "END", EXCEPT FOR
*.... THE CASE WHEN THE WHOLE ODF IS IN MEMORY, IN WHICH CASE THE
*.... FILE IS CLOSED IN LINOPM AFTER THE CONTENTS ARE READ IN.
*
*.... PUT THE WHOLE THING IN A LOOP TO AVOID A GO TO AT THE END
      more = .true.
      lit = 'y'
!!! DO LOOP (1)------------------------------------------------
!-------------------------------------------------------------*
      do while (more)

            
         more = readi0 (mode1, rank, comm)
!

         if((more) .and. (itemp .eq. 0)) then
            read (33, *) nmod
!            read (33,*) nrhox
            read(33, *) structnum, nrhox

            print*, 'nmod = ', nmod

            print*, 'structnum = ', structnum, 'nrhox = ', nrhox

!     WRITE TITLE INFO TO OUTPUT FILES 
!
#ifdef MPI

        if (rank .eq. 0)  call print_summary(mode1)          
        call output(rank)

        if (nmod .gt. 1) then  
           if ((mode1 .eq. 2) .and. (mod(nmod, sizee) .ne. 0)) then
             print*, '# of models can not be distributed amount cores'
             print*, ' applications will be aborted'
             call mpi_abort(comm, error, ierr )
           elseif (mode1 .eq. 2) then 
             nmod = nmod/sizee
            end if
        endif 

#else 
            call print_summary(mode1)
            call output(0)
#endif
       do while( iterat .le. numits) 
         if (rank .eq. 0)  print*, 'this is iteration  =', iterat
         
!.... BEGINNING OF THE MODEL LOOP
!   
! nmod = number of models ------------------------------------!

        do imod = 1, nmod

#ifdef MPI
           if (rank .eq. 0 ) then
               open (unit = 41, file = 'mpsa.wave',
     &        form = 'formatted', 
     &        access = 'sequential')

           endif
      finum = 1
      if ((mode1 .eq. 2) .and. (nmod .gt. 1))  finum = rank*nmod+imod


#else
      if ((iterat .eq. 1) .and. (imod .eq. 1)) then   
      finum = 1
               open (unit = 41, file = 'mpsa.wave',
     &        form = 'formatted', status = 'new',
     &        access = 'sequential')
      endif

#endif

! 
!     READ ATMOSPHERE MODEL
!--------------------------------------------------------------!

      if ((imod .eq. 1) .and. (iterat .eq. 1)) then

        do mm = 1, finum
               read (33, *) (rhox(j), t(j), p(j), xne(j),
     &            abross(j), accrad(j), vturb(j), j = 1, nrhox)
       
        end do

      else if (imod .gt. 1) then 

            read(33, *) structnum, nrhox

            print*, 'structnum = ', structnum, 'nrhox = ', nrhox

            if (nrhox .gt. maxd) then
               write( *,*) ' nrhox = ', nrhox, ', which is > maxd'
               stop
            endif

            read (33, *) (rhox(j), t(j), p(j), xne(j),
     &            abross(j), accrad(j), vturb(j), j = 1, nrhox)

      end if


      if (iterat .eq. 1) then 

!---------- calculate rosseland tau--------------------------------------!
!       call integ (rhox, abross, tauros, (nrhox), (abross(1)*rhox(1)))
!     
!    Check if abross was  PROPERLY CALCULATED  IN THE MODEL ATMOSPHERE
!
               if (abross(nrhox).lt. 1e-08) then
                  do j = 1, nrhox
                     abross(j) = rosstab(t(j), p(j), vturb(j))
                  end do
               end if


*
*.... ADD HERE SOME BOUNDARY PRESSURE INITIALIZATIONS 
*
               pradk0 = asixth * teff**4
!               pradk0 = 0.0d0
               pturb0 = pturb(1)
               pcon   = 0.0d0
               pzero  = pcon + pturb0+  pradk0
*
!     SHOULD PROBABLY CHANGE THE FOLLOWING AND ALLOW FOR PRADK0 NOT BEING GIVEN
!         could use the T^4 prescription similar to readin routine 
!
!Prad is set initially to zero and is calculated here:

               call integ (rhox, accrad, prad, (nrhox),
     &                     (accrad(1) * rhox(1)))
!
!-----------------------------------------------------------------!
!   if we do model calculations, then we need to re-scale the model!
!   --- this is only done on rank = 0, and the solution is passed to all nodes!

      if ( mode1 .eq. 1) then 
       if (rank .eq. 0 ) then 
           do j = 1, krhox
            taustd(j)=exp(2.30258509299405d0*(tau1lg+dble(j-1)*steplg))
           end do
!
            call integ (rhox, abross, tauros, (nrhox),
     &                  (abross(1)*rhox(1)))

! check if the deepest point in the tau grid is to 
! far away from the deepest  point in the initial starting atmosphere 
!  or larger than tau = 100

            if( (taustd(krhox) .gt. 300.0d0) 
     &           .or. (taustd(krhox) .gt. 2.0d0*tauros(nrhox))) then 
               print*, 'adjust tau grid, taustd (krhox):',taustd(krhox) 
#ifdef MPI
               call mpi_abort(comm, ierr) 
#else
               stop
#endif
            end if 


            tauros(1) = min(tauros(1), taustd(1))
            idum = map1(tauros, rhox,   (nrhox), taustd, dum1, (krhox))
            idum = map1(tauros, t,      (nrhox), taustd, dum2, (krhox))
            idum = map1(tauros, p,      (nrhox), taustd, dum3, (krhox))
            idum = map1(tauros, xne,    (nrhox), taustd, dum4, (krhox))
            idum = map1(tauros, abross, (nrhox), taustd, dum5, (krhox))
            idum = map1(tauros, prad,   (nrhox), taustd, dum6, (krhox))
            idum = map1(tauros, vturb,  (nrhox), taustd, dum7, (krhox))
            idum = map1(tauros, bmin,   (nrhox), taustd, dum8, (krhox))
!
            do j = 1, krhox
               rhox(j)   = dum1(j)
               t(j)      = dum2(j)
               p(j)      = dum3(j)
               xne(j)    = dum4(j)
               abross(j) = dum5(j)
               prad(j)   = dum6(j)
               pradk(j)  = prad(j) + pradk0
               vturb(j)  = dum7(j)
               bmin(j)   = dum8(j)
            end do
!
!           
!
            do i = 1, 8
               idum = map1 (tauros, bhyd(1,i), (nrhox), taustd,
     &                      dum1, (krhox))
!
               do j = 1, krhox
                  bhyd(j,i) = dum1(j)
               end do
!
            end do


!
            nrhox = krhox

!
!
!----- before we continue do we need per-tab rosseland mean?

         call calcross
         print*, 'Rosseland mean calulated on the flight'
!-----------------------------------------------
!

!
            if ((teff1 .lt. teff - 1.0)    .or.
     &               (teff1 .gt. teff + 1.0)    .or.
     &               (gnew  .lt. grav * 0.999)  .or.
     &               (gnew  .gt. grav * 1.001)) then
               trat = teff1 / teff
               trat4 = trat ** 4
!
               do j = 1, nrhox
                  t(j)      = t(j) * trat
                  tauros(j) = taustd(j)
                  prad(j)   = prad(j) * trat4
                  pradk(j)  = pradk(j) * trat4
                  pturb(j)  = 0.0d0
               end do
!
               pradk0 = pradk0 * trat4
               pzero  = pcon + pradk0 + pturb0
               teff   = teff1
               flux = sigma / fourpi * teff ** 4
               grav   = gnew
               glog   = log10(grav)
               call ttaup (t, tauros, abross, ptotal, p, prad, pturb,
     &                     vturb, (grav), (nrhox))
!
               do j = 1, nrhox
                  rhox(j)   = ptotal(j) / grav
                  ptotal(j) = ptotal(j) + pzero
               end do

            end if
       end if !----- rank .eq. 0 

!------- now communicate all the atmosphere strucutre info to all core:
!-- only once in the beginning: pradk0, pzero, teff, flux, grav, glog, vturb, nrhox
#ifdef MPI
       if (sizee .gt. 1)  call begin_comm
#endif 
! allocate saved t, p, tau structure
           allocate(stau(nrhox) ) 
           allocate(st(nrhox))
           allocate(sp(nrhox))
           stau =0.0d0
           st = 0.0d0
           sp = 0.0d0
      end if !--- mode1 .eq. 1 
      end if !--- iteration .eq. 1

!------- now communicate all the atmosphere strucutre info to all core: 
! t, rhox, p, ptotal, pturb   
#ifdef MPI
       if (mode1 .eq. 1 .and. sizee .gt.1) then 
          call startiteration
       endif 
#endif 

!-----------------------------------------------------------------------------------       
!      
!
!
!----- set all temperature dependent quantities to the right values for the calculation

        do j = 1, nrhox
            pradk(j)  = prad(j) + pradk0
        end do
!
       do j = 1, nrhox
          tk(j)     = k * t(j)
          hckt(j)   = hc / tk(j)
          hkt(j)    = h / tk(j)
          tkev(j)   = k_ev * t(j)
          tlog(j)   = log(t(j))
          xnatom(j) = p(j) / tk(j) - xne(j)
          rho(j)    = xnatom(j) * wtmole * 1.660d-24
          if (ifturb) pturb(j) = 0.5 * rho(j) * vturb(j) ** 2
       end do
!--- we need to check a different convergence criteria for model calculations 

       if (mode1 .eq. 1) then 

         maxTch2 = 0.0d0
      
         do j = 1, nrhox
          maxTch1 =abs(t(j) - st(j))/t(j)

          if (maxTch2 .lt. maxTch1 ) then
            maxTch2 = maxTch1
            indexT = j
          endif

          st(j) = t(j)
          stau(j) = tauros(j)
          sp(j)  = p(j)


         end do 

         if (maxTch2 .lt. 1.0d-5) then 
           print*, 'model converged, rank, iterat = ',rank,  iterat
           iterat = numits
          end if  

!-----   to check from time to time we write out the atmosphere structure
         if ((mod(iterat ,10) .eq. 1)  .and. (rank .eq. 0)) then
           write(91,*) 'strucutre after iteration =  ',  iterat-1 
           do j = 1, nrhox
             write(91+rank,911) rhox(j), t(j), p(j), xne(j),
     &            tauros(j), prad(j), vturb(j)
911       format(1PE14.7,0PF9.1,1P8E11.3,0PF12.3,F9.3)
           end do
         endif

        endif ! --- close if(model1 .eq. ) 
!--------------------------------------------------------------!
           if(mode1 .eq. 2) then 
               iter  = imod 
           elseif (mode1 .eq. 1) then 
               iter = iterat
           endif 

!
!.... CHANGING ITEMP TELLS THE SUBROUTINES THEY HAVE A NEW TEMPERATURE
!
               itemp = itemp + iter
!  only if pressure it on the hydrostatic eq. is solved
!--------------------------------------------------------------!

               if(ifpres) then
                  call calc_press_hydro
                  ifedns = .false.
               end if

!-----   solve Statistical equilibrium 
              
               call calc_se

!.... ERASE FREQUENCY INTEGRALS

               if (ifcorr) call tcorr
               call ross
               call radiap
               if (nlteon) call stateq
            


#ifdef MPI
!   --- split the frequency loop across the cores for model calculations 
!   ---  of for flux calculations on many cores

!      if((iterat .eq. 1) .and. (sizee .gt. 1) .and. (mode1 .eq. 1)) then
      if((iterat .eq. 1) .and. (sizee .gt. 1) .and. (nmod .eq.1)) then  
         oldnuhi = nuhi
         oldnulo = nulo 
         points = nuhi-nulo+1
         pprank = int(points/sizee)
         nulo = rank*pprank+nulo 
         nuhi = nulo + pprank -1 
         if(rank .eq. sizee-1)  nuhi = oldnuhi  
         if (mode1 .eq. 2) then 
           freqsplit = .true.
           call set_splitf(points, nmu ) 
           nfrcount = points*nmu

         endif 
      end if  
      
#endif
!
!.... FREQUENCY INTEGRATION SECTION 
!           ( big do loop from nulo, nuhi)      !
!
               do nu = nulo, nuhi
!------------ if wave open ---------------------------------------
!----------------------------------------------------------------*
                  if (ifwave) then
!
                     if (wbegin .le. 1.0d10) then
                        wave = wbegin + dble(nu-nulo) * deltaw
                        freq = 2.997925d17 / wave
                        rco  = abs (deltaw / wave * freq)
!
                     else
!
!.... EQUALLY SPACED FREQUENCIES
!
                        freq = wbegin + dble(nu-nulo) * deltaw
                        rco  = deltaw
                     end if
!
                  else
                     freq = freset (nu)
                     rco  = rcoset (nu)


                  end if
!-------------------------------------------------------------------*
! --------  close if wave ------------------------------------------


                  freqlg = log10(freq)
                  freqln = log(freq)
                  freq15 = freq * 1.0d-15
                  waveno = freq / c_cm
                  wave =  2.997925d17 / freq 

                  if (freqsplit) glwave(nu - oldnulo +1) = wave 

                  do j = 1, nrhox
                     ehvkt(j)  = exp(-freq * hkt(j))
                     stim(j)   = 1.0d0 - ehvkt(j)
                     bnu(j)    = 1.47439d-2 * freq15**3 * ehvkt(j) /
     &                           stim(j)
                     dbnudt(j) = bnu(j) * freq * hkt(j) / t(j) / 
     &                           stim(j)
                     if(numnu .eq. 1) dbnudt(j) = 4.0d0 * sigma / pi * 
     &                                            t(j) ** 3
                  end do
*
!----------- opacity switches no 15 and no 16 -------------------
                  if (.not. ifop(15) .and. .not. ifop(16)) then
*
*.... THIS SECTION FOR NO METALLIC LINE BLANKETING
*
                     nsteps  = 1
                     stepwt = 1.0d0
                     call kapp
! 
                     if (ifop(14) .or. ifop(17) .or. ifop(18)) then
*
*.... FOR THE CASE OF HYDROGEN LINE BLANKETING
*
                        if (ifop(14)) call hlinop
*
*.... FOR THE CASE OF EXTRA LINE ABSORPTION
*
                        if (ifop(17)) call xlinop
*
*.... FOR THE CASE OF EXTRA LINE SCATTERING
*
                        if (ifop(18)) call xlisop
*
                        do j = 1, nrhox
                           aline(j) = ahline(j) + alines(j) + axline(j)
                           sline(j) = bnu(j)
                           if (aline(j) .gt. 0.0d0) sline(j) =
     &                       (ahline(j) * shline(j) + 
     &                        alines(j) * bnu(j) + 
     &                        axline(j) * sxline(j)) / aline(j)
                           sigmal(j) = siglin(j) + sigxl(j)
                        end do
*
                     end if
*
                     call josh

*
*.... CHECK FOR ZEROS
*
                     stopfl = .false.
*
                     do j = 1, nrhox
*
                        if (snu(j) .lt. 0.0d0) then
                           write (16, *) 'snu .lt. 0.0 at j =', j
                           stopfl = .true.
                        end if
*
                        if (jnu(j) .lt. 0.0d0) then
                           write (16, *) 'jnu .lt. 0.0 at j =', j
                           stopfl = .true.
                        end if
*
                        if (hnu(j) .lt. 0.0d0) then
                           write (16, *) 'hnu .lt. 0.0 at j =', j
                           stopfl = .true.
                        end if
*
!                        if (stopfl) stop
                        if (stopfl) then
                          print *, ' snu, jnu or hnu < 0 at ', j
                          stop
                        end if
                     end do
*
                     if (ifsurf .eq. 2) then
                        hnu(1) = surfi(1)
                        contin = surfi(1)
*
                     else
                        contin = hnu(1)
                     end if
*
                     put = contin
                     call outpt2(rank) 
                     rcowt = rco * stepwt
*
                     if (ifsurf .eq. 0) then
                        if (ifcorr) call tcorr2 ((rcowt))
                        call radiap2 ((rcowt))
                        call ross2 ((rcowt))
                        if (nlteon) call stateq2 ((rcowt))
                     end if

*
*.... THIS PASSES VALUE OF STEPWT TO OUTPUT
*
                     put  = stepwt
                     iput = nsteps
                     call outpt3(rank, nu-oldnulo+1) 
                     call outpt4(rank) 
*
!----------------------------------------------------
                  else
*
*.... THIS SECTION IS FOR THE LINE BLANKETED MODELS
*.... FIRST THE CONTINUUM OPACITIES
*
                     call kapp
!
*
*.... THIS NEXT LOOP IS NEEDED TO INITIALIZE ALINE, ... EACH NU
*
                     do j = 1, nrhox
                        aline(j)  = 0.0d0  !alines(j)
!--------------------------------------------------------------
! !fudge line opacities to study effect of lines in certain 
! ! wavelength intervals
! 
!            if ((waveno .le. 10000.0) ) then
!                         aline(j) = aline(j)*2.0d0
!                        end if


                        sline(j)  = bnu(j)
                         if (cast_version ) then
                           if (aline(j) .gt. 0.0d0) sline(j) =
     &                       (ahline(j) * shline(j) +
     &                        alines(j) * bnu(j) +
     &                        axline(j) * sxline(j)) / aline(j)
                         end if


                        sigmal(j) = siglin(j) + sigxl(j)
                     end do
*
!                     write(184, *) 'continuum at wave = ', wave
!                     write(184, *) wave

                     call josh

*
                     if (ifsurf .eq. 2) then
                        hnu(1) = surfi(1)
                        contin = surfi(1)
*
                     else
                        contin = hnu(1)
                     end if
*
                     put = contin
                     call outpt2(rank) 
*
                     if (ifsurf .eq. 0) then
*
                        do j = 1, nrhox
                           abtotc(j)  = abtot(j)
                           alphac(j)  = alpha(j)
                           dtaunuc(j) = dtaunu(j)
                           hnuc(j)    = hnu(j)
                           jnuc(j)    = jnu(j)
                           residc(j)  = 0.0d0
                           snuc(j)    = snu(j)
                           taunuc(j)  = taunu(j)
                        end do
*
                     end if
*
                     sumwt = 0.0d0
*
*.... METALLIC LINE BLANKETING
*
                     finish = .false.
                     n = 0
*
                     do while (.not. finish)
                        n = n + 1
!------ --  here are the ODFs read (three different possibilities:
!          a) once and, b) from memory c) several files if microtrubulence varies 
!           throughout atmosphere
!- --------------------------------------------------------------------------
! 
                        if (ifop(15) .and. lodf .eq. 'constant') then
               
                           if (ifkbin ) then 
                             ! old routines used for mpsa.odf in binary
                             ! for that ifkbin has to be set ture

                             if (freqid .eq. 'little') then
                              call linopl ((n), nsteps, stepwt)
                             else
                              call linopb ((n), nsteps, stepwt)
                             end if
                           else 
                             ! .nc files with defined grid were already read in mpsa.read
                             ! here we only interpolate on the T,P grid

                             call interPT(n, stepwt)
                             nsteps = isubbin

                           endif 

                        else if (ifop(15) .and. 
     &                           lodf     .eq.  'memory  ') then
                           call linopm ((n), nsteps, stepwt)
*
                        else if (ifop(15) .and. 
     &                           lodf     .eq.  'variable') then
                           call linopv ((n), nsteps, stepwt)
                        end if

                        do j = 1, nrhox
                           aline(j)  = alines(j)
!--------------------------------------------------------------
! !fudge line opacities to study effect of lines in certain 
! ! wavelength intervals
! !                         
!            if ((waveno .le. 10000.0)) then
!                             aline(j) = aline(j)*2.0d0
!                           end if

                           sline(j)  = bnu(j)

                           if (cast_version ) then 
                             if (aline(j) .gt. 0.0d0) sline(j) =
     &                       (ahline(j) * shline(j) +
     &                        alines(j) * bnu(j) +
     &                        axline(j) * sxline(j)) / aline(j)
                           end if 

                           sigmal(j) = siglin(j) + sigxl(j)
                        end do
!

!                        write(184, *) 'lines at wave = ', wave
!                        write(184, *) 'next subbin number = ', n
!                        write(184, *) wave
!                        write(184, *) n

                        call josh

                        sumwt = sumwt + stepwt
                        rcowt = rco * stepwt
                        residc(1) = 1.0d0
*
                        if (ifsurf .eq. 0) then
*
                           do j = 1, nrhox
                              if (hnuc(j) .gt. 0.0d0) then
                                 residc(j) = hnu(j) / hnuc(j)
                              else
				  residc(j) = 1.0d0
                              end if
                           end do
*
                        else if (ifsurf .eq. 1 .and. 
     &                           contin .gt. 0.0d0) then
                           residc(1) = hnu(1) / contin
*
                        else if (ifsurf .eq. 2) then
                           hnu(1) = surfi(1)
                           if (contin .gt. 0.0d0) residc(1) = 
     &                                            surfi(1) / contin
                        end if
*
                        if (stepwt .gt. 0.0d0 .and. 
     &                      ifsurf .eq. 0) then
                           if (ifcorr) call tcorr2 ((rcowt))
                           call radiap2 ((rcowt))
                           call ross2 ((rcowt))
                           if (nlteon) call stateq2 ((rcowt))
                        end if
*
                        put  = stepwt
                        iput = nsteps 
                        call outpt3(rank, nu-oldnulo+1) 
*
*.... CHECK TO SEE IF THIS FREQUENCY IS FINISHED
!
!       we do not want to ignore any of the sub-bins, not matter what
                        if (ifilter .eq. 1 .or. ifilter .eq. 0) then
                          residc(1) = 0.0d0
                        end if 


                        if (n .eq. nsteps) then
                           finish = .true.
                           n = 0
                        else if (residc(1) .gt. 0.9995d0) then
                           finish = .true.
                           n = 0 
                        else if (residc(1) .gt. 0.998d0 .and. 
     &                           ifsurf    .gt. 0)    then
                           finish = .true.
                           n = 0 
                        else if (ifsurf .eq. 0) then
                           j = 1
                           do while (   j      .le. nrhox .and.
     &                               residc(j) .ge. 0.998d0)
                              j = j + 1
                           end do
                           if (j .gt. nrhox) then 
                             finish = .true.
                             n = 0 
                           end if

                        end if

*
                     end do
*
*.... FINISH OFF THIS FREQUENCY INTERVAL
*
                     if (ifilter .eq. 1 .or. ifilter .eq. 0 ) then 
                        stepwt = 0.0d0 
                     else 
                       residc(1) = 1.0d0
                       stepwt    = 1.0d0 - sumwt
                       if (stepwt .lt. 0.0001d0) stepwt = 0.0d0
                     end if 
*
                     if (ifsurf .eq. 0) then
*
                        do j = 1, nrhox
                           abtot(j)  = abtotc(j)
                           alpha(j)  = alphac(j)
                           dtaunu(j) = dtaunuc(j)
                           hnu(j)    = hnuc(j)
                           jnu(j)    = jnuc(j)
                           snu(j)    = snuc(j)
                           taunu(j)  = taunuc(j)
*
                           if (hnuc(j) .gt. 0.0d0) then
                              residc(j) = hnu(j) / hnuc(j)
                           else
                              residc(j) = 1.0d0
                           end if
*
                        end do
*
                     else if (ifsurf .eq. 1) then
                        hnu(1) = contin
*
                     else
                        surfi(1) = contin
                     end if
*
                     sumwt = sumwt + stepwt
                     rcowt = rco * stepwt
*
                     if (stepwt .gt. 0.0d0 .and. ifsurf .eq. 0) then
                        if (ifcorr) call tcorr2 ((rcowt))
                        call radiap2 ((rcowt))
                        call ross2 ((rcowt))
                        if (nlteon) call stateq2 ((rcowt))
                     end if
*
                     put  = stepwt
                     iput = nsteps 
                     call outpt3(rank, nu-oldnulo+1) 
                     call outpt4(rank) 
                  end if
!----------- end if for the opacity switches ----------------
!------------------------------------------------------------
*
               end do
!
!------------( close big do loop) ---------------------------*
!.... END OF THE FREQUENCY INTEGERATION ---------------------! 
!
#ifdef MPI
! ---- if flux was calculated in parallel, we need to gather it on one core and write it out
         if (freqsplit) then 
          call MPI_Barrier(comm, ierr)
          !    use MPI_allreduce :
          call MPI_allreduce(glsurfin, allsurfin,  nfrcount, 
     &            MPI_double_precision,  MPI_sum, comm, ierr)
      
          call MPI_allreduce(glwave, allwave,  points,
     &            MPI_double_precision,  MPI_sum, comm, ierr)

!-------- first delete all punched files, and mpsa.wave
         
          close (unit = 17, status = 'delete') 
! --- on rank 0 write to final file:
          if (rank .eq.0) then 

            close (unit = 41, status = 'delete')

            open(unit = 41, file='surf_flux.dat',form = 'formatted', 
     &        access = 'sequential')  

! --- if print gt 0 we should write some header! 
            if (ifprnt .gt. 0 ) then 


             write (41, *) 'teff', teff, 'gravity', glog, wlte,
     &                      'title', title
             write(41, *) nmu, 'angles', (angle(j), j = 1, nmu)

            endif 


            do j = 1, points
             write(41, *) allwave(j), (allsurfin(j,nu),nu = 1, nmu) 
            enddo 
            close(unit=41, status= 'keep')
 
          endif 
         endif 

#endif 

!------------------------------------------------------------! 
!.... FINISH OFF THIS MODEL ITERATION
! --- ifsurf is 0 for default and only in model calculations:
!

               if (ifsurf .eq. 0) then 
#ifdef MPI
          if( (mode1 .eq. 1) .and. (sizee .gt. 1)) then 
            call finiteration
          endif 
#endif
! next calls are only made by rank 0 if model is calculated: 
 
                 if (rank .eq. 0) then 
                   print *, 'Final calls to ross '
                   call ross3
                   print *, '	      to high '
                   call high
                   print *, '	  and to radiap '

!
!.... FOLLOWING TWO LINES WERE REVERSED FOLLOWING BOB'S CODE OF 28 JUNE 1995
!
                   call radiap3
                   print *, 'returned from radiap3 for ray  ', imod 
                   if (ifpres .and. ifconv ) call convec
!
                   print *, '     tidy up, tcorr, stateq and output '
                   if (ifcorr) call tcorr3
                   if (nlteon) call stateq3
                   if (ifturb) call turb
                   if (ifprnt .gt. 0) call outpt5(rank) 
                   if (ifprnt .gt. 0) print *, ' We called outpt5'

                  end if 

               end if
!----------------------------------------------------------!
! --------- clean up  punch files -------------------------!

               if(ifpnch(iterat) .gt. 0) then

#ifdef MPI
                  if (rank .eq. 0 ) then

                  if (imod .eq. nmod) then
                     close (unit = 41, status = 'keep')
                  else
                     close (unit = 41, status = 'delete')
                  end if

                  end if

#else
                  if (imod .eq. nmod) then
                     close (unit = 41, status = 'keep')
                  else
                     close (unit = 41, status = 'delete')
                  end if

#endif


               else
                  close (unit = 41, status = 'delete')
               end if
!
            end do
!--END DO (2)---------------------------------------------*
!.... THE END OF LOOP OVER RAYS
           iterat = iterat + 1
#ifdef MPI
      if((sizee .gt. 1) .and. (mode1 .eq. 1)) then
        call MPI_BCAST( iterat, 1, MPI_INTEGER, 0, comm, ierr)
      end if 
#endif 




         end do   ! -- iteration  while loop        
        end if
      end do
!--END DO (1) --------------------------------------------
!
!----- tidy up ------------------------------------------
!--------------------------------------------------------


      if (rank .eq. 0 ) then
      close (unit = 15)
      close (unit = 33)
      end if
!----------------------------------------------------------------
!-  final arrangements if the model calculation was done : 
!-

      if (mode1 .eq. 1) then
        if (rank .eq. 0 ) then 
! ---- last make sure pressure and xne correspond to new model
          call calc_press_hydro
          call pops(0.0d0, 1, xne)        
! ----  write model to the filefluxmodel
          open (unit = 34, file = filefluxmodel, form = 'formatted',
     &       access = 'sequential')
          call outpt6
          do j = 1, nrhox
            hnu(j) = 0.0d0
            jnu(j) = 0.0d0
          end do 
          close (unit = 34)
       
        end if 
     
        call dealloc_comm

        if(allocated(stau)) deallocate(stau)
        if(allocated(st)) deallocate(st)
        if(allocated(sp)) deallocate(sp)

! if we also want to calculate the flux, then close & and delet all previous output
! and go back to the beginning!

        if (flagflux .eq. 1) then         
          mode1 = 2
!    close a couple of files!!!!
          close (unit = 19)
          close (unit = 41, status = 'delete')
          close (unit = 16, status = 'delete')
          close (unit = 17, status = 'delete')
          close (unit = 18, status = 'delete')

          call close_freq
    
          go to 10000

        end if 

      endif


      if(ifprnt .gt. 0 .and. ifprnt .le. 4) then
         close (unit = 16, status = 'keep')
*
      else
         close(unit = 16, status = 'delete')
      end if
*
      if(ifpnch(1) .gt. 0) then
         close (unit = 17, status = 'keep')
*....    close (unit = 41, status = 'keep')
*
      else
         close (unit = 17, status = 'delete')
*....    close (unit = 41, status = 'delete')
      end if
*
      if(ifprnt .eq. 5) then
         close (unit = 18, status = 'keep')
      else
         close (unit = 18, status = 'delete')
      end if
*
      call close_freq
      if (freqsplit) call close_splitf
      end if 
!! end of the flagmodel and flagflux 


#ifdef MPI
      time2 = MPI_WTIME()
      time2 = time2-time1
      if (rank .eq. 0 ) then
       print*, 'time difference: ',time2
      endif

      call MPI_Finalize(ierr)
#endif
      
      print*, 'O.K.'

      call exit(0)  
7000  continue
      print*,'mpsa.control is missing'
      go to 8000
  
7100  continue
      print*, 'can not find  ', INPUTfile
      go to 8000

7200  continue
      print*, 'Can not find ', modelfile
      go to 8000

8000  continue
#ifdef MPI
      print*,'going to abort ..'
      call mpi_abort(comm, ierr)

#else
      print*, 'going to abort ...'
      stop
#endif

      end
*
*********** E N D   O F   T H E  M A I N  P R O G R A M ****************
*
      subroutine deriv (x, f, dfdx, n)
      implicit none
*
*.... ASSUMES THAT ANY ZERO IN X OCCURS AT AN ENDPOINT
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision  dfdx(*), f(*), x(*)
      integer  n
*
!---------------------------- LOCAL VARIABLES --------------------------
*
      double precision  d, d1, s, scale, tan, tan1
      integer j, n1
*
!------------------------------- EXECUTION -----------------------------
*
      dfdx(1) = (f(2) - f(1)) / (x(2) - x(1))
      n1 = n - 1
      dfdx(n) = (f(n) - f(n1)) / (x(n) - x(n1))
*
      if (n .gt. 2) then
         s = abs (x(2) - x(1)) / (x(2) - x(1))
*
         do j = 2, n1
            scale = max(abs(f(j-1)), abs(f(j)), abs(f(j+1))) / 
     &              abs(x(j))
            if (scale .eq. 0.0d0) scale = 1.0d0
            d1 = (f(j+1) - f(j)) / (x(j+1) - x(j)) / scale
            d =  (f(j)   - f(j-1)) / (x(j)   - x(j-1)) / scale
            tan1 = d1 / (s * sqrt(1.0d0 + d1 ** 2) + 1.0d0)
            tan =  d  / (s * sqrt(1.0d0 + d  ** 2) + 1.0d0)
            dfdx(j) = (tan1 + tan) / (1.0d0 - tan1 * tan) * scale
         end do
*
      end if
*
      end
*
*********** E N D   O F   S U B R O U T I N E   D E R I V E ************
*
      double precision function expi (n, x)
      implicit none
*
*.... EXPONENTIAL INTEGRAL FOR POSITIVE ARGUMENTS AFTER CODY AND
*.... THACHER, MATH. OF COMP.,22,641(1968)
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision  x
      integer  n
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision 
     &         a0, a1, a2, a3, a4, a5, 
     &         b0, b1, b2, b3, b4,
     &         c0, c1, c2, c3, c4, c5, c6, 
     &         d1, d2, d3, d4, d5, d6, 
     &         e0, e1, e2, e3, e4, e5, e6, ex, ex1,
     &         f1, f2, f3, f4, f5, f6, 
     &         x1
      integer  i
      save     ex1, x1
*
!--------------------------- INITIALIZATION ---------------------------
*
      data x1 / -1.0d20 /
*
      data a0, a1, a2, a3, a4, a5 /
     &    -44178.5471728217d0,  57721.7247139444d0,  
     &      9938.31388962037d0,  1842.11088668000d0,   
     &       101.093806161906d0,    5.03416184097568d0 /
*
      data b0, b1, b2, b3, b4 /
     &       76537.3323337614d0, 32597.1881290275d0, 
     &        6106.10794245759d0,  635.419418378382d0, 
     &          37.2298352833327d0 /
*
      data c0, c1, c2, c3, c4, c5, c6 /
     &      4.65627107975096d-7, 0.999979577051595d0, 
     &      9.04161556946329d0, 24.3784088791317d0,    
     &     23.0192559391333d0,   6.90522522784444d0,
     &      0.430967839469389d0 /
*
      data d1, d2, d3, d4, d5, d6 /
     &   10.0411643829054d0, 32.4264210695138d0, 41.2807841891424d0,
     &   20.4494785013794d0,  3.31909213593302d0, 0.103400130404874d0/
*
      data e0, e1, e2, e3, e4, e5, e6 /
     &      -0.999999999998447d0, -26.6271060431811d0, 
     &    -241.055827097015d0,   -895.927957772937d0, 
     &   -1298.85688746484d0,    -545.374158883133d0,
     &      -5.66575206533869d0 /
*
      data f1, f2, f3, f4, f5, f6 /
     &      28.6271060422192d0,  292.310039388533d0, 
     &    1332.78537748257d0,   2777.61949509163d0, 
     &    2404.01713225909d0,    631.657483280800d0 /
*
!------------------------------- EXECUTION -----------------------------
*
      if (x .ne. x1) then
         x1 = x
         ex = dexp(-x1)
*
         if (x1 .gt. 4.0d0)  then
            ex1 = (ex + ex * (e0 + (e1 + (e2 + (e3 + (e4 + (e5 + e6 /
     &              x1) / x1) / x1) / x1) / x1) / x1) /
     &       (x1 + f1 +(f2 + (f3 + (f4 + (f5 + f6 / x1) / x1) / x1)
     &               / x1) / x1)) / x1
*
         else if (x1 .gt. 1.0d0) then
            ex1 = ex * (c6 + (c5 + (c4 + (c3 + (c2 + (c1 + c0 * x1) *
     &            x1) * x1) * x1) * x1) * x1) / (d6 + (d5 + (d4 +
     &            (d3 + (d2 + (d1 + x1) * x1) * x1) * x1) * x1) * x1)
*
         else if (x1 .gt. 0.0d0) then
            ex1 = (a0 + (a1 + (a2 + (a3 + (a4 + a5 * x1) * x1) * x1) *
     &               x1) * x1) / (b0 + (b1 + (b2 + (b3 + (b4 + x1) *
     &               x1) * x1) * x1) * x1) -dlog(x1)
*
         else
            ex1 = 0.0d0
         end if
*
      end if
*
      expi = ex1
*
      if (n .gt. 1) then
*
         do i = 1, n - 1
            expi = (ex - x1 * expi) / dble(i)
         end do
*
      end if
*
      end
*
*********** E N D   O F   F U N C T I O N    E X P I *******************
*
      subroutine high
      use types
      use atlcomm
      implicit none
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.abtotb'
      include 'common.height'
      include 'common.rhoxbl'
      include 'common.stateb'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  rhoinv(maxd)
      integer  j
*
!-------------------------------  EXTERNALS ----------------------------
*
      external integ
*
!------------------------------- EXECUTION -----------------------------
*
      do j = 1, nrhox
         rhoinv(j) = 1.0d-5 / rho(j)
      end do
*
      call integ (rhox, rhoinv, height, nrhox, 0.0d0)
      end
*
*********** E N D   O F   F U N C T I O N    H I G H *******************
*
      subroutine integ (x, f, fint, n, start)
      use types
      use atlcomm

      implicit none
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision  f(*), fint(*), start, x(*)
      integer  n
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.parblk'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  db, df, wt
      integer           i, nn
      logical cast_version
      double precision a1(maxd), b1(maxd), c1(maxd)

*
!-------------------------------  EXTERNALS ----------------------------
*
      external parcoe !oparcoe 
*
!------------------------------- EXECUTION -----------------------------
*
!      cast_version = .false.

!      if (cast_version ) then 

!       call oparcoe(f, x, a1, b1, c1, n)
!       fint(1) = start 
!       nn = n -1 
       
!       do i = 1, nn
!        fint(i+1) = fint(i) +(a1(i)+b1(i)/2.*(x(i+1)+x(i))+
!     &     c1(i)/3.*((x(i+1)+x(i))*x(i+1)+x(i)*x(i)))*(x(i+1)-x(i))
!       end do 
!
!      else 

      nn = n
      call parcoe (f, x, (nn))
      fint(1) = start
      df      = x(2) - x(1)
      fint(2) = fint(1) + (a(1,1) + (b(1,1) * 0.5d0 + c(1,1) /
     &          3.0d0 * df) * df) * df
*
      do i = 2, nn - 2
         wt        = 1.0d0
         if (c(i,1) .ne. 0.0d0 .or. 
     &       c(i,2) .ne. 0.0d0 ) wt = abs(c(i,1)) / (abs(c(i,1)) + 
     &      abs(c(i,2)))
         df        = (x(i) - x(i+1))
         db        = -df
         fint(i+1) = fint(i) + wt * (a(i,2) + (b(i,2) * 0.5d0 + c(i,2)/
     &      3.0d0 * db) * db) * db - (1.0d0 - wt) * (a(i,1) + (b(i,1) *
     &      0.5d0 + c(i,1) / 3.0d0 * df) * df) * df
      end do
*
      df       = x(nn) - x(nn - 1)
      fint(nn) = fint(nn - 1) + (a(nn,1)+(b(nn,1)*0.5d0+c(nn,1) /
     &         3.0d0 * df) * df) * df


!      end if

      end
*
*********** E N D   O F   S U B R O U T I N E   I N T E G **************
*
      integer function map1 (xold, fold, nold, xnew, fnew, nnew)
      implicit none
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision  fnew(*), fold(*), xnew(*), xold(*)
      integer  nnew, nold
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  ab, af, bb, bf, cb, cf, d, db, df, wt
      integer  l, ll, k, lm1, lm2, lp1
*
!------------------------------- EXECUTION -----------------------------
*
      l = 2
      ll = 0
*
      do k = 1, nnew
*
         do while (l .le. nold .and. xnew(k) .ge. xold(l))
            l = l + 1
         end do
*
         if (l .gt. nold) l = nold
*
         if (l .gt. 2 .and. l .lt. nold) then
*
*.... PARABOLIC CASE
*
            if (l .ne. ll) then
*
               if (l .gt. 3 .and. l .eq. ll+1) then
                  ab = af
                  bb = bf
                  cb = cf
*
               else
*
*.... MUST COMPUTE THE BACKWARD COEFFICIENTS
*
                  lm1 = l - 1
                  lm2 = l - 2
                  d = (fold(lm1) - fold(lm2)) /
     &                (xold(lm1) - xold(lm2))
                  cb = ((fold(l) - fold(lm1)) /
     &                  (xold(l) - xold(lm1)) - d) /
     &                 (xold(l) - xold(lm2))
                  bb = d + cb * (xold(lm1) - xold(lm2))
                  ab = fold(lm1)
               end if
*
               lp1 = l + 1
               lm1 = l - 1
               d = (fold(l) - fold(lm1)) / (xold(l) - xold(lm1))
               cf = ((fold(lp1) - fold(l)) / 
     &               (xold(lp1) - xold(l)) - d) /
     &              (xold(lp1) - xold(lm1))
               bf = d + cf * (xold(l) - xold(lm1))
               af = fold(l)
               wt = 0.0d0
               if (cf .ne. 0.0d0) wt = abs(cf) / (abs(cf) + abs(cb))
               ll = l
            end if
*
            df = xnew(k) - xold(l)
            db = xnew(k) - xold(lm1)
            fnew(k) = (1.0d0 - wt) * (af + (bf + cf * df) * df) + wt * 
     &                (ab + (bb + cb * db) * db)
*
         else
*
            if (l .ne. ll) then
               ll = l
               lm1 = l - 1
               af = fold(lm1)
               bf = (fold(l) - fold(lm1)) / (xold(l) - xold(lm1))
            end if
*
            fnew(k) = af + bf * (xnew(k) - xold(lm1))
         end if
*
      end do
*
      map1 = ll - 1
      end
*
*********** E N D   O F   F U N C T I O N   M A P 1 ********************
*
      subroutine output(rank) 
      use types
      use atlcomm

      implicit none
*
*.... MODIFIED TO PRINT OUT <JUN(J)> FOR IFPRNT => 1
*.... 1993 JUL - MADE MAXMU A GLOBAL PARAMETER IN COMMON.SIZEBL
*.... 1994 JUL - FIXED BUG WITH DIMENSION AND EQUIVALENCE OF DUMMY
!                AND THE DIMENSION OF ABLOG.  THIS FORCED A CHANGE
!                IN THE WAY THAT THE OPACITIES ARE OUTPUT.
*.... 1995 OCT - EXPLICIT CONVERSION OF LOGICAL ifop TO INTEGER iop
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.abross'
      include 'common.abtotb'
      include 'common.convbl'
      include 'common.depart'
      include 'common.elembl'
      include 'common.fluxbl'
      include 'common.freqbl'
      include 'common.height'
      include 'common.ifblkk'
      include 'common.ifopbl'
      include 'common.ionsbl'
      include 'common.iterbl'
      include 'common.junkbl'
      include 'common.musblk'
      include 'common.opsblk'
      include 'common.optotb'
      include 'common.putblk'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.taushj'
      include 'common.teffbl'
      include 'common.tempbl'
      include 'common.turbpr'
      include 'common.waveyb'
      include 'common.xabund'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      character  a*3, b*3
      double precision  ablog(38), contin, dummy(maxd, 38),
     &                  hlam, hlamlg, hlammg, hnulg, hnumg, hsurf, 
     &                  jbar(maxd), resid, stepwt, surfin(maxmu), 
     &                  tauend, wave
      integer    i, iop(20), iz, j, jtau1, n, nsteps, nu, mu, nfi
      logical    header
      save       contin, header, hsurf, jbar, n, nsteps, nu, 
     &           surfin, wave
      integer    rank
*
!--------------------------- EQUIVALENCE -------------------------------
*
      equivalence (dummy(1,1), aal1(1))
*
!------------------------------ EXECUTION -------------------------------
*
*.... HEADINGS
*
      print*, ' iter ',  iter

      if(ifprnt .gt. 0) then
         header = .false.
!        nu = nulo - 1
         nu = 0
*
         if (ifpnch(iterat) .ge. 2) then
            write (17, 111) 'teff', teff, 'gravity', glog, wlte, 
     &                      'title', title
111         format (a, f8.0, a9, f8.5, 1x, a / a6, a )
            if(ifsurf .eq. 2) then
               write(17, 112) nmu, 'angles',
     &                                       (angle(mu), mu = 1, nmu)
112            format(i3, a7, 10f7.4 / 10x, 10f7.4)
*
            else if (ifsurf .lt. 2) then
               write (17, 114) '#', 'wave', 'freq', 'flux', 
     &                         'cont', 'resid' 
114            format (a8, a10, a12, a14, a13, a13)
            end if
*
         end if
*
         if (ifprnt .eq. 5) write (18, 111) 'teff', teff, 
     &      'gravity', glog, wlte, 'title', title
*
      end if
*
      return
*
!----------------------- INITIALIZE SUMS OVER STEPS --------------------
*
      entry outpt2(rank) 
*
*.... INITIALIZE <JBAR>
*
      do j = 1, nrhox
         jbar(j) = 0.0d0
      end do
*
      hsurf = 0.0d0
      wave  = 2.997925d17 / freq
      nu    = mod(nu, numnu-nulo+1)
      nu    = nu + 1
*
      do mu = 1, nmu
         surfin(mu) = 0.0d0
      end do
*
      n      = 0
      contin = put
      return
*
!---------------- SUM OVER STEPS AND STEP DEPENDENT QUANTITIES ---------
*
      entry outpt3(rank, nfi) 
      n      = n + 1
      nsteps = iput
      stepwt = put
      hsurf  = hsurf + hnu(1) * stepwt
*
*.... SUM UP <JBAR>
*
      do j = 1, nrhox
         jbar(j) = jbar(j) + jnu(j) * stepwt
      end do
*
      do mu = 1, nmu
         surfin(mu) = surfin(mu) + surfi(mu) * stepwt
      end do
      if (freqsplit) glsurfin(nfi ,1:nmu) = surfin(1:nmu)

*
      if (ifprnt .gt. 0) then
*
         if (nsteps .gt. 1) then
*
            if (.not. header) then
               write (16, 101) 'wave(nm)', 'hlambda', 'log h',
     &            'mag', 'frequency', 'hnu', 'log h', 'mag',
     &            'resid', 'tauone taunu'
               header = .true.
            end if
*
            if(n .eq. 1 .and. ifprnt .gt. 1) write (16, 303) '0'
303         format (a)
*
            if(contin .gt. 0.0d0) then
               resid = hnu(1) / contin
*
            else
               resid = 1.0d0
            end if
*
            hnulg = log10(max(hnu(1), 1.0d-50))
            hnumg = -2.5d0 * hnulg
            jtau1 = 1
*
            do while (taunu(jtau1) .lt. 1.0d0)
               jtau1 = jtau1 + 1
            end do
*
            tauend = log10(taunu(nrhox))
            if(ifprnt .gt. 1) write (16, 305) stepwt, hnu(1), 
     &         hnulg, hnumg, resid, jtau1, tauend
305         format (57x, f10.8, 1pe14.4, 0pf12.5, f10.3, f9.5, i5, f8.2)
         end if
*
         if (ifprnt .eq. 3) then
            write (16, 312) 'wavelength(nm)', wave, 'frequency', freq,
     &         'rhox','taunu','abtot','alpha','Bnu','Snu','Jnu','Hnu',
     &         (j, rhox(j), taunu(j), abtot(j), alpha(j), bnu(j),
     &              snu(j), jnu(j), hnu(j), j = 1,nrhox)
312         format ( // 6x, a, f12.3, 3x, a, 1pe13.6 //
     &              12x, a, 7x, 3(a, 6x), 3(a, 8x), a //
     &              (6x, i2, 1p8e11.3))
*
         else if (ifprnt .eq. 4) then
            write (16, 320) 'opacities'
320         format (t55, a, /)
            write (16, 321) 
     &         (' ifop ', i = 1, 20), 
     &         (i, i = 1, 9), (i, i = 10, 20), 
     &         '  H1  ', '  H2+ ', '  H-  ', ' Hray ', '  He1 ', 
     &         '  He2 ', '  He- ', ' Heray', ' Cool ', ' Luke ', 
     &         '  Hot ', ' Elec ', ' H2ray', ' Hline', ' Lines', 
     &         ' Lscat', ' Xline', ' Xlsct', ' Xcont', ' Xscat'
321         format (8x, 20a6 / 5x, 9(5x, i1), 1x, 11(4x, i2) // 
     &              8x, 20a6 /)
*
            do j = 1, nrhox
*
               do i = 1, 38
                  ablog(i) = 0.0d0
                  if(dummy(j,i) .gt. 1.0d-99)ablog(i)=log10(dummy(j,i))
               end do
*
*..... THIS STRANGE ORDER IS NEEDED TO PULL THE CORRECT ENTRIES OUT OF
*..... THE ALPHABETIZED common.opsblk
*
               write (16, 326) j, 
!                   H1         H2+        H-         Hray       He1
     &           ablog(13), ablog(6),  ablog(11), ablog(32), ablog(7),  
!                   He2        He-       Heray       Cool      Luke
     &           ablog(8),  ablog(9),  ablog(34), ablog(3),  ablog(15),  
!                   Hot       Elec       H2ray       Hline     Lines
     &           ablog(12), ablog(31), ablog(33), ablog(10), ablog(14), 
!                  Lscat      Xline      Xlsct       Xcont     Xscat
     &           ablog(35), ablog(18), ablog(37), ablog(17), ablog(36),
     &           j
326            format (i4, t8, 20f6.2, 1x, i3)
            end do
*
         end if
*
      end if
      return
*
!-------------------------- PRINT SUMS OVER STEPS ----------------------
*
      entry outpt4(rank) 
*
*.... OUTPUT <JBAR>
*
      if (ifprnt .eq. 5) then
         write(18, 450) 'frequency', nu, freq, 
     &                        '  wavelength(nm)', wave
450      format( a, i5, 1pe15.6, a, 0pf10.3 /)
         write(18, 451) ('depth', '<Jnu(j)>', j = 1, 4)
451      format(1x, a, 3x, a, 3(4x, a, 3x, a), /)
         write(18, 452) (j, jbar(j), j = 1, nrhox)
452      format((i4, 1pe14.4, 3(i6, e14.4)))
      end if
*
      if(nsteps .eq. 1) contin = hsurf

      if(contin .gt. 0.0d0) resid = hsurf / contin
      if(hsurf .le. 0.0d0) hsurf = 1.0d-50
*
      if(ifprnt .gt. 1) then
         jtau1 = 1
*
         do while (taunu(jtau1) .lt. 1.0d0)
            jtau1 = jtau1 + 1
         end do
*
         tauend = log10 (taunu(nrhox))
*
         if (nsteps .gt. 1) then
            jtau1  = 0
            tauend = 0.0d0
         end if
*
         if (ifsurf .le. 1) then
*
            if (.not. header) then
               header = .true.
               write (16, 101) 'wave(nm)', 'hlambda', 'log h',
     &            'mag', 'frequency', 'hnu', 'log h', 'mag',
     &            'resid', 'tauone taunu'
101            format ( //
     &                 12x, a, 4x, a, 2(7x, a), 4x, a, 8x, a, 10x, a,
     &                 7x, a, 4x, a, 2x, a)
            end if
*
            hlam   = hsurf * freq / wave
            hnulg  = log10 (hsurf)
            hlamlg = log10 (hlam)
            hlammg = -2.5d0 * hlamlg
            hnumg  = -2.5d0 * hnulg
            write (16, 401) nu, wave, hlam, hlamlg, hlammg, freq,
     &         hsurf,hnulg, hnumg, resid, jtau1, tauend, nu
401         format (/, i5, f14.3, 1pe13.4, 0pf12.5, f10.3, 1pe14.6,
     &              e13.4, 0pf12.5, f10.3, f9.5, i5, f8.2, i5)
*
         else if (ifsurf .eq. 2) then
*
            if( .not. header) then
               header = .true.
               write(6,102) 'wave(nm)', 'frequency', 
     &            'tauone taunu', 'mu  intensity', 'mu  intensity', 
     &            'mu  intensity', 'mu  intensity', 'mu  intensity'
102            format( /////10x, a, 5x, a, 3x, a, 5(4x, a))
            end if
*
            write (16, 406) nu, wave, freq, jtau1, tauend, 
     &                     (angle(mu), surfin(mu), mu = 1, nmu)
406         format (/, i5, f10.3, 1pe15.6, i6, 0pf6.2, 
     &              5(0pf7.4, 1pe10.3) / (42x, 5(0pf7.4, 1pe10.3)))
         end if
*
      end if
*
      if(ifpnch(iterat) .ge. 2) then
*
         if (ifsurf .eq. 2) then

          if (rank .eq. 0 ) then
            write (41, 441) wave
441         format (f12.5)
          end if


            write (17, 416) (surfin(mu), mu = 1, nmu)
416         format ((1p8e11.3))

         else if (ifsurf .lt. 2) then
            write (17, 411) 'flux', nu, wave, freq, hsurf, contin, 
     &                       resid
411         format (a, i5, f10.2, 1pe16.6, e13.4, e13.4, 0pf10.5)
            if(nu .eq. nuhi) write (17, 411) 'flux'
*
         end if
*
      end if
      return
*
!-------------------------------- SUMMARIES ----------------------------
*
      entry outpt5(rank) 
*
      if (ifprnt .gt. 0) then
          write (16, *) 'DEBUG(OUT 1) ' 
         write (16, 501) 'rhox', 'ptotal', 'pturb', 'grdadb',
     &      'dltdlp', 'velsnd', 'dlrdlt', 'heatcp', 'hscale', 'vconv',
     &      'flxcnv',
     &      (j, rhox(j), ptotal(j), pturb(j), grdadb(j),
     &           dltdlp(j), velsnd(j), dlrdlt(j), heatcp(j),
     &           hscale(j), vconv(j), flxcnv(j), j = 1, nrhox)
501       format (/ 8x,a,7x,a,5x,a,6x,a,5x,a,4x,a,4(5x, a),6x,a,//
     &          (i3, 1p11e11.3))


!----------------------------------------------------------------
!          write (16, 501) 'rhox', 'ptotal', 'pturb', 'grdadb',
!     &      'dltdlp', 'velsnd', 'dlrdlt', 'heatcp', 'hscale', 'vconv',
!     &      'flxcnv'
!501          format( 11(a, 7x))
!
!          write (16, 5001) (j, rhox(j), ptotal(j), pturb(j), grdadb(j),
!     &           dltdlp(j), velsnd(j), dlrdlt(j), heatcp(j), 
!     &           hscale(j), vconv(j), flxcnv(j), j = 1, nrhox)
!5001     format (i3, 11 (1p11e11.3) //)
!
!-----------------------------------------------------------------

         write (16, *) ' DEBUG( OUT 2)'
         write (16, 502) 'flux', flux
502      format (a, 1pe12.4)
         write (16, 503) 'xnatom', 'raden', 'pradk', 'xnfph1',
     &      'xnfph2', 'xnfphe1', 'xnfphe2', 'xnfphe3', 'vturb', 
     &      (j, xnatom(j), raden(j), pradk(j), xnfph(j,1), xnfph(j,2), 
     &          xnfphe(j,1), xnfphe(j,2), xnfphe(j,3), vturb(j),
     &          flxcnv0(j), flxcnv1(j), j = 1, nrhox)
503      format ( // 4(6x, a), 5x, a, 3(4x, a), 5x, a, //
     &      (i3, 1p11e11.3))
         write (16, 99) 'pradk0', pradk0
99       format (a, 1pe12.4)
         write (16, 541) 'teff', teff, 'log g', glog, title,
     &      'ray ', iter
541      format ( // a5, f8.0, a8, f9.2, 8x, a74, 2x, a9, i3 //)
*
         if (.not. ifcorr) then
*
            do j = 1, nrhox
               flxrad(j) = flux - flxcnv(j)
            end do
*
         end if
*
         do j = 1, nrhox
            flxcnv(j) = flxcnv(j) / (flxcnv(j) + flxrad(j))
         end do
*
         write (16, 542) 
     &      'electron', 'rosseland', 'height', 'rosseland', 'fraction', 
     &                  'radiative', 'per cent',  'flux',
     &      'rhox', 'temp', 'pressure', 'number', 'density', 'mean',
     &              '(km)', 'depth', 'conv flux', 'acceleration', 
     &              'error', 'deriv',
     &      (j, rhox(j), t(j), p(j), xne(j), rho(j), abross(j), 
     &          height(j), tauros(j), flxcnv(j), accrad(j), flxerr(j), 
     &          flxdrv(j), j = 1, nrhox)
542   format (34x, a, 11x, a, 3x, a, 2x, a, 2x, a, 2x, a, 3x, a, 3x, a /
     &        7x, a,  6x, a, 3x, a, 3x, a, 3x, a, 4x, a, 7x, a, 5x, a, 
     &         3x, a,  1x, a, 3x, a, 3x, a //
     &       (i3, 1pe10.3, 0pf9.1, 1p7e10.3, e11.3, 0pf10.2, f8.2))
      end if
*
      if(ifpnch(iterat) .gt. 0) then
*
!----------------------------------- PUNCHOUT --------------------------
*
         if (ifconv) then
            a = 'on'
         else
            a = 'off'
         end if
*
         if (ifturb) then
            b = 'on'
         else
            b = 'off'
         end if
*
         do i = 1, 20
            iop(i) = 0
            if(ifop(i)) iop(i) = 1
         end do
*
         write (17, 552) 'teff ', teff, 'gravity', glog, wlte,
     &      'title ', title,
     &      ' opacity ifop', iop,
     &      ' convection ', a, mixlth, ' turbulence ', b, 
     &      trbfdg, trbpow, trbsnd, trbcon,
     &      'abundance scale ', xscale, ' abundance change', 
     &                                  (iz, abund(iz), iz = 1, 2),
     &      (' abundance change',
     &      (iz, abund(iz), iz = i, i+5), i = 3, 98, 6),
     &       ' abundance change', (iz, abund(iz), iz = 99, 99)
552      format (a, f7.0, 2x, a, f8.5, 1x, a /
     &         a, a /
     &         a, 20i2 /
     &         a, a, f6.2, a, a, 4f6.2 /
     &         a, f9.5, a, 2(i2, f8.5) /
     &        (a, 6(i3, f7.2)))
*
*.... CHANGED 20 JULY 1995 TO CONFORM WITH BOB'S mpsav.for;83 OF 20 JUNE 1995
*
         write (17, 554) 'read deck6', nrhox,
     &      ' rhox, t, p, xne, abross, accrad, vturb',
     &      (rhox(j), t(j), p(j), xne(j), abross(j), accrad(j), 
     &       vturb(j), flxcnv(j), vconv(j), j = 1, nrhox)
554      format (a, i3, a / (1pe15.8, 0pf9.1, 1p7e10.3))
         write (17, 555) 'pradk', pradk0
555      format (a, 1pe11.4)
         if (nlteon) write (17, 556) 
     &      'read departure coefficients', nrhox,
     &      ' rhox  bhyd 1-6  bmin',
     &      (rhox(j), (bhyd(j,i), i = 1, 6), bmin(j) , j = 1, nrhox)
556      format (a, i3, a, /
     &          (1pe11.4, 0p7f9.4))
*
         if (.not. ifwave) then
            write (17, 562) 'read frequencies', numnu,nulo,nuhi,freqid
562         format (a, 3i4, 3x, a)
            write (17, 563) (nu, freset(nu), rcoset(nu), nu = 1, numnu)
563         format (i5, 1p2e17.8, i5, 2e17.8)
         end if
         write (17, 571) 'begin',  'ray ', iter, ' completed '
571      format (a, 20x, a, i3, a)
      end if
*
      return
!----------- new output model for further flux-calculations------------!

!- write model to the filefluxmodel
       entry outpt6

       write(34,*) 1
       write(34,*) nrhox
       do j = 1, nrhox
          write(34,340) rhox(j), t(j), p(j), xne(j),
     &            abross(j), accrad(j), vturb(j)
340       format(1PE14.7,0PF9.1,1P8E11.3,0PF12.3,F9.3)

       end do 

!--- write information on model calculation:
         if (ifconv) then
            a = 'on'
         else
            a = 'off'
         end if
!
         if (ifturb) then
            b = 'on'
         else
            b = 'off'
         end if
!
         do i = 1, 20
            iop(i) = 0
            if(ifop(i)) iop(i) = 1
         end do

      write(34,342) 'teff ', teff, 'gravity', glog, wlte,
     &      ' overshoot =  ', overwt ,
     &      ' opacity ifop', iop,
     &      ' convection ', a, mixlth, ' turbulence ', b,
     &      trbfdg, trbpow, trbsnd, trbcon,
     &      'abundance scale ', xscale, ' abundance change',
     &                                  (iz, abund(iz), iz = 1, 2),
     &      (' abundance change',
     &      (iz, abund(iz), iz = i, i+5), i = 3, 98, 6),
     &       ' abundance change', (iz, abund(iz), iz = 99, 99)
342      format (a, f7.0, 2x, a, f8.5, 1x, a /
     &         a, f6.2 /
     &         a, 20i2 /
     &         a, a, f6.2, a, a, 4f6.2 /
     &         a, f9.5, a, 2(i2, f8.5) /
     &        (a, 6(i3, f7.2)))


      return



      end
*
*********** E N D   O F   S U B R O U T I N E   O U T P U T ************
*
      subroutine parcoe (f, x, n)

      
      use types
      use atlcomm


      implicit none
*
*.... NOW ALWAYS IN DOUBLE PRECISION
*.... THIS HAS BEEN REWRITTEN IN THE SPIRIT OF THE U LONDON SUGGESTION
*.....  VAR(M,1 OR 2)  1=FORWARD INTERPOLATION, 2=BACKWARD INTERPOLATION
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision  f(*), x(*)
      integer           n
*
!----------------------------- COMMONS ---------------------------------
*
      include 'common.parblk'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  d, dx, wt 
      integer           j, n1, jp2, jp1, jm1
*
!------------------------------- EXECUTION -----------------------------
*

      n1      = n - 1
      c(1, 1) = 0.0d0
      c(1, 2) = c(1, 1)

      b(1, 1) = (f(2) - f(1)) / (x(2) - x(1))
      b(1, 2) = b(1, 1)

      a(1, 1) = f(1)
      a(1, 2) = a(1, 1)

      c(n, 1) = 0.0d0
      c(n, 2) = c(n, 1)

      b(n, 1) = (f(n) - f(n1)) / (x(n) - x(n1))
      b(n, 2) = b(n, 1)

      a(n, 1) = f(n1)
      a(n, 2) = a(n, 1)
*
!--------------------
      if (n .gt. 2) then
*
         do j = 2, n1 - 1
            jp2     = j + 2
            jp1     = j + 1
            jm1     = j - 1
            d       = (f(jp1) - f(j)) / (x(jp1) - x(j))
            dx      = x(jp2) - x(j)
            c(j, 1) = (f(jp2) - f(jp1)) / (x(jp2) - x(jp1)) /
     &               dx - d / dx
            b(j, 1) = d + c(j,1) * (x(jp1) - x(j))
            a(j, 1) = f(jp1)
*
            if (j .eq. 2) then
*.... MUST COMPUTE THE BACKWARD COEFFICIENT WHEN J = 2
               d       = (f(2) - f(1)) / (x(2) - x(1))
               dx      = x(3) - x(1)
               c(j, 2) = (f(3) - f(2)) / (x(3) - x(2)) / dx - d / dx
               b(j, 2) = d + c(j,2) * (x(2) - x(1))
               a(j, 2) = f(2)
*
            else
               a(j, 2) = a(jm1, 1)
               b(j, 2) = b(jm1, 1)
               c(j, 2) = c(jm1, 1)
            end if
*
         end do
*
         a(n1, 1) = f(n)
         b(n1, 1) = -d
         c(n1, 1) = 0.0d0
         a(n1, 2) = f(n1)
         b(n1, 2) = b(n1-1, 1)
         c(n1, 2) = c(n1-1, 1)
      end if
*
      end
*
*********** E N D   O F   S U B R O U T I N E   P A R C O E ************
*
      subroutine radiap
    
      use types
      use atlcomm
      implicit none
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*.... ON SUBSEQUENT ENTRIES
*
      double precision  rcowt
*
!----------------------------- COMMONS -------------------------------
*
!      include 'common.sizebl'
      include 'common.constb'
*
      include 'common.abtotb'
      include 'common.fluxbl'
      include 'common.rhoxbl'
      include 'common.taushj'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      integer j 
!      double precision  hflux(maxd)
!      save  hflux
*
!-------------------------------  EXTERNALS ----------------------------
*
      external integ
*
!--------------------------- EXECUTION --------------------------------
*
*.... INITIALIZATION
*
      do j = 1, nrhox
         hflux(j)  = 0.0d0
         raden(j)  = 0.0d0
         accrad(j) = 0.0d0
      end do
*
      pradk0 = 0.0d0
      return
*
!---------------------- INTEGRATION OVER FREQUENCY ---------------------
*
      entry radiap2 (rcowt)
*
      do j = 1, nrhox
         accrad(j) = accrad(j) + abtot(j) * hnu(j) * rcowt
         hflux(j)  = hflux(j) + hnu(j) * rcowt
         raden(j)  = raden(j) + jnu(j) * rcowt
        if (hnu(j) .lt. 0 .or. jnu(j) .lt. 0) then
       !   print *, ' OUPS ', hnu(j), jnu(j), raden(j), rcowt
        end if
      end do
*
      pradk0 = pradk0 + knu(1) * rcowt
      return
*
!---------------------------- FINISH IT OFF ----------------------------
*
      entry radiap3
*
        print *, ' Entry 3 radiap '

      do j = 1, nrhox
         raden(j)  = raden(j)  * fourpi / c_cm
         accrad(j) = accrad(j) * fourpi / c_cm
*
*.... FUDGE TO KEEP MODEL FROM BLOWING UP WITH LARGE FLUX ERRORS
*
         if (hflux(j) .gt. flux) accrad(j) = accrad(j) * flux/hflux(j)
      end do
*

      pradk0 = pradk0 * fourpi / c_cm
      if (hflux(1) .gt. flux) pradk0 = pradk0 * flux / hflux(1)

      call integ (rhox, accrad, prad, (nrhox), (accrad(1) * rhox(1)))
        
*
      print *, ' and calculating new radiation pressure '
      do j = 1, nrhox
         pradk(j) = prad(j) + pradk0
      end do
      print *, ' returning home '
!       print *, ' finish pradk0  and pradk', pradk0, pradk
      return
*
      end
*
*********** E N D   O F   S U B R O U T I N E   R A D I A P ************
*
      subroutine ross

      use types
      use atlcomm

      implicit none
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*.... ON SUBSEQUENT ENTRIES
*
      double precision  rcowt
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.constb'
*
      include 'common.abross'
      include 'common.abtotb'
      include 'common.freqbl'
      include 'common.rhoxbl'
      include 'common.tempbl'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      integer  j
*
!-------------------------------  EXTERNALS ----------------------------
*
      external integ
*
!---------------------------- EXECUTION -------------------------------
*
!--------------------------- INITIALIZATION ---------------------------
*
      do j = 1, nrhox
         abross(j) = 0.0d0
      end do
*
      return
*
!----------------------- FREQUENCY INTEGRATION -------------------------
*
      entry ross2 (rcowt)
*
      do j = 1, nrhox
         abross(j) = abross(j) + dbnudt(j) / abtot(j) * rcowt
      end do
*
      return
*
!-------------------------------- COMPLETION ---------------------------
*
      entry ross3
*
      do j = 1, nrhox
         abross(j) = (4.0d0 * sigma / pi) * t(j) ** 3 / abross(j)
      end do
*
      call integ (rhox, abross, tauros, (nrhox), (abross(1)*rhox(1)))
      end
*
*********** E N D   O F   S U B R O U T I N E   R O S S ****************
*
*
      subroutine solvit (aa, ntrue, n, b)
      use types
      use atlcomm

      implicit none
*
*.... THIS REPLACES THE ROUTINE BOB HAD WHICH WAS BASED ON A
*.... GAUSS-JORDAN ELIMINATION SCHEME
*.... THIS SUBROUTINE INVERTS THE MATRIX A(N,N) AND RETURNS
*.... THE SOLUTION IN THE VECTOR B(N).  IT USES THE METHOD OF
*.... TRIANGULAR DECOMPOSITION DISCUSSED IN CHAPTER 9 OF
*.... "A FIRST COURSE IN NUMERICAL ANALYSIS" BY ANTHONY RALSTON.
*.... THE EQUATIONS HAVE THE FORM AX=B.
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*
      integer, intent(in) ::  ntrue, n
      double precision, intent(in) ::  aa(1:ntrue,1:ntrue)
      double precision, intent(out):: b(*)
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  a(maxd, maxd), dd, stor(maxd)
      integer  l1, i, j, jj, pp, p(maxd)
*
!------------------------------- EXECUTION -----------------------------
*
      l1 = n
*
      do i = 1, l1
*
         do j = 1, l1
            a(j,i) = aa(j,i)
         end do
*
      end do
*
      do j = 1, l1
*
         do i = 1, l1
            stor(i) = a(i,j)
         end do
*
         if (j .gt. 1) then
*
            do jj = 1, j - 1
               pp = p(jj)
               a(jj,j) = stor(pp)
               stor(pp) = stor(jj)
*
               do i = jj + 1, l1
                  stor(i) = stor(i) - a(i,jj) * a(jj,j)
               end do
*
            end do
*
         end if
*
         dd = stor(j)
*
         do i = j, l1
*
            if (dabs(dd) .le. dabs(stor(i))) then
               p(j) = i
               a(j,j) = stor(i)
               dd = stor(i)
            end if
*
         end do
*
         if (j .ne. l1) then
            pp = p(j)
            stor(pp) = stor(j)
*
            do i = j + 1, l1
               a(i,j) = stor(i) / a(j,j)
            end do
*
         end if
*
      end do
*
      do i = 1, l1
         stor(i) = b(i)
      end do
*
      do i = 1, l1
         pp = p(i)
         b(i) = stor(pp)
         stor(pp) = stor(i)
*
         if (i .ne. l1) then
*
            do j = i + 1, l1
               stor(j) = stor(j) - b(i) * a(j,i)
            end do
*
         end if
*
      end do
*
      do i = l1, 1, -1
         b(i) = b(i) / a(i,i)
*
         if (i .gt. 1) then
*
            do j = 1, i - 1
               b(j) = b(j) - b(i) * a(j,i)
            end do
*
         end if
*
      end do
*
      end
*
*********** E N D   O F   S U B R O U T I N E   S O L V I T ************
*
      subroutine stateq
      use types
      use atlcomm
      use kappa_cal 

      implicit none
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*.... ON SUBSEQUENT ENTRIES
*
      double precision  rcowt

!      THE BOUND-BOUND COLLISION RATES WERE DERIVED FROM AN ANALYTIC FIT
!      TO THE CROSS SECTION CALCULATIONS OF BURKE, ORMONDE, AND WHITAKER,
!      PROC. PHYS. SOC., 1968, VOL 92, 319
*
!      THE CROSS SECTION USED (IN UNITS OF PI*A0**2) IS
*
!      QIJ = 4*FIJ*(EH/E0)**2*(LOG(E/E0)/(E/E0)+.148 /(E/E0)**6)
*
!      FIJ = OSCILLATOR STRENGTH
!      EH = GROUND STATE BINDING ENERGY
!      E0 = THRESHOLD ENERGY
!                                                   D M PETERSON MAY 1968
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.constb'
*
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.iterbl'
!      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.taushj'
      include 'common.tempbl'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  a(6,6),  dqrad(maxd,6), dqrd(maxd), dt, 
     &                  expi, f(8,8), fre5, gik, hcont(6), hminbf, hvc, 
     &                  q, qassoc, qcharg, qcoll(8,8), qelect, 
     &                  qradik(maxd,6), qradki(maxd,6), qrdhmk(maxd), 
     &                  qrdkhm(maxd), sqrtt, rfrwt, right(6), rj, rje, 
     &                  rjedt, th, theta, told(maxd), x0, y, z
      integer  i, j, l, n
      save     told, qrdhmk, qrdkhm, dqrd, dqrad, qradki, qradik
*
!-------------------------------  EXTERNALS ----------------------------
*
      external solvit
*
!--------------------------- EQUIVALENCE -------------------------------
*.... TURN OFF FOR NOW
!      equivalence (a(1), taunu(1)), (qcoll(1), hnu(1))
*
!--------------------------- INITIALIZATION ---------------------------
*
      data f / 8*0.0d0, 0.4162d0, 7*0.0d0, 0.07910d0, 0.6408d0, 6*0.0d0, 
     &         0.02899d0, 0.1193d0, 0.8420d0, 5*0.0d0, 0.01394d0, 
     &         0.04467d0, 0.1506d0, 1.038d0,  4*0.0d0, 0.007800d0, 
     &         0.02209d0, 0.05585d0, 0.1794d0, 1.231d0, 3*0.0d0, 
     &         0.004814d0, 0.01271d0, 0.02768d0, 0.06551d0, 0.2070d0, 
     &         1.425d0, 2*0.0d0, 0.003184d0, 0.008037d0, 0.01604d0, 
     &         0.03229d0, 0.07455d0, 0.2340d0, 1.615d0, 0.0d0 /
*
!-------------------------------- EXECUTION ----------------------------
*
*.... ERASE FREQUENCY INTEGRALS
*
      do j = 1, nrhox
         told(j)   = t(j)
         qrdhmk(j) = 0.0d0
         qrdkhm(j) = 0.0d0
         dqrd(j)   = 0.0d0
      end do
*
      do i = 1, 6
*
         do j = 1, nrhox
            dqrad(j,i)  = 0.0d0
            qradki(j,i) = 0.0d0
            qradik(j,i) = 0.0d0
         end do
*
      end do
*
      return
*
!------------------------ FREQUENCY INTEGRALS --------------------------
*
      entry stateq2(rcowt)
      rfrwt = fourpi / h * rcowt / freq
      hvc = 2.0 * h * freq * (freq / c_cm) ** 2
*
      do n = 2, 6
         hcont(n) = coulx ((n), (freq), 1.0d0)
      end do
*
      fre5 = freq * 1.0d-5
*
      if (freq .ge. 2.111d14) then
         hminbf = 6.801d-20 + (5.358d-8 + (1.481d3 +
     &      (-5.519d12 + 4.808d21 / fre5) / fre5) / fre5) / fre5
*
      else if (freq .gt. 1.8259d14) then
         hminbf = 3.695d-16 + (-1.251d-1 + 1.052d13 / freq) / freq
*
      else
         hminbf = 0.0d0
      end if
*
      do j = 1, nrhox
         rj    = rfrwt * jnu(j)
         rje   = rfrwt * ehvkt(j) * (jnu(j) + hvc)
         rjedt = rje   * hkt(j)   * freq / t(j)
*
         do i = 2, 6
            qradik(j,i) = qradik(j,i) + hcont(i) * rj
            dqrad(j,i)  = dqrad(j,i)  + hcont(i) * rjedt
            qradki(j,i) = qradki(j,i) + hcont(i) * rje
         end do
*
         qrdhmk(j) = qrdhmk(j) + hminbf * rj
         dqrd(j) = dqrd(j) + hminbf * rjedt
         qrdkhm(j) = qrdkhm(j) + hminbf * rje
      end do
*
      return
*
!---------------------------------- THIRD ENTRY ------------------------
*
      entry stateq3
      if (ifprnt .gt. 0) write (16,201)
     &   'hminus statistical equilibrium', 'rhox', 'qelect', 'qassoc',
     &   'qcharg', 'qrdkhm', 'qrdhmk', 'bmin'
201    format ( ////// 36x, a, / 10x, a, 7x, a, 4(6x, a), 7x, a)
*
      do j = 1, nrhox
         dt        = t(j) - told(j)
         theta     = 5040.0d0 / t(j)
         qelect    = 10.0d0 ** (-8.7) * theta ** (1.5) * xne(j)
         qassoc    = 10.0d0 ** (-8.7) * 2.0 * bhyd(j,1) * xnfph(j,1)
         qcharg    = 10.0d0 ** (-7.4) * theta ** 0.333333 * xnfph(j,2)
         qrdkhm(j) = qrdkhm(j) + dqrd(j) * dt
         bmin(j)   = (qrdkhm(j) + qelect + qassoc + qcharg) /
     &               (qrdhmk(j) + qelect + qassoc + qcharg)
         write (16, 211) j, rhox(j), qelect, qassoc, qcharg, 
     &                      qrdkhm(j), qrdhmk(j), bmin(j)
211      format (i5, 1p6e12.3, 0pf10.4)
      end do
*
      if (ifprnt .gt. 1) write (16, 31) 
     &   'statistical equilibrium rates',
     &   'rate=sign(log10(max(abs(rate*1.0e20),1.0)),rate) ',
     &   '  rad   1-k   k-1   2-k   k-2   3-k   k-3',
     &        '   4-k   k-4   5-k   k-5   6-k   k-6',
     &   '  coll  1-k   2-k   3-k   4-k   5-k   6-k   5-8   6-8  ',
     &   '  coll  1-2   1-3   1-4   1-5   1-6   1-7   2-3   2-4  ',
     &           '2-5   2-6   2-7   3-4   3-5   3-6   3-7   4-5 ',
     &         '  4-6   4-7   5-6   5-7   6-7 '
31    format ( ////// 30x, a, 4x, a, /
     &        a, a, a, a, /
     &        a, a, a)
*
      do j = 1, nrhox
         dt    = t(j) - told(j)
         th    = hydip / tkev(j)
         sqrtt = sqrt (t(j))
*
         do i = 1, 7
            y = i
            qcoll(i,i) = 2.2d-8 * y**3 / sqrt(th) * exp(-th / y**2) *
     &                   xne(j)
*
*.... QCOLL(I,I) IS THE BOUND FREE RATE
*
            do l = i + 1, 8
               z = l
               gik = 1.0d0 / y**2 - 1.0d0 / z**2
               x0 = th * gik
               q = 2.186d-10 * f(i, l) / gik**2 * x0 * sqrtt *
     &                  (expi(1, (x0)) + 0.148d0 * x0 * expi(5, (x0)))
               qcoll(i, l) = q * xne(j)
               qcoll(l, i) = qcoll(i, l) * (y/z)** 2 * exp(x0)
            end do
*
         end do
*
         y = 8.0
         qcoll(8,8) = 2.2d-8 * y**3 / sqrt(th) * exp(-th / y**2) *
     &                xne(j)
*
         do i = 1, 6
            a(i,i) = qradik(j,i)
            qradki(j,i) = qradki(j,i) + dqrad(j,i) * dt
            right(i) = qradki(j,i) + qcoll(i,i) + qcoll(i,7) +
     &                 qcoll(i,8)
*
            do l = 1, 8
               a(i, i) = a(i, i) + qcoll(i, l)
            end do
*
            if (i .lt. 6) then
*
               do l = i + 1, 6
                  a(i, l) = -qcoll(i, l)
                  a(l, i) = -qcoll(l, i)
               end do
*
            end if
*
         end do
*
         call solvit (a, (6), (6), right)
*
         do i = 1, 6
            bhyd(j,i) = right(i)
         end do
*
         if (ifprnt .gt. 1) then
*
            do i = 1, 6
               qradki(j,i) = sign (log10 (
     &            max(abs(qradki(j,i) * 1.0d20),1.0d0)) , qradki(j,i))
               qradik(j,i) = sign (log10 (
     &            max(abs(qradik(j,i)*1.0d20),1.0d0)) , qradik(j,i))
            end do
*
            do i = 1, 8
*
               do l = 1, 8
                  qcoll(i, l) = sign (log10 (max(
     &               abs(qcoll(i, l) * 1.0d20),1.0d0)) , qcoll(i, l))
               end do
*
            end do
*
            write (16, 100) '0', j, (qradik(j,i), qradki(j,i), i = 1,6),
     &         (qcoll(i,i), i = 1, 6), qcoll(5,8), qcoll(6,8)
100         format (a, i5, 12f6.2, 6x, 8f6.2)
            write (16, 110) 
     &         (qcoll(1,l), l = 2, 7), (qcoll(2,l), l = 3, 7),
     &         (qcoll(3,l), l = 4, 7), (qcoll(4,l), l = 5, 7),
     &         (qcoll(5,l), l = 6, 7),  qcoll(6,7)
110       format (6x, 21f6.2)
         end if
*
      end do
*
      write (16, 170) 'statistical equilibrium for hydrogen',
     &   'rhox', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6',
     &   (j, rhox(j), (bhyd(j,i), i = 1, 6), j = 1, nrhox)
170   format ( ////// 30x, a /
     &       15x, a, 10x, a, 5(8x, a) /
     &       (8x, i2, 1pe11.4 , 1x, 0p6f10.4))
      end
*
*********** E N D   O F   S U B R O U T I N E   S T A T E Q ************
*
      subroutine tcorr

      use types
      use atlcomm
      implicit none
*
*.... 10 JUNE 1993 VERSION
*.... 1995 JUL 20 - CHANGED TO CONFORM TO BOB'S mpsav.for;83 BY
!                   ADDING LOCAL VARIABLE cccccc, AND A LOOP THAT USES IT
*
!--------------------------- DUMMY ARGUMENTS ---------------------------
*.... ON SUBSEQUENT ENTRIES
*
      double precision  rcowt
*
!----------------------------- COMMONS ---------------------------------
*
      include 'common.constb'
*
      include 'common.abross'
      include 'common.abtotb'
      include 'common.convbl'
      include 'common.depart'
      include 'common.fluxbl'
      include 'common.freqbl'
      include 'common.iterbl'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.steplg'
      include 'common.taushj'
      include 'common.tcorrb'
      include 'common.teffbl'
      include 'common.tempbl'
      include 'common.turbpr'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision ab1(maxd), ab2(maxd), cccccc(maxd), cnvfl, 
     &                 cnvflx(maxd), codrhx(maxd), d, dabros(maxd), 
     &                 dabtot(maxd), ddel(maxd), ddlt(maxd), del,
     &                 diagj, drhox(maxd), dtau(maxd), dtdrhx(maxd), 
     &                 dtflux(maxd), dtlamb(maxd), dtsur, dtsurf(maxd), 
     &                 dum(maxd), dum1(maxd), dum2(maxd), dum3(maxd),
     &                 dum4(maxd), dum5(maxd), dum6(maxd), dum7(maxd),
     &                 dum8(maxd), dum9(maxd),ex, expi, fluxco, g(maxd),
     &                 gflux(maxd), hratio(maxd), dum10(maxd), 
     &                 hratio1, oldt1(maxd), p1(maxd), 
     &                 p2(maxd), ppp(maxd), prdnew(maxd), ptot1(maxd), 
     &                 ptot2(maxd), ratioA, 
     &                 rdiagjnu(maxd), rhoxstd, rjminsnu(maxd), 
     &                 rrr(maxd), t1(maxd),
     &                 taustd(maxd), tav, teff25, term1, 
     &                 term2, tinteg(maxd), tnew1(maxd), tnew2(maxd), 
     &                 tone(1), tplus(maxd), ttwo(1), vco , dumm1(1) 
      integer          i, idum, j, map1, ifudge 
      integer          j1smooth, j2smooth

*
      
      double precision pturbnew(maxd), tsmooth(maxd) 
      double precision taub
      double precision, parameter :: errorT = 1.0d-05 
      logical smoothing 
      save     rdiagjnu, rjminsnu

*
!-------------------------------  EXTERNALS ----------------------------
*
      external deriv, integ, map1
*
!----------------------------- INITIALIZATION -------------------------
*
*.... VARIABLES FOR THE MIXING LENGTH FROM HENYEY ET AL
*.... THESE ARE NOT USED HERE FOR NOW.
*
      double precision  nuc, yc
      data  yc, nuc / 0.5d0, 8.0d0 /

! ---- smoothing weigths
      double precision wtj(3)
      data (wtj(i), i=1,3) / 0.33d0, 0.34d0, 0.33d0 / 
*
!------------------------------- EXECUTION -----------------------------
      smoothing = .false. 
 


*
*.... ERASE FREQUENCY INTEGRALS ON THE FIRST ENTRY
*
      do j = 1, nrhox
         flxrad(j) = 0.0d0
         rdabh(j)  = 0.0d0
         rdiagj(j) = 0.0d0
         rjmins(j) = 0.0d0
      end do

*
      return
*
!------------------ THE SECOND ENTRY = FREQUENCY INTEGRATION -----------
*
      entry tcorr2 (rcowt)
      call deriv (rhox, abtot, dabtot, (nrhox))
*
      do j = 1, nrhox
         flxrad(j)   = flxrad(j) + rcowt * hnu(j)
         rdabh(j)    = rdabh(j)  + rcowt * hnu(j) * dabtot(j) / 
     &                 abtot(j)
         rjminsnu(j) = abtot(j) * rcowt * jmins(j)
         rjmins(j)   = rjmins(j) + rjminsnu(j)
      end do
*
      term2 = 0.0d0
*
      do j = 1, nrhox
         term1 = term2
         if(j .ne. nrhox) dtaunu(j) = taunu(j+1) - taunu(j) 
         d = max(dtaunu(j), 1.0d-10)
*
         if (d .gt. 0.01d0) then
            ex = 0.0d0
            if (d .lt. 10.0d0) ex = expi (3, (d))
            term2 = 0.5d0 * (d + ex - 0.5d0) / d
*
         else
            term2 = (0.922784335098467d0 - log(d)) * d / 4.0d0 + 
     &              d ** 2 / 12.0d0 - d ** 3 / 96.0d0 + d ** 4 / 
     &              720.0d0
         end if
*
         diagj = term1 + term2
         rdiagjnu(j) = abtot(j) * (diagj - 1.0d0) / 
     &                 (1.0d0 - alpha(j) * diagj) * 
     &                 (1.0d0 - alpha(j)) * dbnudt(j) * rcowt
         rdiagj(j) = rdiagj(j) + rdiagjnu(j)
      end do
*
      return
*
!------------------------------- ENTRY 3 -------------------------------
*
*.... AVRETT-KROOK TAU CORRECTION MODIFIED FOR CONVECTION
*
      entry tcorr3
      call deriv (rhox, t,      dtdrhx, (nrhox))
      call deriv (rhox, dltdlp, ddlt,   (nrhox))
      call deriv (rhox, abross, dabros, (nrhox))
*
      do j = 1, nrhox
         cnvflx(j) = 0.0d0
      end do
*
      if (ifconv) then
*
         do j = 3, nrhox
            cnvflx(j) = flxcnv(j)
         end do
*
      end if
*
*.... MODIFIED 20 JULY 1995
*.... CONFORMS WITH BOB'S mpsav.for;83 OF 28 JUNE 1995
*
      do j = 2, nrhox-1
         cccccc(j) = 0.25d0 * (cnvflx(j-1) + 2.0d0 * cnvflx(j) +
     &                         cnvflx(j+1))
      end do
!
      cccccc(nrhox)=.25*cnvflx(nrhox-2)+.25*cnvflx(nrhox-1)+
     &               0.5*cnvflx(nrhox)

      do j = 2, nrhox
         cnvflx(j) = cccccc(j)
         flxcnv(j) = cnvflx(j)
      end do
!
! part 1: ------------------------------------------------------!
      do j = 1, nrhox
         rdabh(j) = rdabh(j) - flxrad(j) * dabros(j) / abross(j)
         del      = 1.0d0
         d        = 0.0d0
!------ 
! --------if convflx > 0.0
*
         if (cnvflx(j) .gt. 0.0d0 .and. flxcnv0(j) .gt. 0.0d0) then
            del    = dltdlp(j) - grdadb(j)
            vco    = 0.5d0 * mixlth * sqrt(max(-0.5d0 * ptotal(j) / 
     &               rho(j) * dlrdlt(j), 0.0d0))
            fluxco = 0.5d0 * rho(j) * heatcp(j) * t(j) * mixlth / fourpi
            if (mixlth .gt. 0.0d0 .and. vco .gt. 0.0d0) d = 8.0d0 * 
     &         5.6697d-5 * t(j) ** 4 / (abross(j) * hscale(j) * 
     &         rho(j)) / (fluxco * fourpi) / vco

            taub=abross(j)*rho(j)*mixlth*hscale(j)
            d=d*taub**2/(2.+taub**2)

            d       = d ** 2 * 0.5d0
            ddel(j) = (1.0d0 + d / (d + del)) / del
         end if
!---- end if convflx 

!-- note flxrad equiv to hflux is the eddingtion flux 
         cnvfl = 0.0d0
         if ( cnvflx(j) / flxrad(j) .gt. 0.001d0 .and. 
     &       flxcnv0(j) / flxrad(j) .gt. 0.001d0) cnvfl = cnvflx(j)


         codrhx(j) = (rdabh(j) + cnvfl * (dtdrhx(j) / t(j) * 
     &               (1.0d0 - 9.0d0 * d / (d + del)) + 
     &                1.5d0 * ddlt(j) / del * (1.0d0 + d / 
     &               (d + del)))) / (flxrad(j) + cnvflx(j) * 
     &                1.5d0 * dltdlp(j) / del * (1.0d0 + d / 
     &               (d + del)))
      end do





!--- Part 2: ----------------------------------------------------! 
      codrhx(1) = 0.0d0
      codrhx(2) = 0.0d0
      call integ (rhox, codrhx, g, (nrhox), 0.0d0)
*
      do j = 1, nrhox
         g(j)     = exp(g(j))
         gflux(j) = g(j) * (flxrad(j) + cnvflx(j) - flux) / 
     &              (flxrad(j) + cnvflx(j) * 1.5d0 * dltdlp(j) * 
     &              ddel(j))
*
*BUG!    1.5*DLTDLP(J)/DEL*(1.+D/(D+DEL))) BUG FOUND BY U LONDON OBS
*
      end do
*
      call integ (tauros, gflux, dtau, (nrhox), 0.0d0)
*
      do j = 1, nrhox
         dtau(j)   = dtau(j) / g(j)
!---- since 1rst order migth be not valid for large flux error,
!-- there is a restriction:  
         dtau(j)   = max(-tauros(j)/3.0d0, 
     &                   min(tauros(j)/3.0d0 , dtau(j)))

!---- this was the old restriction, see eq. (7.15) in Bob K. manual
!
!         dtau(j)   = max(-tauros(j) * 0.5d0, 
!     &                   min(tauros(j) * 0.5d0, dtau(j)))
!-------------------------------------------------------------
         dtflux(j) = -dtau(j) * dtdrhx(j) / abross(j)
      end do
!----------------------------------------------------*
!    decrease changes furthermore, in upper regions (where the second approach 
!    should be valid                
!
!      j = 1
!
!      do while (j .le. nrhox .and. tauros(j) .lt. 0.03d0)
!
!         do i = 1, j
!            dtflux(i) = dtflux(i) * 0.5d0
!         end do
!
!         j = j + 1
!      end do
!      dtflux(1) = 0.0d0
!      dtflux(2) = 0.0d0
!-------------------- commented out in the newer version ! **

!-----------------------------------------------------------------------
!
!  PART 3:  ------------------------------------------------



      do j = 1, nrhox
         flxerr(j) = (flxrad(j) + cnvflx(j) - flux) / flux * 100.0d0
      end do
*
      call deriv (tauros, flxerr, flxdrv, (nrhox))
      teff25 = teff / 25.0d0
!
      do j = 1, nrhox
         if (cnvflx(j) / flxrad(j) .lt. 0.001d0) flxdrv(j) = rjmins(j) / 
     &      abross(j) / flux * 100.0d0

         dtlamb(j) = -flxdrv(j) * flux * 0.01d0 / rdiagj(j) * abross(j)
*
         if(cnvflx(j) / flxrad(j) .ge. 1.0d-5 .or. 
     &      tauros(j) .ge. 1.0d0)             then
            dtlamb(j)   = 0.0d0
            dtlamb(j-1) = dtlamb(j-1) * 0.5d0
            dtlamb(j-2) = dtlamb(j-2) * 0.5d0
            dtlamb(j-3) = dtlamb(j-3) * 0.5d0
            dtlamb(j-4) = dtlamb(j-4) * 0.5d0
            dtlamb(j-5) = dtlamb(j-5) * 0.5d0
         end if
*
*.... FUDGE TO AVOID VERY LARGE TEMPERATURE CORRECTIONS
*
         dtlamb(j) = max(-teff25, min(teff25, dtlamb(j)))
      end do

!---PART 4: surface correction ---------------------------*
      dtsur = (flux - flxrad(1)) / flux * 0.25d0 * t(1)
      dtsur = max(-teff25, min(teff25, dtsur))
! 
      do j = 1, nrhox
         dum(j) = dtflux(j) + dtlamb(j)
      end do
!

      call integ (tauros, dum, tinteg, (nrhox), 0.0d0)
      dumm1(1) = 0.1d0
      idum = map1 (tauros, tinteg, (nrhox), dumm1, tone, 1)

      dumm1(1) = 2.0d0
      idum = map1 (tauros, tinteg, (nrhox), dumm1, ttwo, 1)
!--- get T_av 
      tav = (ttwo(1)  - tone(1)) * 0.5d0
      if (dtsur * tav .le. 0.0d0) tav = 0.0d0
      if (abs(tav) .gt. abs(dtsur)) tav = dtsur
      dtsur = dtsur - tav
      dtsur = 0.0d0 
*
      do j = 1, nrhox
         dtsurf(j) = dtsur
         hratio(j) = cnvflx(j) / (cnvflx(j) + flxrad(j))
      end do

!--------- commented out in newer version ----------------
!*
!      idum = map1 (tauros, hratio, (nrhox), 1.0d0, hratio1, 1)
!*
!      do j = 1, nrhox
!
!         dtsurf(j) = dtsur * (1.0d0 - hratio1)
!      end do
!---------------------------------------------------------*


      do j = 1, nrhox
*
*.... REDUCE SURFACE CORRECTION TO THE FRACTION OF RADIATIVE FLUX
*
!         dtsurf(j) = dtsurf(j) * (1.0d0 - hratio(j))
         t1(j) = dtflux(j) + dtlamb(j) + dtsurf(j)
      end do


!------ print stuff ------------------------------------!
      if (ifprnt .gt. 0)
     &   write (16, 100) 'rhox', 't', 'dtlamb', 'dtsurf',
     &      'dtflux', 't1', 'conv/total',  'error', 'deriv',
     &      (j, rhox(j), t(j), dtlamb(j), dtsurf(j), dtflux(j),
     &          t1(j), hratio(j), flxerr(j), flxdrv(j),
     &          j = 1, nrhox)
100     format ( /// 8x, a, 9x, a, 6x, a, 2(3x, a), 6x, a, 3x, a, 6x, //
     &   a, 5x, a, //
     &   (i3, 1pe12.4, 0pf10.1, 4f9.1, 1x, 1pe11.3, 1x, 0p2f10.3))
!---------------------------------------------------------------------!*

! -- PART 5: -- if DeltaT changes sign it is muliplied 0.5 otherwise it is multiplied 1.25
      dum(1) = 0d0
      do j = 1, nrhox

         if(.not. ifconv .or. hratio(j) .le. 0.0d0) then
*
            if(.not. ifconv .or. j .lt. nrhox / 3) then
!            if(.not. ifconv .or. j .lt. nrhox / 2) then
*
               if (iter .gt. 1) then
*
                  if (oldt1(j) * t1(j) .gt. 0.0d0 .and.
     &                abs(oldt1(j)) .gt. abs(t1(j)) ) t1(j) = t1(j) * 
     &                1.25d0
                  if(oldt1(j) * t1(j) .lt. 0.0d0) t1(j) = t1(j) * 
     &                0.5d0
               end if
*
            end if
*
         end if
         if (abs(t1(j))/t(j) .gt. abs(dum(1))) dum(1) = abs(t1(j))/t(j)
*
        oldt1(j) = t1(j)
      end do

      if (dum(1) .lt. errorT) iterat = numits
      print*, 'iteration, max dT/T = ', iterat, dum(1) 
!      if (iterat .eq. numits) print*, 'final max dT/T= ',iterat,dum(1) 
      if (iterat .eq. numits) then 
      end if
*
*.... DETERMINE RHOX CORRECTION TO MAINTAIN CONSTANT TAUROS
*
      do j = 1, nrhox
         tplus(j)  = t(j) + t1(j)
         taustd(j) = 10.0d0 ** (tau1lg + dble(j-1) * steplg)
      end do
   

!--- Prevent extrapolation outward of the surface:-
!*
!      if(tauros(1) .gt. taustd(1)) then
!         tauros(1) = taustd(1)
!         rhoxstd   = min(taustd(1) / abross(1), rhox(1))
!         ratioA     = rhoxstd / rhox(1)
!         p(1)      = ratioA * p(1)
!         prad(1)   = ratioA * prad(1)
!         xne(1)    = ratioA * xne(1)
!         rhox(1)   = rhoxstd
!      end if
!-------------------------------------------------------------
*
      idum = map1 (tauros, t,    (nrhox), taustd, tnew1,  (nrhox))
      idum = map1 (tauros, prad, (nrhox), taustd, prdnew, (nrhox))
      idum = map1 (TAUROS, Pturb,(nrhox), taustd, pturbnew,(nrhox))
      idum = map1 (tauros, tplus, (nrhox), taustd, tnew2, (nrhox))
!
! 
      call ttaup (tnew1, taustd, ab1, ptot1, p1, prdnew, pturbnew, 
     &            vturb, (grav), (nrhox))
!
      call ttaup (tnew2, taustd, ab2, ptot2, p2, prdnew, pturbnew, 
     &            vturb, (grav), (nrhox))
*
      do j = 1, nrhox
         ppp(j) = (ptot2(j) - ptot1(j)) / ptot1(j)
      end do
*
      idum = map1 (taustd, ppp, (nrhox), tauros, rrr, (nrhox))
*
      do j = 1, nrhox
         drhox(j) = rrr(j) * rhox(j)
!         print*, j, drhox(j) 
      end do
!----- commented out in new version !*
!      drhox(1) = 0.0d0
!      drhox(2) = 0.0d0
!      drhox(3) = 0.0d0
!      drhox(4) = drhox(4) / 16.0d0
!      drhox(5) = drhox(5) /  8.0d0
!      drhox(6) = drhox(6) /  4.0d0
!      drhox(7) = drhox(7) /  2.0d0
!--------------------------------*
!

!----- finally correct temperature ------------------------! 
      do j = 1, nrhox
         t(j) = t(j) + t1(j)
      end do

!---------------------------------------------------------
!.... FORCE MONOTONICITY
!
      do j = nrhox - 1, 1, -1
         t(j) = min(t(j), t(j+1) - 1.0d0)
      end do
!*
!*.... MINIMUM TEMPERATURE Set to minimum in ODFs
!*
      do j = 1, nrhox
         t(j) = max(t(j), 1995.0d0)
      end do
!-------------------------------------------
      ifudge = 0
      do j=1, nrhox
        tk(j)=1.38054d-16*t(j)
        hkt(j)=6.6256d-27/tk(j)
        hckt(j)=hkt(j)*2.99792458d10
        tkev(j)=8.6171d-5*t(j)
        tlog(j)=log(t(j))
 
      end do 

      if (ifudge .eq. 1) return
!*
!*.... CHANGE RHOX TO MAINTAIN CONSTANT TAUROS and interpolate on the right scale
!*
      do j = 1, nrhox
         rhox(j) = rhox(j) + drhox(j)
      end do
*
      idum = map1(tauros, rhox,   (nrhox), taustd, dum1, (nrhox))
      idum = map1(tauros, t,      (nrhox), taustd, dum2, (nrhox))
      idum = map1(tauros, p,      (nrhox), taustd, dum3, (nrhox))
      idum = map1(tauros, xne,    (nrhox), taustd, dum4, (nrhox))
      idum = map1(tauros, abross, (nrhox), taustd, dum5, (nrhox))
      idum = map1(tauros, prad,   (nrhox), taustd, dum6, (nrhox))
      idum = map1(tauros, vturb,  (nrhox), taustd, dum7, (nrhox))
      idum = map1(tauros, bmin,   (nrhox), taustd, dum8, (nrhox))
      idum = map1(tauros,pturb,   (nrhox), taustd ,dum9, (nrhox))
        idum=map1(tauros,accrad,  (nrhox), taustd,dum10, (nrhox))
        do  j=1,nrhox
          if(tauros(1).gt.taustd(j))then
            if(dum5(j).lt.0)dum5(j)=abross(1)
          end if
        end do 


      do j = 1, nrhox
         rhox(j)   = dum1(j)
         t(j)      = dum2(j)
         tk(j)     = k * t(j)
         hckt(j)   = hc / tk(j)
         hkt(j)    = h / tk(j)
         tkev(j)   = k_ev * t(j)
         tlog(j)   = log(t(j))
         p(j)      = dum3(j)
         xne(j)    = dum4(j)
         abross(j) = dum5(j)
         prad(j)   = dum6(j)
         pradk(j)  = prad(j) + pradk0
         vturb(j)  = dum7(j)
         bmin(j)   = dum8(j)
         pturb(j)  = dum9(j)
      end do
!
      do i = 1, 6
         idum = map1(tauros, bhyd(1,i), (nrhox), taustd, dum1, 
     &               (nrhox))
!
         do j = 1, nrhox
            bhyd(j, i) = dum1(j)
         end do
!
      end do
!------ smoothing -----------------------------------
      if (smoothing) then 
 
        if ( dum(1) .lt. errorT*1000.0d0 ) then 
           j1smooth = 2
           j2smooth = nrhox-20
        end if 

        if (j1smooth .gt. 0 ) then
           do  j=j1smooth, j2smooth
             t(j)=wtj(1)*t(j-1)+wtj(2)*t(j)+wtj(3)*t(j+1)
           end do
        end if

      end if 

! ---------------------------------------------------!
      do j = 1, nrhox
         tauros(j) = taustd(j)
      end do
*
      end
*
*********** E N D   O F   S U B R O U T I N E   T C O R R **************
*
      subroutine turb
      use types
      use atlcomm
      implicit none
*
!----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.convbl'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.turbpr'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      integer j
*
!------------------------------- EXECUTION -----------------------------

      do j = 1, nrhox
         vturb(j) = (trbfdg * rho(j) ** trbpow + trbsnd * velsnd(j) /
     &      1.0d5 + trbcon) * 1.0d5
         pturb(j) = rho(j) * vturb(j) ** 2 * 0.5d0
      end do
*
      end
*
*********** E N D   O F   S U B R O U T I N E   T U R B ****************
*
      logical function isinf(a)
       implicit none 
       double precision  a

      if ((a*0.d0).ne.0) then
       isinf = .true.
      else
      isinf = .false.
      end if
      return
      end

***********************************************************************

*
      SUBROUTINE DFSSOLVIT(A,N,B,IPIVOT)
      IMPLICIT REAL(KIND=8) (A-H,O-Z)
      PARAMETER (kw=99)
!C     SOLVES LINEAR EQUATIONS
!C     A IS A COMPLETELY FILLED N BY N ARRAY WHICH IS DESTROYED.
!C     B IS THE RIGHT SIDE VECTOR OF LENGTH N AND RETURNS AS THE SOLUTION
!C     IPIVOT IS A SCRATCH AREA OF LENGTH N.
      REAL(kind=8)  A(N*N), B(N), IPIVOT(N)
!      DIMENSION A(1),B(1),IPIVOT(1)
      EQUIVALENCE(AMAX,SWAP,PIVOT,T)
      DO J=1,N
        IPIVOT(J)=0
      END DO

      DO 550 I=1,N
      AMAX=0.

      DO J=1,N
        IF(IPIVOT(J).EQ.1)GO TO 105
        JK=J-N
        DO K=1,N
         JK=JK+N
         IF(IPIVOT(K).EQ.1)GO TO 100
         AA=ABS(A(JK))
         IF(AMAX.GE.AA)GO TO 100
         IROW=J
         ICOLUM=K
         AMAX=AA
  100 Continue   
        END DO 
  105 continue
      end do     

      IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
      IF(IROW.EQ.ICOLUM)GO TO 260
      IRL=IROW-N
      ICL=ICOLUM-N

      DO L=1,N
       IRL=IRL+N
       SWAP=A(IRL)
       ICL=ICL+N
       A(IRL)=A(ICL)
       A(ICL)=SWAP
      end do
      SWAP=B(IROW)
      B(IROW)=B(ICOLUM)
      B(ICOLUM)=SWAP
  260 ICIC=ICOLUM*N+ICOLUM-N
      PIVOT=A(ICIC)
      A(ICIC)=1.
      ICL=ICOLUM-N
      DO 350 L=1,N
      ICL=ICL+N
  350 A(ICL)=A(ICL)/PIVOT
      B(ICOLUM)=B(ICOLUM)/PIVOT
      L1IC=ICOLUM*N-N
      DO 550 L1=1,N
      L1IC=L1IC+1
      IF(L1.EQ.ICOLUM)GO TO 550
      T=A(L1IC)
      A(L1IC)=0.
      L1L=L1-N
      ICL=ICOLUM-N
      DO  L=1,N
       L1L=L1L+N
       ICL=ICL+N
       A(L1L)=A(L1L)-A(ICL)*T
      END DO

      B(L1)=B(L1)-B(ICOLUM)*T
  550 CONTINUE
      RETURN
      END
!---- end of test other solvit function


      subroutine iparcoe(f, x, a, b, c, n)

      use types
      use atlcomm


      implicit none


!.... 2007 MAY - FIXED PROBLEM WHEN n = 3

      real(kind=dp) f(*), x(*)
      real(kind=dp) a(*), b(*), c(*)


!------------------------ o_parcoe VARIABLES -------------------------

      integer:: j, n 

      real(kind=dp) :: d
      real(kind=dp) :: wt

!------------------------ o_parcoe EXECUTION -------------------------

         c(1) = 0.0d0
         b(1) = (f(2) - f(1)) / (x(2) - x(1))
         a(1) = f(1) - x(1) * b(1)

         c(n) = 0.0d0
         b(n) = (f(n) - f(n-1)) / (x(n) - x(n-1))
         a(n) = f(n) - x(n) * b(n)

         if(n .gt. 2) then
            c(2) = 0.0d0
            b(2) = (f(3) - f(2)) / (x(3) - x(2))
            a(2) = f(2) - x(2) * b(2)

            if(n .gt. 3) then
               c(3) = 0.0d0
               b(3) = (f(4) - f(3)) / (x(4) - x(3))
               a(3) = f(3) - x(3) * b(3)
            end if

            do j = 4, n-1 ! THIS IS SKIPPED WHEN N .LE. 4
               d = (f(j) - f(j-1)) / (x(j) - x(j-1))
               c(j) = f(j+1) / ((x(j+1) - x(j)) * (x(j+1) - x(j-1))) -
     &                f(j) / ((x(j) - x(j-1)) * (x(j+1) - x(j))) +
     &                f(j-1) / ((x(j) - x(j-1)) * (x(j+1) - x(j-1)))
               b(j) = d - (x(j) + x(j-1)) * c(j)
               a(j) = f(j-1) - x(j-1) * d + x(j) * x(j-1) * c(j)
            end do

            do j = 4, n-1 ! THIS IS SKIPPED WHEN N .LE. 4

               if(c(j) .ne. 0.0d0) then
                  wt = abs(c(j+1)) / (abs(c(j+1)) + abs(c(j)))
                  a(j) = a(j+1) + wt * (a(j) - a(j+1))
                  b(j) = b(j+1) + wt * (b(j) - b(j+1))
                  c(j) = c(j+1) + wt * (c(j) - c(j+1))
               end if

            end do

            a(n-1) = a(n)
            b(n-1) = b(n)
            c(n-1) = c(n)
         end if

      end subroutine iparcoe

      subroutine integmod  (x, f, fint, n, start, nstart)
      use types
      use atlcomm

      implicit none

!--------------------------- DUMMY ARGUMENTS ---------------------------
!
      double precision  f(*), fint(*), start, x(*)
      integer  n, nstart 
!
!----------------------------- COMMONS ---------------------------------
*
      include 'common.parblk'
*
!--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  db, df, wt
      integer           i, nn, nc 
      double precision a1(maxd), b1(maxd), c1(maxd), d 

*
!-------------------------------  EXTERNALS ----------------------------
*
      external parcoe 
*
!------------------------------- EXECUTION -----------------------------
      nn = n
      nc = nstart
!    make sure nstart is < nn-2
      if(nstart .ge. (nn-2) ) nc=nn-3

      call parcoe (f, x, (nn))

!--- modify the a,b,c coefficients for linear interpolation:


        do i =nc-1, nn-2 
!!LINEAR INTEROPOLATION FORMULA !!LSA
            c(i, 1) = 0.0d0
            c(i, 2) = c(i, 1)

            b(i, 1) = (f(i+1) - f(i)) / (x(i+1) - x(i))
            b(i, 2) = b(i, 1)

            a(i, 1) = f(i)
            a(i, 2) = a(i, 1)
                  d = (f(i+1) - f(i)) / (x(i+1) - x(i)) 
         end do
!
         a(n-1, 1) = f(n)
         b(n-1, 1) = -d  
         c(n-1, 1) = 0.0d0
         a(n-1, 2) = f(n-1)
         b(n-1, 2) = b(n-2, 1)
         c(n-1, 2) = c(n-2, 1)

!---------------------- modification of a, b, c finished --------!





      fint(1) = start
      df      = x(2) - x(1)
      fint(2) = fint(1) + (a(1,1) + (b(1,1) * 0.5d0 + c(1,1) /
     &          3.0d0 * df) * df) * df

!--- do quadratic interpolation until nstart -1

      do i = 2, nc-2 
         wt        = 1.0d0
         if (c(i,1) .ne. 0.0d0 .or.
     &       c(i,2) .ne. 0.0d0 ) wt = abs(c(i,1)) / (abs(c(i,1)) +
     &      abs(c(i,2)))
         df        = (x(i) - x(i+1))
         db        = -df

!QUADRATIC INTERPOLATION
          fint(i+1) = fint(i) + wt * (a(i,2) + (b(i,2) * 0.5d0 + c(i,2)/
     &      3.0d0 * db) * db) * db - (1.0d0 - wt) * (a(i,1) + (b(i,1) *
     &      0.5d0 + c(i,1) / 3.0d0 * df) * df) * df

      end do 

      do i = nc-1, nn-2

!LINEAR INTEGRATION !!LSA
          df        = (x(i+1) - x(i))
       fint(i+1) = fint(i) + (a(i,1) + (b(i,1) * 0.5d0 + c(i,1) /
     &          3.0d0 * df) * df) * df

      end do
!
       df       = x(nn) - x(nn - 1)
       fint(nn) = fint(nn - 1) + (a(nn,1)+(b(nn,1)*0.5d0+c(nn,1) /
     &         3.0d0 * df) * df) * df



      end subroutine integmod

      subroutine cut_ray_out_of_slice(modelfile)

      implicit none

      character(30), intent(in) :: modelfile

      logical :: rn_exists

      integer :: j, k

      integer :: rn, num, ND, nrays, one

      real*8, dimension(:), allocatable :: N, T, p

      real*8 :: x

      one = 1

      inquire(file = './INPUT/rn.inp', exist = rn_exists)

      if (rn_exists) then

          open(unit = 1539, file = './INPUT/rn.inp', action =  'read')

          read(1539, *) rn

          close(1539)

          open(unit = 1540, file = './INPUT/slice.inp', action = 'read')

          read(1540, *) nrays

          do j = 1, nrays

             read(1540, *) num, ND

             if (num == rn) then

                 allocate(N(ND))
                 allocate(T(ND))
                 allocate(p(ND))

                 read(1540, *) (N(k), T(k), p(k), x, x, x, x, k = 1, ND)

                 open(1541, file = trim(adjustl(modelfile)),
     &                      action = 'write')

                 write(1541, '(I1)') one

                 write(1541, '(2(I4,2x))') num, ND

                 do k = 1, ND

                     write(1541, '(7(2x,e15.7))')
     &                     N(k), T(k), p(k), x, x, x, x

                 enddo

                 close(1541)

                 deallocate(N)
                 deallocate(T)
                 deallocate(p)

                 exit

             else

                 read(1540, *) (x, x, x, x, x, x, x, k = 1, ND)

             endif

          enddo

          close(1540)

      endif

      end subroutine
