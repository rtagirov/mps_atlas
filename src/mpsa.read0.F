      logical function readi0 (mode, myrank, comm)
      use types 
      use atlcomm
!
      implicit none
!
!     MODE = 1 COMPUTE A MODEL
!     MODE = 2 READ A PREVIOUSLY CALCULATED MODEL FOR SOME APPLICATION
!     MODE = 20 SAME AS 2 BUT ON ENCOUNTERING END RETURN WITH NRHOX = 0
!
!     FILE = 11 IS OPENED TO READ NON-SOLAR OPACITY TABLES IF THE 
!               "READ KAPPA" INSTRUCTION IS ENCOUNTERED.  IT ASSUMES THAT THE 
!               APPROPRIATE TABLE HAS BEEN MOVED TO A FILE CALLED "KAPPA.ROSS"
!
!     FILE = 12 IS OPENED IN READMOL TO GET MOLECULAR DATA
!
!     FILE = 13 IS OPENED TO READ IN MODEL STRUCTURES UP TO BUT NOT 
!               INCLUDING THE "BEGIN" IF THE "READ PUNCH" OPTION IS 
!               ENCOUNTERED.  ITS GENERIC NAME IS "MPSA.PUNCHED.INPUT"
!
! --- ODFs will be read from .nc files ---(Oct 2018) -----------------!
!     for now only one ODF file with a single value
!    of the microturbulence velocity is read. The name is specified in 
!    the altlas9.control file
!---------------------------------------------------------------------!
!
!  the following has not been changed and can be still used:
!
!               IF MICROTURBULENCE IS TO BE ALLOWED TO VARY WITH DEPTH, 
!               THEN A SEPARATE ODF IS NEEDED FOR EACH VELOCITY.  
!               THIS LEADS TO:
!
!!     FILE 20 = ODF AT 0 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]0
!!               IT GENERIC NAME IS "MPSA.ODF.0"
!!     FILE 21 = ODF AT 1 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]1
!!               IT GENERIC NAME IS "MPSA.ODF.1"
!!     FILE 22 = ODF AT 2 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]2
!!               IT GENERIC NAME IS "MPSA.ODF.2"
!!     FILE 24 = ODF AT 4 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]4
!!               IT GENERIC NAME IS "MPSA.ODF.4"
!!     FILE 28 = ODF AT 8 KM/S MICROTURBULENCE - BDFXXX[BIG/LIT]8
!!               IT GENERIC NAME IS "MPSA.ODF.8"
!!
!!---------------------------- DUMMY VARIABLES --------------------------
!!
      integer mode
!
!------------------------------- COMMONS -------------------------------
!
      include 'common.constb'
      include 'common.abross'
      include 'common.convbl'
      include 'common.depart'
      include 'common.edenbl'
      include 'common.elembl'
      include 'common.fluxbl'
      include 'common.height'
      include 'common.ifblkk'
      include 'common.ifopbl'
      include 'common.iterbl'
      include 'common.junkbl'
      include 'common.musblk'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.steplg'
      
      include 'common.tabptv'

      include 'common.teffbl'
      include 'common.tempbl'
      include 'common.turbpr'
      include 'common.waveyb'
      include 'common.xabund'

#ifdef MPI
      include 'mpif.h'
      integer ierr
#endif 

!
!------------------------------- EXTERNALS -----------------------------
!
      external freeff, integ, map1, readmol, ttaup, wordf, rosstab
      external init_freset 
!
!-------------------------- LOCAL VARIABLES --------------------------
!
      character(3) a, b
      character card*80, ifopc(20)*5, word*8, wordf*80
! --- hardcoded input pahts:
      character(50) inputpunch, kaprosfile, baseODF    
      character(2) varend(5)

!---- ---------------------------------------------------------------

      double precision dummy, exp10, freeff,  rhoxa(maxd), sum, 
     &                 taustd(maxd),  tmin, trat, trat4, 
     &                 tsmooth(maxd), vnew,  wtj, wtjm1, wtjp1, 
     &                 x, rosstab
      double precision deltafr
      
      double precision wavlimlo, wavlimhi     

      integer  i, idum, il, iop(20), ip, it, iterpun, iv, iz,  
     &         j, j1, j2, map1, mu, nn, nnew, nnn, nu, place
      integer  dum1

      integer ncid, ier, myrank, comm, nbin 
      integer ifbigfreq 

      logical  begin, iswch, null, op11, op13, op19, op20,
     &         op21, op22, op24, op28, stopf, ifstandart 
      logical wavelimit

!
!----------------------------- INITIALIZATION --------------------------
!
      data ifopc / 
     &   'H1   ', 'H2+  ', 'H-   ', 'Hray ', 'He1  ', 'He2  ', 'He-',
     &   'Heray', 'Cool ', 'Luke ', 'Hot  ', 'Elec',  'H2ray', 'Hline', 
     &   'Lines', 'Lscat', 'Xline', 'Xlsct', 'Xcont', 'Xscat' /

      data varend / '.0', '.1', '.2', '.4', '.8' /

      exp10(x) = exp(x * 2.30258509299405)

!
!------------------------------- EXECUTION -----------------------------
!

!--- hardcoded file pahts -> should be changed later!

      kaprosfile = './INPUT/kappa.ross'

! ---- initialise: 

! ---- assume that we are lte, unless set otherwise

      nlteon = .false.
      wlte= 'lte' 
      do j = 1, maxd
        do i = 1, 8
           bhyd(j, i) = 1.0
        end do
        bmin(j) = 1.0
      end do

!---- initialise ifop array:
      
      do i = 1,13 
          iop(i) =  1 
          ifop(i) = .true.
      end do
      if (mode .eq. 0 ) ifop(15) = .false.
      if (mode .ne. 0) ifop(15) = .true. 

!     ---------------------

      if (mode .eq. 1) then
!
!       if model atmosphere is calculated the default is: 
         ifcorr = .true.
         ifpres = .true.
         numits = 15 
      else 

!  flux or ODF calculations some of the settings are pre-defined: 
       numits = 1
!      temperature corrections are off: 

       ifcorr = .false.
       ifpres = .false.
!
       do j = 1, maxd
         flxerr(j) = 0.0d0
         flxdrv(j) = 0.0d0
       end do
      end if 

!-------------------------------------
      ifbigfreq = 0
      wavelimit = .false.

      iftabk = .false.
      begin = .false.
      ameanODF = 0

!    ------------------------------------------------------------- 


!---- start reading the input file and change settings      ------!
      do while (.not. begin)
         read (15, 100) card
100      format (a)
         place  = 0
!
         if (index(card(:5), 'abun') .ne. 0) then
!
!.... ABUNDANCE INFORMATION OR CHANGES
!
            if (index(card, 'scal') .ne. 0) then
!
!.... SCALE ALL ABUNDANCES HEAVIER THAN HELIUM
!
               xscale = freeff(card, place)
!
!.... CHECK TO SEE IF THERE IS MORE ABUNDANCE INFORMATION ON THIS CARD
!
               if (index(card, 'chan') .ne. 0) then
                  iz = freeff(card, place)
!
                  do while (place .ne. 80)
                     abund(iz) = freeff(card, place)
                     iz = freeff(card, place)
                  end do
!
               end if
!
            else if (index(card, 'chan') .ne. 0) then
!
!.... CHANGE PARTICULAR ELEMENTS
!
               iz = freeff(card, place)
!
               do while (place .ne. 80)
                  abund(iz) = freeff(card, place)
                  iz = freeff(card, place)
               end do
!
            else if (index(card, 'abso') .ne. 0) then
!
!.... FINAL ABUNDANCES AFTER SCALING
!
               iz = freeff(card, place)
!
               do while (place .ne. 80)
                  abund(iz) = freeff(card, place)
                  if (iz .gt. 2) abund(iz) = abund(iz) / xscale
                  iz = freeff(card, place)
               end do
!
            else
               write (16, 9002) '  in abundance i cannot read ', card
               write (*,  9002) '  in abundance i cannot read ', card
9002           format (a, a)
               stop
            end if
! are the ODFs with filters?

         else if (index(card, 'filter') .ne. 0) then 

            if (index(card, ' on') .ne. 0 ) then

               ifilter = 1  
               print*, 'filter ODFs are going to be calculated'

             else if (index(card, ' off')  .ne. 0 ) then

               ifilter = 0 
               print*, ' no filters are used to calculate ODF'

             end if





! new reading flags for dynamical bins

         else if (index(card, 'user defined') .ne. 0) then 

            if (index(card, ' off') .ne. 0 ) then

               ifkbin = .true. 
               print*, 'the standard Kuruzc bin grid is used'

             else if (index(card, ' on')  .ne. 0 ) then

               ifkbin  = .false.
               print*, 'a  user defined bin grid is used '

             end if

         else if (index(card, 'convect') .ne. 0) then 
!
             if (index(card, 'off') .ne. 0 ) then 
!     
             else if(index(card, 'on') .ne. 0) then  
               ifconv = .true.
               mixlth = freeff(card, place) 
             else 
               print*, ' convection settings unreadable '
               go to 8000
             endif 
!
         else if (index(card, 'turbulence') .ne. 0) then 

             if (index(card, 'off') .ne. 0) then 
              ifturb = .false.
             else if (index(card, 'on') .ne. 0 ) then 
              ifturb = .true. 
              trbfdg = freeff(card, place)
              trbpow = freeff(card, place)
              trbsnd = freeff(card, place)  
              trbcon = freeff(card, place)
 
             endif 
 
         else if (index(card, 'velocit') .ne. 0) then
!
!.... Vturb velosities and how many are needed 
!
            if (index(card, 'valu') .ne. 0) then
               numvt  = freeff(card, place)
!
               do iv = 1, numvt 
                  ivt(iv) = freeff(card, place)
!
                  if (place .eq. 80 .and. iv .lt. numvt ) then
                     read (15, 100) card
                     ivt(iv) = freeff(card, place)
                  end if
!
               end do
!
            else
               write (16, 9001) '  in turb veloc i cannot read ', card
               write ( *, 9001) '  in turb veloc i cannot read ', card
               stop
            end if
9001           format (a,a)
!
!  which of the vturb velocitites should be use for further calculations?

         else if (index(card, 'odfnumb') .ne. 0) then

             neediv = freeff(card,place)

!   do we want to recalculate the electron number density?

         else if (index(card, 'recalxne') .ne. 0 ) then 
             
             if (index(card, ' on') .ne. 0 ) then 

               recalxne = .true. 
               print*, 'recalxne = true '

             else if (index(card, ' off')  .ne. 0 ) then 
           
               recalxne = .false.
               print*, 'recalxne = false '

             end if  


!
         else if (index(card(:5), 'begi') .ne. 0) then
!
!.... BEGIN THE CALCULATION
!
            begin = .true.
!
         else if (index(card(:5), 'end') .ne. 0) then
!
!.... END OF THE CALCULATION
!
            if (mode .eq. 20) then
               nrhox = 0
               readi0 = .true.
!
            else
!
!.... CLOSE THE TAPE FILE IF USING LINE OPACITY
!
               if (ifop(15) .and. (lodf .eq. 'constant' .or. 
     &                             lodf .eq. 'memory  ')) then
                  close (unit=19, status = 'keep')
!
               else if (ifop(15) .and. lodf .eq. 'variable') then
                  close (unit=20, status = 'keep')
                  close (unit=21, status = 'keep')
                  close (unit=22, status = 'keep')
                  close (unit=24, status = 'keep')
                  close (unit=28, status = 'keep')
               end if
!
               readi0 = .false.
            end if
!
            return
!
         else if (index(card(:5), 'freq') .ne. 0) then

!.... GET  FREQUENCY GRIDS
!
            freqid = wordf(card, place)

            if (index(freqid, 'big') .ne. 0) then
!.... FREQUENCY SET BIG.  CURRENT SIZE = 337.  OLD SIZE = 331
!
               numnu = 328
               nulo =1
               nuhi = 328
               ifbigfreq = 1
            else if (index(freqid, 'little') .ne. 0) then
!.... FREQUENCY SET LITTLE.  CURRENT SIZE = 1221.  OLD SIZE = 1215
!
               ifbigfreq = 2
               numnu = 1221 
               nulo = 1 
               nuhi = 1221 

            else
               write (16, 9001) '  in frequency i cannot read ', card
               write ( *, 9001) '  in frequency i cannot read ', card
               stop
            end if

!  set starting wavelength and finishing wavelength for emergent spectra calculations
! -- key word is wavlimit
 
         else if (index(card, 'wave limit') .ne. 0) then  
              wavelimit = .true.
              wavlimlo = freeff(card, place)
              wavlimhi = freeff(card, place) 
!
!
         else if (index(card(:5), 'grav') .ne. 0) then
!
!.... GRAVITY
!
            grav = freeff(card, place)
            if (grav .lt. 10.0d0) grav = exp10(grav)
            glog = log10(grav)
!
         else if (index(card(:5), 'lodf') .ne. 0) then
            lodf = wordf(card, place)
         print *, ' lodf ', lodf
!
         else if (index(card(:5), 'mole') .ne. 0) then
!
!.... MOLECULES
!
            if (index(card, ' on') .ne. 0) then
               ifmol = .true.
!.... READ MOLECULES
!
               call readmol
!
            else if (index(card, ' off') .ne. 0) then
               ifmol = .false.
            else
               write (16, 9001) '  in molecules i cannot read ', card
               write ( *, 9001) '  in molecules i cannot read ', card
            end if
!
         else if (index(card(:5), 'opac') .ne. 0) then
!
!.... OPACITY
!
            if (index(card, 'ifop') .ne. 0) then
!
               do i = 1,20
                  iop(i) = freeff(card, place)
                  if(iop(i) .eq. 0) ifop(i) = .false.
                  if(iop(i) .eq. 1) ifop(i) = .true.
               end do
!
            else
!
               if (index(card, ' on') .ne. 0) then
                  iswch = .true.
!
               else if (index(card, ' off') .ne. 0) then
                  iswch = .false.
!
               else
                  write (16, 9001) '  in opacity i cannot read ', card
                  write ( *, 9001) '  in opacity i cannot read ', card
                  stop
               end if
!
               null = .false.
!
               do while (.not. null)
                  word = wordf(card, place)
!
                  if (word .eq. ' ') then
                     null = .true.
!
                  else
!
                     do i = 1, 20
                        if (word .eq. ifopc(i)) ifop(i) = iswch
                     end do
!
                  end if
!
               end do
!
            end if
!
         else if (index(card(:5), 'pres') .ne. 0) then
!
!.... PRESSURE IS TO BE COMPUTED
!
            if (index(card, ' on') .ne. 0) then
               ifpres = .true.
!
            else if (index(card, ' off') .ne. 0) then
               ifpres = .false.
!
            else
               write (16, 9001) '  in pressure i cannot read ', card
               write ( *, 9001) '  in pressure i cannot read ', card
               stop
            end if
!
         else if (index(card(:5), 'prin') .ne. 0) then
!
!.... PRINTING INSTRUCTIONS
!
            ifprnt = freeff(card, place)
!
         else if ( (index(card(:10), 'punc') .ne. 0)   .and.
     &             (index(card(:10), 'read') .eq. 0)) then
!
!.... "PUNCHING" INSTRUCTIONS
!
            if (allocated(ifpnch)) then  
               ifpnch(1) = freeff(card, place)
             else 
               allocate(ifpnch(1)) 
               ifpnch(1) = freeff(card, place)
            endif 

!            ifpnch = freeff(card, place)
!
         else if (index(card(:5), 'read') .ne. 0) then
!
!.... READ IN VARIOUS KINDS OF INPUT
!
            if (index(card, 'kapp') .ne. 0) then
!
!.... READ NON-SOLAR ROSSELAND OPACITY TABLE
!
               inquire(unit = 11, opened = op11)
!
               if(.not. op11) then
                  open (unit = 11, file = kaprosfile , status = 'old', 
     &                  form = 'formatted', err=7000)
!
               else
                  rewind (unit = 11)
               end if
!
               read (11, 100) card
               write (16, *) card
               read (11, 100)
!
               do it = 1, 56
!
                  read(11,*) 
                  read(11,*)
                  read(11,*)
                  read(11,*) 
                  do ip = 1, 38
                     read (11, 1152) (tabkap(it, ip, iv), iv = 1, 5)
1152                 format (5x, 5x, 5f7.3)
                  end do
!
               end do
!
               iftabk = .true.
               close (unit = 11)

            else 
               write (16, 9001) ' in starting i cannot read ', card
               write ( *, 9001) ' in starting i cannot read ', card
            end if

!-------- scale model if model calculations ------------------------
!----------------------------------------------------------------------
!
         else if (index(card(:5), 'scal') .ne. 0) then
!
!.... SCALE THE T(TAU) FROM A MODEL WITH A DIFFERENT EFFECTIVE
!.... TEMPERATURE AND SURFACE GRAVITY
!
!    three option to scale the atmosphere: 
!    a) tau min, tau max in lg10 scale, and # of points -> steplg will be determined
!    b) tau min, tau max in lg10 scale, and steplg --> number of points will be determined
!    c) as initially implemented: tau min in lg10 scale, steplg, and points -> tau max


            if (index(card, 'x') .ne. 0 ) then 
               
              krhox  = freeff(card, place)
              tau1lg = freeff(card, place)
              taulalg = freeff(card, place)
              teff1  = freeff(card, place)
              gnew   = freeff(card, place)
!   determin steplg
              steplg = (taulalg - tau1lg)/dble(krhox) 
              dum1 = int(steplg*1000.d0)
              steplg = dble(dum1)*0.001

              print*,'case a: tau1, step, taumax, krhox =',
     &               10**tau1lg, steplg, 10**taulalg, krhox
 
            else if (index(card, 'b') .ne. 0 ) then 

              tau1lg = freeff(card, place)
              steplg = freeff(card, place)
              taulalg = freeff(card, place)
              teff1  = freeff(card, place)
              gnew   = freeff(card, place)
!    determin krhox
              krhox = int((taulalg - tau1lg)/steplg) +1 
              print*,'case b: tau1, step, taumax, krhox =',
     &               10**tau1lg, steplg, 10**taulalg, krhox

            else 

              krhox  = freeff(card, place)
              tau1lg = freeff(card, place)
              steplg = freeff(card, place)
              teff1  = freeff(card, place)
              gnew   = freeff(card, place)
              taulalg = tau1lg+dble(krhox)*steplg
              print*,'case c: tau1, step, taumax, krhox =',
     &               10**tau1lg, steplg, 10**taulalg, krhox

! 
            end if 



            if (gnew .lt. 10.0d0) gnew = exp10(gnew)
!
!.... TEST THAT KRHOX <= MAXD
!
            if(krhox .gt. maxd) then
               write(16,*) ' krhox = ', krhox, ', which is > maxd'
               write( *,*) ' krhox = ', krhox, ', which is > maxd'
               stop
            end if
!
!
!
            if (teff .le. 0.0d0) then
               write (16, *) ' must specify teff before scaling'
               write ( *, *) ' must specify teff before scaling'
               stop
!
            end if
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!----- read if arithmetic or geometric mean should be calculated in ODF
!
        else if (index(card(:5), 'mean') .ne. 0) then 



            if (index(card, 'arith')  .ne. 0) then 
              print*, 'ODF will be calculated with arithmetic mean'
              ameanODF = 1
            else 
              print*, 'ODF will be calculated with geometric mean'

              ameanODF = 0

            endif 


         else if(index(card(:5), 'overs') .ne. 0 ) then 
           
              if (index(card, ' on') .ne. 0 ) then

               overwt = 1.0d0 
               print*, 'overshoot is on'

             else if (index(card, ' off')  .ne. 0 ) then

               overwt  = 0.0d0
               print*, ' overshoot is off'

             end if 
             



         else if (index(card(:5), 'surf') .ne. 0) then
!
!.... SURFACE FLUXES OR INTENSITIES ARE DESIRED
!
            if (index(card, 'inte') .ne. 0) then
               nmu = freeff(card, place)
!
               do mu = 1, nmu
                  angle(mu) = freeff(card, place)
!
                  if (place .eq. 80 .and. mu .lt. nmu) then
                     read (15, 100) card
                     angle(mu) = freeff(card, place)
                  end if
!
               end do
!
               ifsurf = 2
!
            else if (index(card, 'flux') .ne. 0) then
               ifsurf = 1
!
            else if (index(card, ' off') .ne. 0) then
               ifsurf = 0
!
            else
               write (16, 9001) '  in surface i cannot read ', card
               write ( *, 9001) '  in surface i cannot read ', card
               stop
            end if
!
         else if (index(card, 'starting m') .ne. 0) then 
!---- read in teff and surface gravity of starting model 

          if (index(card, 'teff') .ne. 0) then
!.... TEFF
            teff = freeff(card, place)
            flux = sigma / fourpi * teff ** 4
!
!.... CHECK TO SEE IF GRAVITY SPECIFICATION IS ON THIS LINE
!
            if (index(card, 'grav') .ne. 0) then
               grav = freeff(card, place)
               if (grav .lt. 10.0d0) grav = exp10(grav)
               glog = log10(grav)
            end if
!
           end if 

         else if (index(card(:5), 'iter') .ne. 0) then 

           numits = freeff(card, place) 
           if (allocated(ifpnch)) then 
            ! check if the size is right... if not
             if (size(ifpnch,1) .lt. numits) then 

                deallocate(ifpnch) 
                allocate(ifpnch(numits))
             endif 

           else 
                allocate(ifpnch(numits))  
           endif 
 
                ifpnch = 0
                ifpnch(numits) = 1
                print*,'We will do iterations = ', numits 

         else if (index(card(:5), 'titl') .ne. 0) then
!
!.... TITLE OF THIS RUN
!
            i = index(card, ' ')
            title = card(i:)
!
         end if
!
      end do
!------------------------------------------------------------
!  end of reading in the input file ( main  while loop) 
!--------------------------------------------------------------!





!--------------------------------------------------------------!
!.... IF LINE OPACITY IS TO BE USED, OPEN THE FILE and read the ODF 
!
      if (mode .eq.0) ODFfile = trim(fodfodf)
!   ODFs are to be read 
!   use 'ODFfile' to get: frequency grid, p, T grid, and ODFs
      if (ifop(15)) then

!--- the old ODFs in binary format are used, then ifkbin 
!     has to be set true the file will be open here, the gird will be set further below
!     and the opacities will be read in mpsa.kapp

!---- set the right name of the ODF to read : 

        if (mode .eq.1) ODFfile = trim(fmodelodf)
        if (mode .eq.2) ODFfile = trim(ffluxodf)
           
        if (lodf .eq. 'variable') then 
          baseODF = trim(ODFfile)
          ODFfile = trim(ODFfile)//'.0'
        endif 
 
!   First open the ODF tables: 
 
        if (ifkbin) then 
  
         if (lodf .eq. 'variable') then
           
           inquire(unit = 20, opened = op20)
           if (.not. op20) open(unit = 20, file = trim(ODFfile),
     &      status = 'old', form = 'unformatted', err=7500)
        
           ODFfile = trim(baseODF)//'.1' 

           inquire(unit = 21, opened = op21)
           if (.not. op21) open(unit = 21, file =  trim(ODFfile),
     &      status = 'old', form = 'unformatted', err=7500)
         
           ODFfile = trim(baseODF)//'.2' 

           inquire(unit = 22, opened = op22)
           if (.not. op22) open(unit = 22, file = trim(ODFfile),
     &      status = 'old', form = 'unformatted', err=7500)

           ODFfile = trim(baseODF)//'.4' 
           inquire(unit = 24, opened = op24)
           if (.not. op24) open(unit = 24, file = trim(ODFfile),
     &      status = 'old', form = 'unformatted', err=7500)
 
           ODFfile = trim(baseODF)//'.8' 
           inquire(unit = 28, opened = op28)
           if (.not. op28) open(unit = 28, file = trim(ODFfile),
     &      status = 'old', form = 'unformatted', err=7500)


         else 

          inquire(unit = 19, opened = op19)

          if (.not. op19) open(unit = 19, file = ODFfile, 
     &      status = 'old', form = 'unformatted', err=7200)

         endif 


!     else if the new nc files are to be read, it will be open here, the
!     frequencies will be allocated and the opacities will be read here
 

        else 
          ifstandart = .false.
!       ibin is the size of the array with the frequency boundaries, thus it is ibin+1
          if (lodf .eq. 'variable' ) then 

           call OpenODF(ncid, trim(ODFfile),  ibin, isubbin,
     &               numpres, numt, vturb,ier)
!       ----- allocate the right frequency grid 
            call init_freset(isubbin,ibin,numt, numpres, ifstandart,1)
!  tabt, tabp, and sbwt need to be read as well!
            call ReadODF(ncid,myrank,comm, kapwvar,inifreset,tabt,tabp,
     &               sbwt,  ibin, numpres, numt, isubbin, ier)
            call CloseNetCDF( myrank, ncid ,ier )
    
             
            do i = 1, 4
              ODFfile =trim(baseODF)//trim(varend(i+1))
              call OpenODF(ncid, trim(ODFfile),  ibin, isubbin,
     &               numpres, numt, vturb,ier)
              call ReadODF(ncid,myrank,comm, kapwvar(1,1,1,1,1+i), 
     &         inifreset,tabt,tabp, sbwt, ibin,numpres,numt,isubbin,ier)
              call CloseNetCDF( myrank, ncid ,ier )
            end do 


          else 
            call OpenODF(ncid, trim(ODFfile),  ibin, isubbin,  
     &               numpres, numt, vturb,ier)
!       ----- allocate the right frequency grid 
            call init_freset(isubbin,ibin,numt, numpres, ifstandart,0)
!  tabt, tabp, and sbwt need to be read as well!
            call ReadODF(ncid,myrank,comm, kapw, inifreset, tabt, tabp,
     &               sbwt,  ibin, numpres, numt, isubbin, ier)
            call CloseNetCDF( myrank, ncid ,ier )
          endif 




!!  frequency grid from bin boundaries needs to be calculated!
          do nu = 1, ibin
             deltafr = inifreset(nu+1) - inifreset(nu)
             freset(nu) = c_nm / (inifreset(nu) +(0.5d0*deltafr))
          end do
!
          do nu = 2, ibin -1
             rcoset(nu) = (freset(nu-1) - freset(nu+1)) * 0.5d0
          end do
!
            rcoset(1) = (c_nm / inifreset(1) - freset(2)) * 0.5
            rcoset(ibin) = freset(ibin-1) * 0.5
        end if 

      end if   

! -----------------------------------------------------------------------------------
!    For the standard binning, set the frequency grid: 
! ----------------------------------------------------------------------------------  
      if ((.not. ifop(15)) .or. ifkbin) then
             ifstandart = .true. 
             numt = 57
             numpres = 25
             isubbin = 12 
             call init_freset(isubbin,1563,numt,numpres,ifstandart, 0)

       if (ifbigfreq .eq. 1 ) then
               ibin = 328
               do nu = 1, numnu
                  freset(nu) = c_nm / inifreset(nu+5)
               end do
!
               do nu = 2, numnu - 1
                  rcoset(nu) = (freset(nu-1) - freset(nu+1)) * 0.5
               end do
!
               rcoset(1) = (c_nm / 8.97666 - freset(2)) * 0.5
               rcoset(numnu) = freset(numnu-1) * 0.5
       else if (ifbigfreq .eq. 2) then 
               ibin = 1221
               do nu = 1, numnu
                  freset(nu) = c_nm / inifreset(nu+342)
               end do
!
               do nu = 2, numnu - 1
                  rcoset(nu) = (freset(nu-1) - freset(nu+1)) * 0.5
               end do
!
               rcoset(1) = (c_nm / 8.97666 - freset(2)) * 0.5
               rcoset(numnu) = freset(numnu-1) * 0.5

       end if 
      end if 
! ----------------------------------------------------------------------------------
!-- check if we calculate the whole range of bins or if there are wavelength limits
!-----------------------------------------------------------------------------------
      if (mode .ne. 0 ) then 
      if (wavelimit) then
          nu = 1
          do while ((freset(nu) .gt. c_nm/wavlimlo) .and. nu .lt. ibin)
            nu = nu +1
          end do

          if ( nu .eq.  ibin) then
            print*,'something wrong,setting nulo and nuhi to default'
            numnu = ibin
            nulo = 1
            nuhi = ibin
          else
            nulo = nu
            nu = 1
            do while((freset(nu) .ge. c_nm/wavlimhi) .and. nu .lt. ibin)
              nu = nu +1
            end do

            if (nu .eq. nulo) nu = nu +1
               nuhi = nu
               numnu = (nuhi - nulo) +1
          end if
          print*, 'wav1,wav2=', c_nm/freset(nulo), c_nm/freset(nuhi)
!    if no bounds are set take the whole range!
      else
          numnu = ibin
          nulo = 1
          nuhi = ibin
          print*, 'wav1,wav2 =', c_nm/freset(nulo), c_nm/freset(nuhi)
      end if


      end if ! --- end of the if for not ODF

!----   NOW everything that needed to be read in is read in ----- ! 

!-------------------------------------------------------------
      if (mode .eq. 1) then
!
!.... PROCESS THE BEGINNING OF THE MODEL
!
         stopf = .false.
!
         if (numnu .eq. 0) then
            write (16, *) '  how many frequencies?'
            write ( *, *) '  how many frequencies?'
            stopf = .true.
         end if
!
         if (teff .eq. 0.0d0) then
            write (16, *) '  what teff?'
            write ( *, *) '  what teff?'
            stopf = .true.
         end if
!
         if (grav .eq. 0.0d0) then
            write (16, *) '  what gravity?'
            write ( *, *) '  what gravity?'
            stopf = .true.
         end if
!
         if (stopf) stop
!
      end if 
!------------------------------------------------------!


!------------------------------------------------------!
      if (abund(1) .lt. 0.0d0) abund(1) = exp10(abund(1))
      if (abund(2) .lt. 0.0d0) abund(2) = exp10(abund(2))
      xabund(1) = abund(1)
      xabund(2) = abund(2)
!
      do iz = 3, 99
         if (abund(iz) .gt. 0.0d0) then
            xabund(iz) = abund(iz) * xscale
            abund(iz)  = log10(abund(iz))
!
         else
            xabund(iz) = exp10(abund(iz)) * xscale
         end if
!
      end do
!
      wtmole = 0.0d0
!
      do iz = 1, 99
         wtmole = wtmole + xabund(iz) * atmass(iz)
      end do
!
      write (16, 1512) 'teff =', teff, 'log g =', glog, wlte,
     &   'title: ', title,
     &   'xscale =', xscale, (elem(iz), abund(iz), iz = 1, 2),
     &    (elem(iz), abund(iz), iz = 3, 99)
1512  format ( a, f7.0, 3x, a, f10.4, 3x, a /
     &         a, a /
     &         a, f10.6, 2(3x, a2, f8.5) /
     &        (7(a2, f7.2, 2x)))
!
!
      write (16, 1517) (ifopc(i), ' =', ifop(i), i = 1, 20)
1517  format (/, 8(a, a, l2, 1x))





!------ Write out all other settings -------------------------------!
      write (16, 1518) 'ifcorr =',  ifcorr,  'ifpres =', ifpres,
     &                 'ifmol  =',  ifmol,   'ifsurf =', ifsurf,   
     &                 'ifconv =', ifconv,   'mixlth =', mixlth, 
     &                 'overshoot =', overwt,
     &                 'ifturb =', ifturb,   'trbfdg =', trbfdg,  
     &                 'trbpow =', trbpow,   'trbsnd =', trbsnd,   
     &                 'trbcon =', trbcon,
     &                 'ifmasks =', ifmasks, 'lodf = ', lodf
1518  format (/ 3(a, l2, 2x), a, i2 /
     &            a, l2, 2x, a, f6.2, 2x, a, f6.2 /
     &            a, l2, 2x,  4(a, f6.2, 2x) /
     &            a, l2, 1x, a, a)


!-------  only for model calculations:
!----------------------------------------------------------!
      if (mode .eq. 1) then
         write (16, *) '  ifprnt', ifprnt,
     &                   '  ifpnch', ifpnch
!
         if (ifwave) then
            write (16, 1536) 'wbegin in ATLAS',wbegin,' deltaw', deltaw,
     &                      '   numnu', numnu
1536        format (a, f11.4, a, f7.4, a, i5)
!
         else
            write (16, 1561) 'freqid = ', freqid, '   numnu =', numnu,
     &                      '   nulo =', nulo,   '   nuhi =', nuhi
1561        format ( / a, a, 3(a, i4) /)
            nn  = (numnu + 3) / 4
            nnn = numnu - nn * 3
            if (nnn .lt. 1) nn=1
            if (nnn .lt. 1) nnn=1
            write (16, 1563) ((nu, freset(nu), rcoset(nu),
     &                         nu = i, numnu, nn), i = 1, nnn)
1563        format ((4(i4, 1p2e13.6)))
!
            if (nn .ne. nnn) then
               nnn = nnn + 1
               write (16, 1564) ((nu, freset(nu), rcoset(nu),
     &                            nu = i, numnu, nn), i = nnn, nn)
1564           format ((3(i4,1p2d13.6)))
            end if
!
            sum = 0.0d0
!
            do nu = nulo, nuhi
               sum = sum + rcoset(nu)
            end do
!
            write(16, 9501) 'sum of the integration weights =', sum
9501        format(a, 1pe13.6)
         end if
!
      end if
!------------------end for model calculations ------------!
!
      readi0 = .true.
      if(myrank .eq. 0) then 

      print*, ' readi0 finished with ', readi0
      print*, ' -------------------------------'
      print*, ' -------------------------------'
      print*, ' -------------------------------'
      print*, '                                '
      end if 


!--------------------------------------------------------!


      return 


!---- deal with not found files: 

7000  continue
      print*, 'can not find file ', kaprosfile 
      go to 8000 

7100  continue
      print*, 'can not find file ', inputpunch
      go to 8000

7200  continue
      print*, 'can not find file ', ODFfile 
      go to 8000 

7500  continue
      print*, 'one of the variable odfs is missing' 
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
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   R E A D I N !!!!!!!!!!!!
!
      subroutine initialise_param(comm, flagODF, flagmodel, flagflux)
       use atlcomm 

       implicit none
#ifdef MPI
      include 'mpif.h'
      integer ierr
#endif 

       integer, intent(in) :: comm
       integer, intent(out) ::  flagODF, flagmodel, flagflux

       integer place, odfcount, fluxcount, modelcount 
       character(len=80) card 
       logical fileend 
        
!------------------------------- EXTERNALS -----------------------------
!
       external freeff,  wordf
       real(kind=8) freeff
       character(len=80) wordf
!  initialise:

       flagODF=0 
       flagflux= 0 
       flagmodel = 0 

       odfcount = 0 
       fluxcount = 0 
       modelcount = 0 


       accrad = 0.0
       prad   = 0.0

       fileend = .false. 
! ---- open control file 

       open (unit = 2, file = 'mpsa.control', form = 'formatted',
     &      status = 'old', err=8500)

!
       do while (.not. fileend )
         read (2, 100) card
100      format (a)
         place  = 0
! ---- first check for input file names: 

         if(index(card, 'odf_input') .ne. 0) then
            fileODFinput = wordf(card,place)
            odfcount = odfcount +1

         else if(index(card, 'flux_input') .ne. 0) then 
            filefluxinput = wordf(card,place)
            fluxcount = fluxcount +1          
         else if(index(card, 'model_input') .ne. 0) then 
            filemodelinput = wordf(card,place)
! ------- second check for model and grid names:
            modelcount = modelcount +1 
         else if(index(card, 'odf_tp_grid') .ne. 0) then 
            filePTgrid = wordf(card,place)
            odfcount = odfcount +1
         else if(index(card, 'model_start') .ne. 0) then 
            filestartmodel = wordf(card,place)
            modelcount = modelcount +1
         else if(index(card, 'flux_model') .ne. 0) then 
            filefluxmodel = wordf(card,place)
            fluxcount = fluxcount +1
!-------- third check for odf files:
! --- NOTE this needs to be changed a bit later:
 
         else if (index(card, 'odf_output')  .ne. 0) then 
            fodfodf  = wordf(card,place)
            odfcount = odfcount +1
         else if (index(card, 'model_odf') .ne. 0) then 
            fmodelodf  = wordf(card, place)
            modelcount=modelcount +1 
         else if (index(card, 'flux_odf') .ne. 0) then 
            ffluxodf = wordf(card,place) 
            fluxcount = fluxcount +1
         else  if(index(card, 'calc_odf') .ne. 0) then 
             flagodf = 1 
             if (index(card, 'off') .ne. 0) flagodf = 0  

         else if(index(card, 'calc_flux') .ne. 0) then 
             flagflux = 1 
             if (index(card, 'off') .ne. 0) flagflux = 0 

         else if (index(card, 'calc_model') .ne. 0) then 
             flagmodel =1 
             if (index(card, 'off') .ne. 0) flagmodel = 0 

!----- check if control file is finished :

         else if (index(card, 'endfile') .ne. 0) then
            fileend  = .true.

         endif



       end do
!------- paranoid check: 
       close(unit =2)


!      first check that one of the modules is set to be calculated:
       if (max(max(flagodf, flagflux), flagmodel) .lt. 1 ) then 
        print*, ' No module is set to be calculated. ' 
        go to 9000
!----  check if input file for the modules that needs to be calculated are set:

       else  
     
         if(flagodf .eq. 1) then 
           if (odfcount .lt. 3) then 
            print*, 'not enough files are specified to calculate ODFs'
            go to 9000
           endif 
 
         endif 

         if (flagflux .eq. 1) then 
           if (fluxcount .lt. 3) then
            print*, 'not enough files are specified to calculate  flux '
            go to 9000
           endif
  
         endif 

         if (flagmodel .eq. 1) then 
           if (modelcount .lt. 3) then
            print*,'not enough files are specified to calculate a model'
            go to 9000
           endif

         endif 
 
 
       end if 
 
      
! - ------------------------------------------------------------

       return 

8500   continue
      print*, 'can not control file '

9000   continue
#ifdef MPI 
      print*,'going to abort ..'
      call mpi_abort(comm, ierr)

#else
      print*, 'going to abort ...'
      stop
#endif 

      end subroutine initialise_param


      subroutine get_data_location(mpif, datafolder) 
       implicit none 

!------------------------------- EXTERNALS -----------------------------
!
       external freeff,  wordf
       real(kind=8) :: freeff
       character(len=80) :: wordf
       
! ----   internal 
       character(256), intent(out) :: datafolder
       character(len=80) :: card

       logical :: fileend, readinloc 
       integer, intent(in):: mpif
       integer :: place 
! ------------------------------------------------------- 

       readinloc = .false. 
       fileend = .false.
! ---- open control file 
       open (unit = 2,file='mpsa.control',form='formatted',status='old')


       do while (.not. fileend )
         read (2, 100) card
100      format (a)
         place  = 0
! ---- first check for input file names: 

         if(index(card, 'data_folder') .ne. 0) then
            datafolder = wordf(card,place)
            readinloc = .true.
            fileend = .true. 

         else if (index(card, 'endfile') .ne. 0) then
            fileend  = .true.
  
         endif 

       end do 

       if (.not. readinloc) then 
         if (mpif .eq. 0) then 
          datafolder = '../' 
         elseif (mpif .eq. 1) then 
          datafolder = '../../'
         else 
           print*, ' problems locating the data/linelists ' 
           print*, ' going to abort ' 
           stop
          endif 
 
       endif 



       close(unit =2 )
       return   


      end subroutine get_data_location




! --------------------------------------------------------------------
      double precision function freeff (card, place)
      implicit none
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      character(len=80)   card
      integer  place
!
!-------------------------- LOCAL VARIABLES --------------------------
!
      character  copy*80, nbr*13
      logical num
      integer p, w, wc, i
!
!----------------------------- INITIALIZATION --------------------------
!
      save  copy, nbr, p
      data  copy /' '/,  nbr /'0123456789+-.'/
!
!------------------------------- EXECUTION -----------------------------
!
      if (copy .ne. card) then
         copy = card
         p = 0
         place = p
      end if
!
!.... LOCATE THE BEGINNING OF THE NEXT NUMBER
!
      num = .false.
      p = place
!
      do while (p .lt. 80 .and. .not. num)
         p = p + 1
!
         do  i = 1,13
            if (nbr(i:i) .eq. copy(p:p)) num = .true.
         end do
!
      end do
!
      if (p .eq. 80) then
         place = p
         freeff = 0
!
      else
!
!....   LOCATE THE END OF THE NUMBER.  TERMINATED BY COMMA OR BLANK
!
         w  = index (copy(p:), ' ') - 1
         wc = index (copy(p:), ',') - 1
!
         if (w .gt. 0 .and. wc .gt. 0) then
            w = min (w, wc)
!
         else
            w = max (w, wc)
         end if
!
         read (copy(p:p+w-1), *) freeff
!
!....    SET UP FOR THE NEXT CALL WITH THIS CARD
!
         p = p + w
         place = p
      end if
!
      end
!
!!!!!!!!!!! E N D   O F   F U N C T I O N   F R E E F F !!!!!!!!!!!!!!!!
!
      character*80 function wordf (card, place)
c     implicit undefined (a-z)
      implicit none
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      character*80 card
      integer  place
!
!-------------------------- LOCAL VARIABLES --------------------------
!
      character copy*80, letter*52
      integer i, p, w, wc
      logical start
!
!---------------------------- INITIALIZATION --------------------------
!
      save      copy, letter, p
      data      copy /' '/,
     &          letter(:26) / 'abcdefghijklmnopqrstuvwxyz' /
     &          letter(27:) / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
!
!------------------------------- EXECUTION -----------------------------
!
      if (copy .ne. card) then
         copy = card
         p = 1
         place = p
      endif
!
!.... FIND THE FIRST BLANK OR COMMA
!
      start = .false.
      p = place
!
      do while (p .lt. 80 .and. .not. start)
         p = p + 1
         if (copy(p:p) .eq. ' ' .or. copy(p:p) .eq. ',') start = .true.
         if (copy(p:p) .eq. '=' ) then 
          start = .true.
          p = p -1
         endif 

      end do
!
!.... FIND THE NEXT LETTER
!
      start = .false.
!
      do while (p .lt. 80 .and. .not. start)
         p = p + 1
!
         do i = 1, 52
            if (letter(i:i) .eq. copy(p:p)) start = .true.
         end do
!----- this is needed for the mpsa.control file to find filenames

         if (copy(p:p) .eq. '=') then 
           start =.true. 
           p = p+1
           do while (copy(p:p) .eq. ' ') 
            p = p+1
           end do 
         endif 
! --------------------------------------------------------------!
      end do
!
!.... IF AT THE END OF THE CARD RETURN WITH A NULL WORD
!
      if (p .eq. 80) then
         wordf = ' '
!
      else
!
!....   LOCATE THE END OF THE WORD.  TERMINATED BY COMMA OR BLANK
!
         w  = index (copy(p:), ' ') - 1
         wc = index (copy(p:), ',') - 1
!
         if (w .gt. 0 .and. wc .gt. 0) then
            w = min (w, wc)
!
         else
            w = max (w, wc)
         end if
!
         read (copy(p:p+w-1), 100) wordf
100      format (a)
!
!....    SET UP FOR THE NEXT CALL WITH THIS CARD
!
         p = p + w - 1
         place = p
      end if
!
      end
!
!!!!!!!!!!! E N D   O F   F U N C T I O N   W O R D F !!!!!!!!!!!!!!!!!!
!
!---- I think this has to be removed:

#ifdef MODF

      double precision function xfreef(xfreef1)
      implicit none
      include 'common.freebl'
     
      external freeff     
!---------------------------- DUMMY VARIABLES --------------------------
      character(len=80)  card 
      dimension xfreef(30)
      double precision xfreef1(30)
 




      end


#endif



!
      block data  blocke
c     implicit undefined (a-z)
      implicit none
!
!------------------------------- COMMONS -------------------------------
!
      include 'common.elembl'
!
      data atmass/ 1.008, 4.003, 6.939, 9.013, 10.81, 12.01, 14.01,
     &             16.00, 19.00, 20.18, 22.99, 24.31, 26.98, 28.09,
     &             30.98, 32.07, 35.45, 39.95, 39.10, 40.08, 44.96,
     &             47.90, 50.94, 52.00, 54.94, 55.85, 58.94, 58.71,
     &             63.55, 65.37, 69.72, 72.60, 74.92, 78.96, 79.91,
     &             83.80, 85.48, 87.63, 88.91, 91.22, 92.91, 95.95,
     &             99.00, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8,
     &             118.7, 121.8, 127.6, 126.9, 131.3, 132.9, 137.4,
     &             138.9, 140.1, 140.9, 144.3, 147.0, 150.4, 152.0,
     &             157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0,
     &             175.0, 178.5, 181.0, 183.9, 186.3, 190.2, 192.2,
     &             195.1, 197.0, 200.6, 204.4, 207.2, 209.0, 210.0,
     &             211.0, 222.0, 223.0, 226.1, 227.1, 232.0, 231.0,
     &             238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0,
     &             254.0/
!
      data elem/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ',
     &           'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',
     &           'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co',
     &           'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
     &           'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
     &           'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe',
     &           'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
     &           'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
     &           'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
     &           'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     &           'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es'/
      end
!
!!!!!!!!!!! E N D   O F   B L O C K   D A T A   B L O C K E !!!!!!!!!!!!
!

     


      block data  blockr
      use types
      use atlcomm

!     implicit undefined (a-z)
      implicit none
!
!------------------------------- COMMONS -------------------------------
!
      include 'common.convbl'
      include 'common.depart'
      include 'common.fluxbl'
      include 'common.ifblkk'
      include 'common.ifopbl'
      include 'common.ionsbl'
      include 'common.iterbl'
      include 'common.junkbl'
      include 'common.musblk'
      include 'common.rhoxbl'
      include 'common.steplg'
      include 'common.teffbl'
      include 'common.turbpr'
      include 'common.waveyb'
!
!--------------------------- INITIALIZATION --------------------------
!
!
      data  bal1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, 
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, 
     &             maxd*1.0d0/
!
      data  bal2 / maxd*1.0d0 /
!
      data  bc1 /  maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0/
!
      data  bc2 /  maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0/
!
      data  bca1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bca2 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bfe1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bhe1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0/
*
      data  bhe2 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, 
     &             maxd*1.0d0, maxd*1.0d0/
!
      data  bhyd / maxd*1.0, maxd*1.0, maxd*1.0, maxd*1.0, 
     &             maxd*1.0, maxd*1.0, maxd*1.0, maxd*1.0/
!
      data  bmin / maxd*1.0 /
!
      data  bk1 /  maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bmg1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bmg2 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0/                            
!
      data  bna1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bo1 /  maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0/
!
      data  bo2 /  maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
*
      data  bsi1 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0/
!
      data  bsi2 / maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0, maxd*1.0d0, maxd*1.0d0,
     &             maxd*1.0d0, maxd*1.0d0/
!
      data  dlrdlt,   dltdlp,   flxcnv,   grdadb,
     &      heatcp,   hscale,   velsnd,   vconv    / 
     &      maxd*0.0, maxd*0.0, maxd*0.0, maxd*0.0, 
     &      maxd*0.0, maxd*0.0, maxd*0.0, maxd*0.0 /
!
      data  flux / 0.0 /
*
      data  flxerr, flxdrv / maxd*0.0, maxd*0.0 /
!
      data  grav / 0.0 /
!
      data  ifconv, mixlth, overwt / .false., 1.0d0, 0.0d0 /
!
      data  ifcorr,  ifmasks,  ifmol,   ifpres, ifsurf /
     &      .true.,  .false., .false.,  .true.,    0   /
!
      data  ifop /
!
!             H1         H2+        H-         Hray       He1     
     &       .true.,    .true.,    .true.,    .true.,    .true.,  
!
!             He2        He-        Heray      Cool       Luke    
     &       .true.,    .true.,    .true.,    .true.,    .true.,
!
!             Hot        Elec       H2ray      Hline      Lines   
     &       .true.,    .true.,    .true.,    .false.,   .true.,
!
!             Lscat      Xline      Xlsct      Xcont      Xscat   
     &       .false.,   .false.,   .false.,   .false.,   .false. /
!
      data  ifprnt / 2 /
!
      data  ifturb, trbcon, trbfdg, trbpow, trbsnd /
     &      .false., 0.0d0,  0.0d0,  0.0d0,  0.0d0  /
!
      data  ifwave / .false. /
!
      data  lodf / 'constant' /
!
      data  nlteon / .false. /
!
      data  krhox / 0 /
!
      data  nrhox / 0 /
!
!      data  numnu / 0 /
!
      data  pturb / maxd*0.0d0 /
!
      data  surfi / maxmu*0.0d0 /,  nmu / 1 /
!
      data  steplg / 0.125d0 /
!
      data  tau1lg / -5.875 /
!
      data  teff / 0.0 /
!
      data  title / ' ' /
!
      data  wlte / 'lte ' /
!
      data  vturb / maxd*0.0d0 /
!
      data  xnfph / maxd*0.0d0, maxd*0.0d0/
!
      data  xnfphe / maxd*0.0d0, maxd*0.0d0, maxd*0.0d0 /
!
      data  xscale / 1.0d0 /
      end
*
!!!!!!!!!!! E N D   O F   B L O C K   D A T A   B L O C K R  !!!!!!!!!!!
!
      subroutine readmol
      use types
      use atlcomm

!     implicit undefined (a-z)
      implicit none
!
!.... MAKE TABLE OF ALL COMPONENTS OF ALL MOLECULES
!
!.... FILE = 12 IS THE FILE OF MOLECULAR EQUILIBRIUM CONSTANTS
!
!.... SAMPLE CODES FOR ATOMS AND MOLECULES
!                             EXTERNAL CODE    INTERNAL COMPONENTS
!          CARBON DIOXIDE...    60808.0             6,8,8
!          HMINUS...........    100.0               1,100
!          NEUTRAL IRON.....    26.0                26
!          H2PLUS...........    101.01              1,1,101
!          HYDROGEN ION.....    1.01                1,101
!          SILICON 3+.......    14.03               14,101,101,101
!
!------------------------------- COMMONS -------------------------------
!
!      include 'common.sizebl'
      include 'common.ifequa'
!
!------------------------- LOCAL VARIABLES ----------------------------
!
      double precision  c, e1, e2, e3, e4, e5, e6, e7, x, xcode(8)
      integer           i, id, iequa, ii, ion, jmol, kloc
      logical           mormol, op12
! --- hardcoded file names
      character(30) inputmolec 

!
!-------------------------- INITIALIZATION ----------------------------
!
!      data  ifequa / 101 * 0 /
!
      data  xcode / 1.0d14, 1.0d12, 1.0d10, 1.0d8, 
     &              1.0d6,  1.0d4,  1.0d2,  1.0d0  /

!
*------------------------------- EXECUTION -----------------------------
!
      inputmolec = './INPUT/molecules.dat'


      do i = 1,101
       ifequa(i) = 0
      end do 


      do jmol = 1, maxmol
         code(jmol) = 0.0d0
      end do
!
      inquire(unit = 12, opened = op12)
      if (.not. op12) open (unit = 12, file = trim(inputmolec), 
     &   status = 'old', form = 'formatted', err=7600)
      write (16, 10) '1molecules input'
   10 format (a)
!
!.... IF IFEQUA=1 AN EQUATION MUST BE SET UP FOR ELEMENT I
!
      jmol    = 0
      kloc    = 1
      locj(1) = 1
      mormol  = .true.
!
      do while (mormol)
         jmol = jmol + 1
!
         if (jmol .gt. maxmol) then
            write (16, 10) ' jmol .gt. maxmol'
            mormol = .false.
!
         else
            read (12, 13) c, e1, e2, e3, e4, e5, e6, e7
13          format (f18.2, f7.3, 6e11.4)
!
            if (c .eq. 0.0d0) then
               mormol = .false.
!
            else
               write (16, 14) jmol, c, e1, e2, e3, e4, e5, e6, e7
14             format (i5, f18.2, f7.3, 1p6e11.4)
!
               ii = 1
!
               do while (ii .le. 8 .and. c .lt. xcode(ii))
                  ii = ii + 1
               end do
!
               if (ii .gt. 8) then
                  close (unit = 12)
                  stop ' in readmol ii .gt. 8'
               end if
!
               x = c
!
               do i = ii, 8
                  id = x / xcode(i) + 0.5
                  x  = x - dble(id) * xcode(i)
                  if (id .eq. 0) id = 100
                  ifequa (id) = 1
                  kcomps (kloc) = id
                  kloc = kloc + 1
               end do
!
               ion = x * 100.0d0 + 0.5
!
               if (ion .ge. 1) then
                  ifequa(100) = 1
                  ifequa(101) = 1
!
                  do i = 1, ion
                     kcomps(kloc) = 101
                     kloc = kloc + 1
                  end do
!
                  if (kloc .gt. maxloc) then
                     write (16, 10) ' kloc .gt. maxloc'
                     mormol = .false.
                  end if
!
               end if
!
               locj (jmol + 1) = kloc
               code (jmol) = c
               equil (1, jmol) = e1
               equil (2, jmol) = e2
               equil (3, jmol) = e3
               equil (4, jmol) = e4
               equil (5, jmol) = e5
               equil (6, jmol) = e6
               equil (7, jmol) = e7
            end if
!
         end if
!
      end do
!
      nummol = jmol - 1
      nloc   = kloc - 1
!
!.... ASSIGN AN EQUATION NUMBER TO EACH COMPONENT
!.... THE FIRST EQUATION IS FOR THE TOTAL NUMBER OF PARTICLES
!.... THE FIRST VARIABLE IS XNATOM
!.... IF ANY COMPONENT IS 100 OR 101 VARIABLE NEQUA IS XNE
!....     AND EQUATION NEQUA IS CHARGE CONSERVATION
!.... DIMENSIONS ARE SET FOR A MAXIMUM 25 EQUATIONS, IN COMMON.SIZEBL
!
      iequa = 1
!
      do i = 1, 100
!
         if (ifequa(i) .gt. 0) then
            iequa = iequa + 1
            ifequa(i) = iequa
            idequa(iequa) = i
         end if
!
      end do
!
      nequa  = iequa
      ifequa(101) = nequa + 1
!
      do kloc = 1, nloc
         id = kcomps(kloc)
         kcomps(kloc) = ifequa(id)
      end do
!
      write (16, 29) 
     &   ' number of molecules used  = ', nummol, ' max = ', maxmol,
     &   ' number of components used = ', nloc,   ' max = ', maxloc,
     &   ' number of equations used  = ', nequa,  ' max = ', maxmeq
   29 format (a, i4, a, i4 /
     &        a, i4, a, i4 /
     &        a, i4, a, i4)
      close (unit = 12)


      return 

7600  continue
      print*, 'Could not find molecules.dat file= ', inputmolec
      print*, 'Going to abort ...'
      stop 
  


      end
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   R E A D M O L !!!!!!!!!!
!
      double precision function rosstab (t, p, v)

      implicit none
!     CHECK if rosstab still works with the local delaration !
!.... ROSSELAND OPACITY FOR SDSC 1X DISTRIBUTION FUNCTIONS
!.... VERSION OF 11 OCTOBER 1992
!
!----------------------------- DUMMY VARIABLES -------------------------
!
      double precision p, t, v
!
!------------------------------- CONSTANTS -----------------------------
!
      integer nnp, nnt, nnv
      parameter (nnp = 38, nnt = 56, nnv = 5)

      include 'common.tabptv'
!
!
!-------------------------- LOCAL VARIABLES --------------------------
!
!*....
      double precision  tabp(38), tabt(56)
!*
!***************** E N D  O F   C O M M O N   T A B P T V *****************
!*



      double precision  ablog, exp10, pkap(nnp), pkapv(nnv, nnp), tlog, 
     &                  tsave, vsave, x
      integer idum, ip, it, iv, ktab(nnt, nnp, nnv), map1
!
!---------------------------- external ---------------------------------
!
      external map1
!
!----------------------------- INITIALIZATION --------------------------
!
!
      data (ktab(it, 1, 1), it = 1, nnt) /                              -2.000
     &  -5201,-5133,-5045,-4888,-4694,-4508,-4360,-4257,-4197,-4170,    -2.000
     &  -4154,-4122,-4029,-3840,-3576,-3270,-2942,-2604,-2265,-1932,    -2.000
     &  -1450,-1011, -685, -549, -519, -520, -526, -531, -527, -513,    -2.000
     &   -492, -482, -482, -487, -482, -472, -448, -440, -435, -434,    -2.000
     &   -435, -431, -428, -428, -429, -428, -428, -430, -433, -434,    -2.000
     &   -432, -427, -418, -405, -386, -360/                            -2.000
!
      data (ktab(it, 2, 1), it = 1, nnt) /                              -1.500
     &  -5234,-5146,-5078,-4967,-4796,-4606,-4433,-4302,-4213,-4164,    -1.500
     &  -4144,-4125,-4076,-3947,-3728,-3450,-3142,-2818,-2490,-2162,    -1.500
     &  -1681,-1227, -834, -592, -515, -504, -511, -520, -525, -518,    -1.500
     &   -490, -478, -473, -480, -479, -470, -451, -437, -432, -429,    -1.500
     &   -432, -429, -426, -425, -427, -427, -425, -427, -432, -433,    -1.500
     &   -432, -427, -418, -405, -386, -360/                            -1.500
!
      data (ktab(it, 3, 1), it = 1, nnt) /                              -1.000
     &  -5278,-5167,-5089,-5010,-4878,-4701,-4517,-4358,-4241,-4164,    -1.000
     &  -4127,-4111,-4087,-4014,-3851,-3611,-3326,-3020,-2703,-2383,    -1.000
     &  -1907,-1446,-1019, -682, -517, -479, -483, -497, -509, -513,    -1.000
     &   -489, -474, -463, -469, -473, -467, -453, -433, -429, -424,    -1.000
     &   -427, -427, -424, -422, -424, -425, -423, -424, -430, -433,    -1.000
     &   -432, -427, -418, -405, -386, -360/                            -1.000
!
      data (ktab(it, 4, 1), it = 1, nnt) /                              -0.500
     &  -5322,-5193,-5093,-5019,-4924,-4775,-4594,-4417,-4271,-4167,    -0.500
     &  -4104,-4078,-4065,-4032,-3931,-3740,-3487,-3201,-2899,-2589,    -0.500
     &  -2121,-1658,-1212, -815, -542, -439, -429, -449, -472, -489,    -0.500
     &   -485, -466, -453, -453, -463, -460, -449, -430, -424, -419,    -0.500
     &   -421, -423, -422, -419, -420, -422, -422, -421, -426, -431,    -0.500
     &   -431, -426, -418, -405, -386, -360/                            -0.500
!
      data (ktab(it, 5, 1), it = 1, nnt) /                               0.000
     &  -5347,-5214,-5093,-5002,-4923,-4808,-4644,-4459,-4288,-4154,     0.000
     &  -4063,-4015,-4002,-3993,-3947,-3818,-3606,-3347,-3063,-2767,     0.000
     &  -2311,-1850,-1394, -960, -599, -392, -339, -356, -393, -429,     0.000
     &   -461, -449, -440, -433, -444, -448, -440, -425, -416, -413,     0.000
     &   -413, -418, -418, -416, -415, -419, -419, -417, -421, -428,     0.000
     &   -430, -426, -418, -405, -386, -360/                             0.000
!
      data (ktab(it, 6, 1), it = 1, nnt) /                               0.500
     &  -5350,-5206,-5071,-4959,-4872,-4779,-4641,-4462,-4274,-4108,     0.500
     &  -3985,-3908,-3878,-3879,-3872,-3808,-3654,-3430,-3170,-2892,     0.500
     &  -2455,-2002,-1542,-1086, -667, -352, -210, -200, -244, -303,     0.500
     &   -388, -408, -412, -406, -416, -426, -424, -413, -404, -403,     0.500
     &   -402, -409, -413, -412, -409, -413, -416, -414, -415, -424,     0.500
     &   -428, -426, -418, -405, -386, -360/                             0.500
!
      data (ktab(it, 7, 1), it = 1, nnt) /                               1.000
     &  -5342,-5161,-5008,-4874,-4767,-4677,-4564,-4403,-4212,-4021,     1.000
     &  -3862,-3748,-3686,-3673,-3684,-3674,-3590,-3415,-3184,-2926,     1.000
     &  -2515,-2080,-1628,-1167, -716, -323,  -63,   18,  -11,  -84,     1.000
     &   -228, -313, -346, -360, -372, -392, -395, -389, -382, -387,     1.000
     &   -388, -396, -405, -407, -405, -406, -411, -410, -409, -418,     1.000
     &   -425, -424, -417, -405, -386, -360/                             1.000
!
      data (ktab(it, 8, 1), it = 1, nnt) /                               1.250
     &  -5336,-5130,-4960,-4812,-4692,-4595,-4491,-4343,-4156,-3959,     1.250
     &  -3784,-3651,-3569,-3540,-3549,-3555,-3505,-3361,-3148,-2902,     1.250
     &  -2504,-2083,-1640,-1182, -726, -307,    4,  142,  134,   61,     1.250
     &   -106, -232, -289, -322, -341, -366, -374, -371, -366, -373,     1.250
     &   -378, -387, -398, -403, -402, -402, -408, -408, -406, -414,     1.250
     &   -423, -423, -417, -405, -386, -360/                             1.250
!
      data (ktab(it, 9, 1), it = 1, nnt) /                               1.500
     &  -5330,-5098,-4903,-4739,-4602,-4495,-4394,-4259,-4080,-3882,     1.500
     &  -3694,-3543,-3440,-3391,-3390,-3406,-3385,-3276,-3085,-2851,     1.500
     &  -2467,-2059,-1628,-1179, -722, -288,   66,  266,  296,  230,     1.500
     &     44, -124, -211, -268, -297, -333, -349, -347, -344, -354,     1.500
     &   -365, -376, -391, -398, -398, -397, -404, -406, -403, -409,     1.500
     &   -420, -422, -417, -405, -386, -360/                             1.500
!
      data (ktab(it, 10, 1), it = 1, nnt) /                              1.750
     &  -5319,-5067,-4843,-4657,-4501,-4379,-4274,-4149,-3983,-3787,     1.750
     &  -3592,-3425,-3304,-3234,-3217,-3232,-3232,-3158,-2995,-2776,     1.750
     &  -2406,-2011,-1594,-1157, -706, -264,  121,  379,  461,  415,     1.750
     &    222,   15, -109, -193, -241, -288, -316, -318, -316, -327,     1.750
     &   -347, -362, -380, -391, -394, -392, -398, -403, -400, -404,     1.750
     &   -416, -420, -416, -405, -386, -360/                             1.750
!
      data (ktab(it, 11, 1), it = 1, nnt) /                              2.000
     &  -5305,-5039,-4785,-4571,-4392,-4249,-4137,-4018,-3865,-3677,     2.000
     &  -3478,-3299,-3160,-3070,-3034,-3040,-3052,-3012,-2880,-2680,     2.000
     &  -2326,-1943,-1541,-1117, -676, -233,  172,  479,  621,  609,     2.000
     &    423,  183,   15,  -96, -168, -228, -273, -283, -280, -291,     2.000
     &   -321, -342, -366, -382, -388, -386, -392, -399, -397, -398,     2.000
     &   -412, -418, -415, -404, -386, -360/                             2.000
!
      data (ktab(it, 12, 1), it = 1, nnt) /                              2.250
     &  -5286,-5010,-4731,-4484,-4281,-4115,-3987,-3868,-3727,-3550,     2.250
     &  -3354,-3167,-3013,-2905,-2848,-2840,-2854,-2840,-2743,-2566,     2.250
     &  -2229,-1859,-1470,-1062, -634, -196,  219,  562,  764,  799,     2.250
     &    636,  379,  168,   27,  -75, -151, -214, -237, -237, -246,     2.250
     &   -283, -315, -346, -369, -380, -381, -384, -393, -393, -393,     2.250
     &   -406, -415, -414, -404, -386, -360/                             2.250
!
      data (ktab(it, 13, 1), it = 1, nnt) /                              2.500
     &  -5261,-4981,-4682,-4404,-4171,-3980,-3831,-3707,-3574,-3410,     2.500
     &  -3221,-3030,-2864,-2738,-2662,-2636,-2644,-2647,-2585,-2436,     2.500
     &  -2121,-1764,-1388, -994, -580, -152,  266,  633,  887,  980,     2.500
     &    857,  595,  344,  173,   39,  -55, -140, -180, -186, -191,     2.500
     &   -233, -278, -319, -352, -369, -373, -375, -386, -388, -387,     2.500
     &   -399, -412, -412, -403, -385, -360/                             2.500
!
      data (ktab(it, 14, 1), it = 1, nnt) /                              2.750
     &  -5229,-4950,-4637,-4331,-4068,-3851,-3678,-3540,-3409,-3258,     2.750
     &  -3079,-2889,-2714,-2574,-2480,-2435,-2431,-2440,-2407,-2291,     2.750
     &  -2003,-1659,-1295, -915, -517, -101,  313,  694,  987, 1138,     2.750
     &   1076,  824,  546,  344,  179,   60,  -43, -106, -123, -127,     2.750
     &   -166, -226, -281, -326, -353, -363, -365, -377, -382, -381,     2.750
     &   -391, -407, -410, -402, -385, -360/                             2.750
!
      data (ktab(it, 15, 1), it = 1, nnt) /                              3.000
     &  -5190,-4915,-4594,-4268,-3974,-3729,-3530,-3373,-3238,-3096,     3.000
     &  -2931,-2746,-2566,-2413,-2302,-2237,-2218,-2225,-2214,-2131,     3.000
     &  -1876,-1548,-1196, -829, -445,  -44,  363,  750, 1071, 1275,     3.000
     &   1289, 1059,  765,  532,  342,  197,   73,  -14,  -47,  -53,     3.000
     &    -84, -158, -228, -291, -331, -350, -353, -365, -375, -374,     3.000
     &   -382, -400, -407, -401, -384, -359/                             3.000
!
      data (ktab(it, 16, 1), it = 1, nnt) /                              3.200
     &  -5152,-4882,-4559,-4221,-3907,-3640,-3420,-3245,-3101,-2962,     3.200
     &  -2806,-2627,-2447,-2287,-2163,-2085,-2051,-2051,-2049,-1992,     3.200
     &  -1767,-1454,-1111, -754, -382,    6,  406,  792, 1128, 1364,     3.200
     &   1447, 1246,  952,  697,  493,  324,  185,   77,   27,   14,     3.200
     &     -9,  -88, -173, -251, -305, -334, -341, -353, -366, -367,     3.200
     &   -374, -394, -403, -399, -384, -359/                             3.200
!
      data (ktab(it, 17, 1), it = 1, nnt) /                              3.400
     &  -5109,-4844,-4523,-4178,-3847,-3559,-3318,-3122,-2965,-2825,     3.400
     &  -2677,-2507,-2328,-2163,-2028,-1936,-1888,-1877,-1879,-1844,     3.400
     &  -1653,-1357,-1024, -677, -315,   61,  451,  834, 1178, 1441,     3.400
     &   1595, 1431, 1145,  870,  653,  465,  311,  184,  115,   92,     3.400
     &     74,   -6, -104, -200, -271, -313, -327, -338, -356, -360,     3.400
     &   -365, -386, -399, -397, -383, -359/                             3.400
!
      data (ktab(it, 18, 1), it = 1, nnt) /                              3.600
     &  -5060,-4802,-4483,-4136,-3794,-3486,-3224,-3008,-2835,-2689,     3.600
     &  -2545,-2384,-2208,-2040,-1897,-1791,-1729,-1707,-1707,-1688,     3.600
     &  -1533,-1256, -934, -596, -245,  120,  498,  876, 1224, 1507,     3.600
     &   1726, 1611, 1339, 1053,  825,  621,  452,  308,  218,  182,     3.600
     &    165,   90,  -21, -135, -227, -284, -309, -321, -342, -350,     3.600
     &   -355, -376, -393, -394, -381, -358/                             3.600
!
      data (ktab(it, 19, 1), it = 1, nnt) /                              3.800
     &  -5005,-4753,-4440,-4093,-3744,-3420,-3139,-2903,-2711,-2554,     3.800
     &  -2412,-2258,-2088,-1918,-1768,-1651,-1575,-1540,-1533,-1524,     3.800
     &  -1404,-1151, -842, -513, -172,  182,  549,  919, 1268, 1563,     3.800
     &   1840, 1784, 1535, 1243, 1004,  790,  608,  450,  339,  286,     3.800
     &    265,  200,   76,  -53, -168, -245, -285, -301, -324, -338,     3.800
     &   -343, -364, -386, -390, -380, -357/                             3.800
!
      data (ktab(it, 20, 1), it = 1, nnt) /                              4.000
     &  -4943,-4700,-4392,-4050,-3698,-3362,-3063,-2808,-2596,-2425,     4.000
     &  -2279,-2130,-1966,-1798,-1642,-1515,-1426,-1377,-1362,-1356,     4.000
     &  -1270,-1042, -746, -427,  -95,  247,  603,  964, 1310, 1614,     4.000
     &   1939, 1946, 1726, 1437, 1187,  967,  774,  605,  475,  404,     4.000
     &    376,  320,  190,   43,  -95, -196, -253, -276, -303, -323,     4.000
     &   -330, -350, -377, -386, -377, -356/                             4.000
!
      data (ktab(it, 21, 1), it = 1, nnt) /                              4.200
     &  -4886,-4650,-4348,-4011,-3658,-3313,-2999,-2726,-2495,-2306,     4.200
     &  -2150,-2003,-1844,-1677,-1519,-1384,-1282,-1220,-1193,-1185,     4.200
     &  -1125, -927, -647, -338,  -16,  316,  661, 1011, 1352, 1659,     4.200
     &   2021, 2094, 1913, 1635, 1377, 1155,  953,  776,  632,  542,     4.200
     &    500,  452,  320,  157,   -3, -130, -209, -245, -275, -302,     4.200
     &   -314, -334, -365, -379, -373, -354/                             4.200
!
      data (ktab(it, 22, 1), it = 1, nnt) /                              4.400
     &  -4822,-4596,-4299,-3969,-3619,-3269,-2942,-2652,-2403,-2195,     4.400
     &  -2026,-1877,-1722,-1558,-1397,-1255,-1143,-1067,-1028,-1015,     4.400
     &   -974, -807, -545, -247,   65,  387,  721, 1061, 1395, 1702,     4.400
     &   2092, 2227, 2094, 1834, 1569, 1346, 1138,  956,  800,  693,     4.400
     &    637,  593,  464,  288,  105,  -49, -153, -206, -241, -276,     4.400
     &   -294, -314, -350, -371, -369, -351/                             4.400
!
      data (ktab(it, 23, 1), it = 1, nnt) /                              4.600
     &  -4753,-4539,-4247,-3925,-3580,-3229,-2893,-2588,-2321,-2096,     4.600
     &  -1911,-1754,-1601,-1439,-1278,-1130,-1009, -920, -868, -846,     4.600
     &   -816, -680, -438, -152,  149,  462,  784, 1113, 1440, 1744,     4.600
     &   2152, 2343, 2267, 2031, 1766, 1541, 1331, 1145,  980,  860,     4.600
     &    788,  741,  622,  435,  232,   50,  -82, -156, -200, -243,     4.600
     &   -268, -290, -331, -359, -362, -347/                             4.600
!
      data (ktab(it, 24, 1), it = 1, nnt) /                              4.800
     &  -4679,-4481,-4193,-3879,-3541,-3192,-2850,-2532,-2250,-2007,     4.800
     &  -1805,-1637,-1481,-1321,-1159,-1008, -878, -779, -713, -680,     4.800
     &   -652, -546, -327,  -55,  237,  539,  850, 1168, 1486, 1786,     4.800
     &   2205, 2443, 2428, 2225, 1967, 1739, 1530, 1340, 1171, 1039,     4.800
     &    953,  898,  790,  596,  377,  169,    6,  -93, -150, -201,     4.800
     &   -236, -262, -306, -344, -353, -342/                             4.800
!
      data (ktab(it, 25, 1), it = 1, nnt) /                              5.000
     &  -4596,-4421,-4135,-3828,-3498,-3154,-2810,-2483,-2187,-1928,     5.000
     &  -1709,-1526,-1365,-1205,-1043, -889, -752, -642, -564, -519,     5.000
     &   -487, -406, -212,   45,  327,  619,  919, 1226, 1535, 1828,     5.000
     &   2254, 2529, 2575, 2413, 2168, 1937, 1730, 1539, 1367, 1224,     5.000
     &   1126, 1061,  964,  767,  536,  304,  112,  -17,  -90, -150,     5.000
     &   -197, -228, -275, -324, -341, -335/                             5.000
!
      data (ktab(it, 26, 1), it = 1, nnt) /                              5.200
     &  -4505,-4358,-4074,-3774,-3452,-3114,-2771,-2439,-2132,-1859,     5.200
     &  -1624,-1426,-1255,-1092, -929, -772, -630, -511, -421, -364,     5.200
     &   -323, -258,  -90,  151,  421,  703,  992, 1288, 1586, 1872,     5.200
     &   2299, 2601, 2705, 2594, 2371, 2140, 1937, 1745, 1572, 1424,     5.200
     &   1312, 1234, 1142,  950,  708,  459,  240,   79,  -16,  -88,     5.200
     &   -146, -186, -237, -295, -324, -324/                             5.200
!
      data (ktab(it, 27, 1), it = 1, nnt) /                              5.400
     &  -4403,-4293,-4011,-3715,-3403,-3072,-2733,-2399,-2083,-1798,     5.400
     &  -1548,-1334,-1151, -983, -818, -658, -511, -384, -284, -214,     5.400
     &   -159, -105,   36,  259,  518,  789, 1068, 1353, 1641, 1918,     5.400
     &   2342, 2664, 2819, 2764, 2569, 2344, 2143, 1952, 1778, 1625,     5.400
     &   1503, 1412, 1322, 1138,  889,  626,  383,  193,   70,  -14,     5.400
     &    -85, -135, -190, -260, -301, -309/                             5.400
!
      data (ktab(it, 28, 1), it = 1, nnt) /                              5.600
     &  -4287,-4226,-3945,-3654,-3349,-3027,-2694,-2360,-2039,-1743,     5.600
     &  -1480,-1252,-1056, -880, -712, -549, -397, -264, -154,  -72,     5.600
     &     -2,   48,  168,  372,  617,  878, 1146, 1421, 1698, 1967,     5.600
     &   2385, 2721, 2917, 2920, 2762, 2549, 2349, 2162, 1988, 1832,     5.600
     &   1701, 1598, 1503, 1330, 1076,  804,  542,  326,  174,   73,     5.600
     &    -11,  -74, -134, -213, -270, -289/                             5.600
!
      data (ktab(it, 29, 1), it = 1, nnt) /                              5.800
     &  -4156,-4157,-3878,-3590,-3291,-2978,-2652,-2321,-1998,-1695,     5.800
     &  -1421,-1180, -970, -785, -612, -445, -289, -149,  -31,   61,     5.800
     &    147,  202,  303,  488,  720,  970, 1228, 1492, 1759, 2019,     5.800
     &   2429, 2772, 3001, 3060, 2948, 2754, 2557, 2373, 2200, 2043,     5.800
     &   1906, 1791, 1686, 1523, 1269,  990,  716,  477,  297,  175,     5.800
     &     76,    0,  -68, -155, -229, -262/                             5.800
!
      data (ktab(it, 30, 1), it = 1, nnt) /                              6.000
     &  -4006,-4087,-3809,-3524,-3231,-2926,-2607,-2282,-1959,-1652,     6.000
     &  -1369,-1116, -893, -698, -518, -347, -186,  -39,   86,  189,     6.000
     &    293,  355,  442,  607,  826, 1065, 1314, 1567, 1824, 2073,     6.000
     &   2474, 2821, 3075, 3185, 3123, 2954, 2763, 2584, 2412, 2255,     6.000
     &   2111, 1986, 1871, 1715, 1465, 1183,  898,  641,  435,  292,     6.000
     &    178,   85,    9,  -85, -177, -227/                             6.000
!
      data (ktab(it, 31, 1), it = 1, nnt) /                              6.250
     &  -3785,-3999,-3721,-3439,-3151,-2855,-2547,-2230,-1911,-1602,     6.250
     &  -1311,-1046, -810, -601, -412, -234,  -66,   87,  222,  336,     6.250
     &    459,  534,  611,  753,  955, 1182, 1420, 1662, 1907, 2145,     6.250
     &   2533, 2880, 3155, 3315, 3321, 3196, 3025, 2852, 2684, 2524,     6.250
     &   2378, 2242, 2111, 1956, 1713, 1425, 1135,  865,  636,  464,     6.250
     &    330,  217,  123,   20,  -90, -166/                             6.250
!
      data (ktab(it, 32, 1), it = 1, nnt) /                              6.500
     &  -3518,-3908,-3631,-3351,-3068,-2779,-2481,-2174,-1862,-1554,     6.500
     &  -1259, -986, -738, -516, -316, -132,   41,  203,  348,  472,     6.500
     &    614,  707,  782,  904, 1087, 1303, 1531, 1762, 1995, 2222,     6.500
     &   2595, 2939, 3227, 3428, 3496, 3426, 3279, 3117, 2954, 2793,     6.500
     &   2643, 2499, 2355, 2195, 1962, 1671, 1377, 1100,  853,  655,     6.500
     &    503,  370,  259,  148,   19,  -84/                             6.500
!
      data (ktab(it, 33, 1), it = 1, nnt) /                              6.750
     &  -3195,-3818,-3541,-3263,-2983,-2700,-2410,-2113,-1810,-1506,     6.750
     &  -1212, -933, -676, -443, -233,  -41,  138,  308,  462,  597,     6.750
     &    757,  868,  949, 1056, 1221, 1423, 1640, 1861, 2084, 2300,     6.750
     &   2659, 2999, 3295, 3522, 3643, 3633, 3525, 3384, 3230, 3071,     6.750
     &   2915, 2761, 2604, 2436, 2210, 1919, 1619, 1334, 1075,  867,     6.750
     &    699,  552,  420,  296,  155,   24/                             6.750
!
      data (ktab(it, 34, 1), it = 1, nnt) /                              7.000
     &  -2989,-3727,-3450,-3174,-2897,-2618,-2336,-2048,-1754,-1457,     7.000
     &  -1166, -885, -622, -380, -160,   39,  225,  400,  563,  708,     7.000
     &    887, 1018, 1111, 1211, 1358, 1546, 1752, 1964, 2176, 2384,     7.000
     &   2729, 3061, 3361, 3607, 3771, 3820, 3758, 3640, 3501, 3349,     7.000
     &   3191, 3028, 2857, 2666, 2447, 2168, 1867, 1581, 1315, 1085,     7.000
     &    899,  746,  602,  466,  315,  159/                             7.000
!
      data (ktab(it, 35, 1), it = 1, nnt) /                              7.250
     &  -2880,-3636,-3360,-3084,-2810,-2535,-2258,-1978,-1693,-1406,     7.250
     &  -1119, -840, -574, -326,  -98,  109,  300,  480,  650,  805,     7.250
     &   1002, 1152, 1262, 1362, 1497, 1673, 1868, 2071, 2274, 2472,     7.250
     &   2803, 3128, 3428, 3686, 3880, 3978, 3969, 3885, 3770, 3629,     7.250
     &   3472, 3303, 3120, 2914, 2688, 2411, 2108, 1821, 1556, 1320,     7.250
     &   1124,  959,  799,  647,  494,  322/                             7.250
!
      data (ktab(it, 36, 1), it = 1, nnt) /                              7.500
     &  -2775,-3545,-3270,-2995,-2722,-2451,-2179,-1906,-1629,-1351,     7.500
     &  -1072, -796, -531, -279,  -45,  169,  367,  551,  725,  887,     7.500
     &   1100, 1272, 1401, 1509, 1636, 1801, 1987, 2181, 2375, 2565,     7.500
     &   2884, 3200, 3498, 3764, 3979, 4120, 4162, 4117, 4025, 3902,     7.500
     &   3751, 3578, 3383, 3163, 2929, 2658, 2354, 2066, 1800, 1561,     7.500
     &   1356, 1186, 1017,  853,  688,  502/                             7.500
!
      data (ktab(it, 37, 1), it = 1, nnt) /                              7.750
     &  -2674,-3454,-3180,-2907,-2636,-2367,-2099,-1831,-1563,-1292,     7.750
     &  -1021, -752, -489, -236,    1,  221,  423,  612,  789,  955,     7.750
     &   1181, 1371, 1521, 1643, 1770, 1928, 2107, 2294, 2481, 2664,     7.750
     &   2971, 3277, 3571, 3841, 4071, 4240, 4327, 4326, 4265, 4166,     7.750
     &   4027, 3853, 3650, 3416, 3169, 2898, 2596, 2305, 2040, 1801,     7.750
     &   1592, 1416, 1244, 1071,  897,  704/                             7.750
!
      data (ktab(it, 38, 1), it = 1, nnt) /                              8.000
     &  -2576,-3364,-3092,-2820,-2550,-2283,-2019,-1756,-1494,-1232,     8.000
     &   -969, -706, -448, -197,   41,  265,  472,  664,  844, 1013,     8.000
     &   1247, 1453, 1624, 1764, 1898, 2053, 2228, 2410, 2589, 2766,     8.000
     &   3063, 3361, 3650, 3920, 4160, 4353, 4478, 4520, 4492, 4415,     8.000
     &   4293, 4124, 3914, 3669, 3410, 3136, 2840, 2548, 2282, 2043,     8.000
     &   1831, 1650, 1476, 1297, 1116,  920/                             8.000
!
      data (ktab(it, 1, 2), it = 1, nnt) /                              -2.000
     &  -5189,-5123,-5036,-4880,-4687,-4501,-4354,-4252,-4192,-4166,    -2.000
     &  -4151,-4120,-4027,-3838,-3574,-3268,-2940,-2602,-2264,-1930,    -2.000
     &  -1449,-1009, -684, -547, -517, -518, -525, -530, -526, -512,    -2.000
     &   -491, -480, -481, -486, -482, -472, -448, -440, -435, -434,    -2.000
     &   -434, -431, -428, -428, -429, -428, -427, -430, -433, -434,    -2.000
     &   -432, -427, -418, -405, -386, -360/                            -2.000
!
      data (ktab(it, 2, 2), it = 1, nnt) /                              -1.500
     &  -5220,-5134,-5068,-4957,-4788,-4598,-4426,-4295,-4207,-4159,    -1.500
     &  -4139,-4122,-4074,-3946,-3726,-3448,-3140,-2816,-2487,-2160,    -1.500
     &  -1680,-1225, -832, -590, -512, -502, -509, -518, -523, -517,    -1.500
     &   -489, -477, -471, -479, -479, -470, -451, -437, -432, -429,    -1.500
     &   -432, -429, -426, -425, -427, -427, -425, -427, -432, -433,    -1.500
     &   -432, -427, -418, -405, -386, -360/                            -1.500
!
      data (ktab(it, 3, 2), it = 1, nnt) /                              -1.000
     &  -5263,-5152,-5076,-5000,-4869,-4693,-4509,-4350,-4233,-4158,    -1.000
     &  -4121,-4107,-4084,-4012,-3849,-3609,-3324,-3017,-2700,-2380,    -1.000
     &  -1905,-1444,-1017, -680, -515, -476, -480, -495, -507, -511,    -1.000
     &   -488, -473, -461, -468, -473, -466, -452, -433, -429, -424,    -1.000
     &   -427, -427, -424, -421, -424, -425, -423, -424, -429, -432,    -1.000
     &   -431, -426, -418, -405, -386, -360/                            -1.000
!
      data (ktab(it, 4, 2), it = 1, nnt) /                              -0.500
     &  -5306,-5177,-5079,-5007,-4914,-4765,-4585,-4408,-4262,-4159,    -0.500
     &  -4096,-4072,-4061,-4029,-3928,-3738,-3484,-3198,-2895,-2586,    -0.500
     &  -2118,-1656,-1210, -813, -540, -437, -427, -446, -470, -487,    -0.500
     &   -483, -465, -451, -452, -462, -459, -448, -429, -423, -419,    -0.500
     &   -421, -423, -422, -419, -420, -422, -421, -421, -426, -431,    -0.500
     &   -431, -426, -418, -405, -386, -360/                            -0.500
!
      data (ktab(it, 5, 2), it = 1, nnt) /                               0.000
     &  -5332,-5197,-5076,-4988,-4910,-4796,-4633,-4449,-4278,-4144,     0.000
     &  -4053,-4007,-3996,-3988,-3943,-3815,-3603,-3344,-3060,-2763,     0.000
     &  -2308,-1847,-1392, -958, -597, -390, -336, -352, -389, -426,     0.000
     &   -459, -447, -438, -431, -443, -447, -439, -424, -415, -412,     0.000
     &   -412, -418, -418, -416, -414, -418, -419, -417, -421, -428,     0.000
     &   -430, -426, -418, -405, -386, -360/                             0.000
!
      data (ktab(it, 6, 2), it = 1, nnt) /                               0.500
     &  -5335,-5189,-5053,-4942,-4857,-4766,-4628,-4450,-4262,-4097,     0.500
     &  -3974,-3897,-3869,-3872,-3866,-3804,-3651,-3426,-3166,-2887,     0.500
     &  -2450,-1998,-1539,-1084, -665, -350, -208, -197, -240, -299,     0.500
     &   -385, -406, -409, -403, -414, -425, -423, -412, -403, -403,     0.500
     &   -402, -409, -413, -412, -409, -413, -416, -414, -415, -424,     0.500
     &   -428, -425, -418, -405, -386, -360/                             0.500
!
      data (ktab(it, 7, 2), it = 1, nnt) /                               1.000
     &  -5327,-5144,-4990,-4855,-4750,-4662,-4550,-4390,-4199,-4009,     1.000
     &  -3849,-3736,-3674,-3663,-3677,-3669,-3585,-3410,-3179,-2921,     1.000
     &  -2509,-2076,-1624,-1163, -714, -320,  -61,   21,   -8,  -80,     1.000
     &   -225, -310, -343, -357, -370, -390, -394, -388, -381, -386,     1.000
     &   -388, -396, -404, -407, -404, -406, -411, -410, -409, -418,     1.000
     &   -425, -424, -417, -405, -386, -360/                             1.000
!
      data (ktab(it, 8, 2), it = 1, nnt) /                               1.250
     &  -5321,-5113,-4941,-4793,-4674,-4580,-4476,-4328,-4142,-3946,     1.250
     &  -3771,-3639,-3557,-3530,-3540,-3549,-3500,-3356,-3143,-2897,     1.250
     &  -2498,-2077,-1636,-1179, -723, -305,    6,  144,  138,   65,     1.250
     &   -102, -228, -286, -319, -338, -364, -373, -369, -365, -372,     1.250
     &   -378, -387, -398, -402, -401, -401, -407, -408, -406, -414,     1.250
     &   -423, -423, -417, -405, -386, -360/                             1.250
!
      data (ktab(it, 9, 2), it = 1, nnt) /                               1.500
     &  -5315,-5081,-4885,-4719,-4583,-4479,-4379,-4245,-4066,-3868,     1.500
     &  -3681,-3530,-3428,-3380,-3381,-3399,-3379,-3270,-3080,-2845,     1.500
     &  -2461,-2053,-1624,-1175, -719, -285,   68,  268,  299,  234,     1.500
     &     48, -120, -208, -264, -294, -331, -348, -346, -342, -353,     1.500
     &   -365, -375, -390, -397, -398, -396, -403, -406, -403, -409,     1.500
     &   -420, -422, -417, -405, -386, -360/                             1.500
!
      data (ktab(it, 10, 2), it = 1, nnt) /                              1.750
     &  -5304,-5050,-4825,-4638,-4482,-4363,-4259,-4135,-3969,-3774,     1.750
     &  -3579,-3412,-3292,-3223,-3208,-3224,-3225,-3153,-2989,-2770,     1.750
     &  -2400,-2005,-1589,-1153, -703, -261,  123,  381,  464,  418,     1.750
     &    226,   19, -105, -189, -238, -285, -314, -316, -314, -326,     1.750
     &   -346, -361, -379, -390, -394, -392, -398, -402, -400, -404,     1.750
     &   -416, -420, -416, -405, -386, -360/                             1.750
!
      data (ktab(it, 11, 2), it = 1, nnt) /                              2.000
     &  -5291,-5022,-4767,-4551,-4372,-4233,-4122,-4005,-3851,-3664,     2.000
     &  -3466,-3287,-3148,-3059,-3024,-3032,-3046,-3006,-2874,-2674,     2.000
     &  -2319,-1937,-1535,-1113, -673, -230,  174,  481,  624,  612,     2.000
     &    427,  187,   19,  -92, -163, -225, -271, -281, -278, -290,     2.000
     &   -320, -341, -365, -382, -388, -386, -392, -398, -396, -398,     2.000
     &   -412, -418, -415, -404, -386, -360/                             2.000
!
      data (ktab(it, 12, 2), it = 1, nnt) /                              2.250
     &  -5272,-4994,-4713,-4466,-4263,-4099,-3973,-3856,-3714,-3538,     2.250
     &  -3343,-3156,-3001,-2894,-2839,-2831,-2847,-2833,-2736,-2559,     2.250
     &  -2223,-1853,-1465,-1058, -630, -194,  221,  564,  766,  802,     2.250
     &    640,  383,  172,   31,  -70, -148, -212, -235, -235, -244,     2.250
     &   -282, -314, -345, -369, -380, -380, -384, -393, -392, -393,     2.250
     &   -406, -415, -413, -404, -385, -360/                             2.250
!
      data (ktab(it, 13, 2), it = 1, nnt) /                              2.500
     &  -5248,-4965,-4665,-4386,-4154,-3964,-3818,-3695,-3562,-3399,     2.500
     &  -3211,-3020,-2853,-2728,-2653,-2627,-2637,-2641,-2578,-2430,     2.500
     &  -2114,-1757,-1382, -990, -577, -149,  268,  635,  889,  982,     2.500
     &    861,  599,  348,  177,   43,  -50, -137, -178, -183, -189,     2.500
     &   -231, -277, -318, -351, -369, -373, -375, -386, -388, -387,     2.500
     &   -399, -412, -412, -403, -385, -360/                             2.500
!
      data (ktab(it, 14, 2), it = 1, nnt) /                              2.750
     &  -5216,-4934,-4620,-4315,-4052,-3837,-3666,-3529,-3399,-3248,     2.750
     &  -3070,-2880,-2705,-2565,-2471,-2427,-2423,-2433,-2400,-2284,     2.750
     &  -1996,-1653,-1290, -911, -513,  -98,  316,  696,  990, 1141,     2.750
     &   1079,  828,  550,  348,  184,   64,  -40, -103, -120, -125,     2.750
     &   -164, -225, -280, -326, -353, -363, -365, -376, -382, -380,     2.750
     &   -391, -406, -409, -402, -385, -360/                             2.750
!
      data (ktab(it, 15, 2), it = 1, nnt) /                              3.000
     &  -5178,-4900,-4579,-4253,-3960,-3717,-3519,-3364,-3229,-3088,     3.000
     &  -2922,-2738,-2558,-2405,-2294,-2229,-2210,-2218,-2207,-2125,     3.000
     &  -1869,-1541,-1190, -824, -441,  -41,  366,  752, 1074, 1277,     3.000
     &   1292, 1063,  769,  535,  347,  201,   77,  -11,  -44,  -50,     3.000
     &    -82, -157, -227, -290, -330, -349, -352, -365, -375, -373,     3.000
     &   -381, -400, -406, -400, -384, -359/                             3.000
!
      data (ktab(it, 16, 2), it = 1, nnt) /                              3.200
     &  -5141,-4868,-4545,-4207,-3895,-3630,-3411,-3236,-3093,-2954,     3.200
     &  -2798,-2620,-2440,-2280,-2156,-2078,-2044,-2044,-2042,-1985,     3.200
     &  -1760,-1448,-1106, -750, -378,    9,  408,  794, 1130, 1366,     3.200
     &   1450, 1250,  956,  701,  497,  328,  188,   80,   30,   17,     3.200
     &     -7,  -87, -171, -250, -304, -333, -341, -352, -366, -367,     3.200
     &   -373, -393, -403, -399, -384, -359/                             3.200
!
      data (ktab(it, 17, 2), it = 1, nnt) /                              3.400
     &  -5099,-4832,-4510,-4166,-3837,-3550,-3310,-3115,-2959,-2818,     3.400
     &  -2670,-2501,-2322,-2156,-2022,-1929,-1881,-1871,-1873,-1837,     3.400
     &  -1646,-1351,-1019, -672, -311,   64,  453,  836, 1180, 1443,     3.400
     &   1597, 1434, 1148,  874,  657,  468,  315,  187,  117,   95,     3.400
     &     76,   -4, -103, -199, -271, -312, -327, -338, -355, -359,     3.400
     &   -364, -385, -399, -397, -383, -359/                             3.400
!
      data (ktab(it, 18, 2), it = 1, nnt) /                              3.600
     &  -5050,-4791,-4472,-4125,-3784,-3477,-3217,-3002,-2829,-2683,     3.600
     &  -2539,-2378,-2203,-2034,-1891,-1785,-1723,-1700,-1700,-1681,     3.600
     &  -1526,-1250, -929, -591, -241,  123,  501,  878, 1226, 1509,     3.600
     &   1728, 1614, 1343, 1056,  828,  624,  456,  311,  221,  185,     3.600
     &    167,   92,  -20, -133, -226, -284, -308, -321, -342, -350,     3.600
     &   -354, -376, -393, -394, -381, -358/                             3.600
!
      data (ktab(it, 19, 2), it = 1, nnt) /                              3.800
     &  -4996,-4744,-4430,-4084,-3736,-3413,-3133,-2898,-2706,-2550,     3.800
     &  -2407,-2253,-2083,-1913,-1763,-1646,-1570,-1534,-1527,-1518,     3.800
     &  -1398,-1145, -836, -508, -167,  185,  552,  921, 1270, 1565,     3.800
     &   1842, 1786, 1538, 1246, 1007,  793,  611,  453,  342,  289,     3.800
     &    267,  201,   78,  -52, -167, -245, -284, -300, -324, -338,     3.800
     &   -343, -364, -386, -390, -379, -357/                             3.800
!
      data (ktab(it, 20, 2), it = 1, nnt) /                              4.000
     &  -4937,-4692,-4384,-4042,-3691,-3356,-3058,-2804,-2592,-2421,     4.000
     &  -2274,-2126,-1962,-1793,-1638,-1511,-1421,-1372,-1356,-1350,     4.000
     &  -1263,-1036, -741, -422,  -91,  251,  606,  967, 1312, 1616,     4.000
     &   1941, 1948, 1728, 1439, 1189,  969,  777,  607,  477,  407,     4.000
     &    378,  322,  191,   44,  -94, -195, -252, -275, -302, -322,     4.000
     &   -329, -350, -377, -386, -377, -356/                             4.000
!
      data (ktab(it, 21, 2), it = 1, nnt) /                              4.200
     &  -4879,-4642,-4340,-4004,-3652,-3309,-2995,-2722,-2492,-2303,     4.200
     &  -2147,-2000,-1841,-1673,-1515,-1379,-1277,-1214,-1187,-1179,     4.200
     &  -1119, -921, -641, -333,  -12,  320,  663, 1014, 1355, 1661,     4.200
     &   2023, 2096, 1915, 1638, 1379, 1157,  955,  778,  634,  545,     4.200
     &    502,  454,  322,  159,   -2, -129, -208, -244, -274, -301,     4.200
     &   -313, -333, -365, -379, -373, -354/                             4.200
!
      data (ktab(it, 22, 2), it = 1, nnt) /                              4.400
     &  -4815,-4589,-4292,-3962,-3613,-3265,-2939,-2649,-2401,-2193,     4.400
     &  -2024,-1874,-1718,-1554,-1393,-1250,-1138,-1061,-1022,-1008,     4.400
     &   -967, -801, -539, -242,   69,  391,  724, 1063, 1397, 1704,     4.400
     &   2093, 2229, 2096, 1836, 1571, 1348, 1140,  958,  801,  695,     4.400
     &    639,  594,  466,  290,  106,  -48, -153, -205, -241, -275,     4.400
     &   -293, -314, -350, -371, -369, -351/                             4.400
!
      data (ktab(it, 23, 2), it = 1, nnt) /                              4.600
     &  -4747,-4533,-4240,-3919,-3575,-3226,-2890,-2586,-2320,-2094,     4.600
     &  -1909,-1752,-1597,-1435,-1274,-1126,-1004, -915, -863, -840,     4.600
     &   -809, -674, -433, -148,  153,  465,  787, 1116, 1442, 1746,     4.600
     &   2154, 2345, 2268, 2033, 1768, 1543, 1333, 1146,  982,  861,     4.600
     &    790,  743,  624,  437,  234,   51,  -81, -156, -200, -242,     4.600
     &   -268, -290, -331, -359, -362, -347/                             4.600
!
      data (ktab(it, 24, 2), it = 1, nnt) /                              4.800
     &  -4672,-4475,-4186,-3873,-3536,-3188,-2847,-2530,-2249,-2006,     4.800
     &  -1804,-1635,-1479,-1319,-1157,-1005, -875, -775, -709, -676,     4.800
     &   -648, -541, -322,  -50,  241,  543,  853, 1171, 1488, 1787,     4.800
     &   2207, 2444, 2429, 2226, 1968, 1741, 1532, 1341, 1173, 1040,     4.800
     &    955,  899,  792,  597,  378,  170,    7,  -93, -150, -200,     4.800
     &   -236, -262, -306, -343, -353, -342/                             4.800
!
      data (ktab(it, 25, 2), it = 1, nnt) /                              5.000
     &  -4590,-4414,-4129,-3822,-3494,-3150,-2807,-2481,-2185,-1927,     5.000
     &  -1707,-1524,-1363,-1203,-1041, -886, -750, -639, -561, -515,     5.000
     &   -484, -401, -207,   49,  331,  623,  922, 1229, 1537, 1830,     5.000
     &   2255, 2530, 2576, 2414, 2169, 1938, 1732, 1540, 1368, 1225,     5.000
     &   1127, 1063,  965,  768,  537,  305,  113,  -17,  -90, -150,     5.000
     &   -196, -227, -275, -323, -341, -335/                             5.000
!
      data (ktab(it, 26, 2), it = 1, nnt) /                              5.200
     &  -4499,-4351,-4069,-3768,-3448,-3111,-2769,-2437,-2130,-1857,     5.200
     &  -1622,-1424,-1253,-1090, -927, -770, -628, -509, -419, -361,     5.200
     &   -319, -254,  -86,  154,  424,  706,  995, 1290, 1588, 1874,     5.200
     &   2300, 2602, 2706, 2595, 2372, 2141, 1938, 1746, 1573, 1424,     5.200
     &   1313, 1235, 1143,  951,  709,  460,  240,   80,  -16,  -87,     5.200
     &   -145, -185, -237, -295, -324, -324/                             5.200
!
      data (ktab(it, 27, 2), it = 1, nnt) /                              5.400
     &  -4397,-4287,-4005,-3711,-3398,-3069,-2731,-2397,-2081,-1796,     5.400
     &  -1546,-1332,-1150, -982, -817, -657, -509, -382, -282, -211,     5.400
     &   -157, -102,   40,  263,  521,  792, 1070, 1355, 1643, 1920,     5.400
     &   2343, 2665, 2820, 2765, 2570, 2345, 2143, 1953, 1779, 1626,     5.400
     &   1504, 1413, 1323, 1139,  890,  627,  383,  194,   70,  -13,     5.400
     &    -85, -135, -189, -259, -301, -309/                             5.400
!
      data (ktab(it, 28, 2), it = 1, nnt) /                              5.600
     &  -4282,-4220,-3940,-3650,-3345,-3024,-2691,-2358,-2037,-1742,     5.600
     &  -1479,-1251,-1055, -879, -711, -548, -396, -262, -152,  -70,     5.600
     &      0,   50,  170,  375,  620,  881, 1149, 1423, 1700, 1969,     5.600
     &   2386, 2721, 2918, 2921, 2763, 2550, 2350, 2162, 1989, 1832,     5.600
     &   1702, 1599, 1503, 1331, 1077,  805,  543,  326,  175,   73,     5.600
     &    -11,  -73, -133, -213, -270, -289/                             5.600
!
      data (ktab(it, 29, 2), it = 1, nnt) /                              5.800
     &  -4151,-4152,-3873,-3586,-3288,-2975,-2650,-2319,-1996,-1694,     5.800
     &  -1420,-1179, -969, -784, -611, -444, -287, -147,  -29,   63,     5.800
     &    149,  204,  306,  491,  723,  973, 1231, 1494, 1761, 2020,     5.800
     &   2430, 2773, 3002, 3061, 2948, 2754, 2557, 2374, 2200, 2044,     5.800
     &   1906, 1791, 1686, 1524, 1270,  991,  716,  478,  298,  176,     5.800
     &     76,    0,  -67, -155, -229, -262/                             5.800
!
      data (ktab(it, 30, 2), it = 1, nnt) /                              6.000
     &  -4002,-4082,-3804,-3520,-3227,-2923,-2605,-2280,-1958,-1650,     6.000
     &  -1367,-1114, -892, -696, -517, -346, -184,  -38,   88,  190,     6.000
     &    295,  357,  444,  610,  828, 1068, 1316, 1569, 1825, 2075,     6.000
     &   2475, 2822, 3076, 3186, 3124, 2954, 2764, 2584, 2412, 2255,     6.000
     &   2112, 1987, 1871, 1716, 1465, 1183,  899,  642,  435,  293,     6.000
     &    178,   85,   10,  -84, -176, -227/                             6.000
!
      data (ktab(it, 31, 2), it = 1, nnt) /                              6.250
     &  -3782,-3993,-3716,-3435,-3148,-2852,-2545,-2228,-1910,-1600,     6.250
     &  -1310,-1045, -809, -600, -411, -233,  -65,   88,  224,  337,     6.250
     &    460,  536,  613,  756,  957, 1185, 1422, 1664, 1909, 2147,     6.250
     &   2534, 2881, 3155, 3316, 3321, 3197, 3025, 2853, 2684, 2525,     6.250
     &   2378, 2243, 2112, 1956, 1713, 1426, 1136,  866,  636,  464,     6.250
     &    330,  217,  124,   21,  -90, -166/                             6.250
!
      data (ktab(it, 32, 2), it = 1, nnt) /                              6.500
     &  -3516,-3904,-3627,-3348,-3066,-2777,-2479,-2172,-1860,-1553,     6.500
     &  -1258, -985, -737, -515, -315, -131,   42,  205,  350,  474,     6.500
     &    616,  709,  784,  906, 1089, 1305, 1532, 1764, 1997, 2223,     6.500
     &   2596, 2940, 3228, 3428, 3496, 3426, 3279, 3117, 2954, 2793,     6.500
     &   2643, 2499, 2355, 2195, 1963, 1672, 1378, 1101,  853,  655,     6.500
     &    504,  371,  259,  149,   19,  -84/                             6.500
!
      data (ktab(it, 33, 2), it = 1, nnt) /                              6.750
     &  -3193,-3813,-3537,-3260,-2981,-2698,-2409,-2111,-1808,-1505,     6.750
     &  -1211, -932, -675, -442, -232,  -40,  139,  309,  463,  599,     6.750
     &    759,  870,  951, 1059, 1223, 1425, 1642, 1863, 2085, 2302,     6.750
     &   2660, 2999, 3295, 3522, 3643, 3633, 3525, 3384, 3230, 3071,     6.750
     &   2915, 2761, 2605, 2436, 2210, 1919, 1620, 1335, 1076,  867,     6.750
     &    700,  552,  421,  297,  155,   24/                             6.750
!
      data (ktab(it, 34, 2), it = 1, nnt) /                              7.000
     &  -2988,-3722,-3447,-3171,-2894,-2616,-2334,-2046,-1752,-1456,     7.000
     &  -1165, -884, -621, -379, -159,   40,  226,  401,  564,  710,     7.000
     &    889, 1020, 1113, 1213, 1360, 1549, 1754, 1966, 2178, 2385,     7.000
     &   2729, 3062, 3361, 3608, 3771, 3820, 3758, 3640, 3502, 3349,     7.000
     &   3191, 3029, 2857, 2666, 2447, 2168, 1867, 1581, 1316, 1085,     7.000
     &    900,  746,  602,  466,  315,  159/                             7.000
!
      data (ktab(it, 35, 2), it = 1, nnt) /                              7.250
     &  -2878,-3632,-3357,-3082,-2807,-2533,-2257,-1977,-1692,-1404,     7.250
     &  -1118, -839, -573, -325,  -97,  110,  302,  482,  651,  806,     7.250
     &   1003, 1154, 1264, 1364, 1499, 1675, 1870, 2073, 2275, 2473,     7.250
     &   2804, 3129, 3428, 3687, 3880, 3978, 3969, 3885, 3770, 3629,     7.250
     &   3472, 3303, 3120, 2914, 2688, 2411, 2108, 1822, 1556, 1321,     7.250
     &   1124,  960,  799,  647,  494,  323/                             7.250
!
      data (ktab(it, 36, 2), it = 1, nnt) /                              7.500
     &  -2773,-3541,-3267,-2993,-2720,-2449,-2177,-1904,-1628,-1349,     7.500
     &  -1070, -795, -529, -278,  -44,  170,  368,  552,  726,  888,     7.500
     &   1102, 1273, 1403, 1511, 1638, 1803, 1989, 2183, 2377, 2567,     7.500
     &   2885, 3200, 3498, 3764, 3980, 4120, 4162, 4117, 4025, 3902,     7.500
     &   3751, 3578, 3383, 3163, 2929, 2658, 2355, 2066, 1800, 1561,     7.500
     &   1357, 1187, 1017,  853,  688,  503/                             7.500
!
      data (ktab(it, 37, 2), it = 1, nnt) /                              7.750
     &  -2672,-3450,-3177,-2904,-2634,-2365,-2097,-1830,-1561,-1291,     7.750
     &  -1020, -751, -488, -235,    2,  222,  425,  613,  790,  957,     7.750
     &   1183, 1373, 1523, 1645, 1772, 1930, 2109, 2296, 2482, 2665,     7.750
     &   2972, 3278, 3572, 3841, 4071, 4241, 4327, 4327, 4265, 4166,     7.750
     &   4027, 3853, 3650, 3416, 3169, 2898, 2596, 2305, 2040, 1801,     7.750
     &   1592, 1416, 1244, 1071,  897,  705/                             7.750
!
      data (ktab(it, 38, 2), it = 1, nnt) /                              8.000
     &  -2574,-3360,-3089,-2817,-2548,-2282,-2017,-1755,-1493,-1230,     8.000
     &   -968, -705, -447, -196,   42,  266,  473,  665,  845, 1015,     8.000
     &   1249, 1455, 1626, 1766, 1900, 2055, 2230, 2411, 2591, 2768,     8.000
     &   3064, 3361, 3650, 3920, 4160, 4353, 4478, 4520, 4492, 4415,     8.000
     &   4293, 4124, 3914, 3669, 3410, 3136, 2840, 2548, 2282, 2043,     8.000
     &   1831, 1650, 1476, 1297, 1116,  920/                             8.000
!
      data (ktab(it, 1, 3), it = 1, nnt) /                              -2.000
     &  -5171,-5106,-5020,-4865,-4673,-4489,-4342,-4241,-4182,-4158,    -2.000
     &  -4146,-4116,-4025,-3836,-3572,-3266,-2937,-2599,-2261,-1927,    -2.000
     &  -1446,-1006, -680, -543, -513, -514, -522, -527, -524, -510,    -2.000
     &   -489, -478, -479, -485, -481, -471, -447, -439, -434, -433,    -2.000
     &   -434, -430, -428, -428, -429, -428, -427, -430, -433, -434,    -2.000
     &   -432, -426, -418, -405, -386, -360/                            -2.000
!
      data (ktab(it, 2, 3), it = 1, nnt) /                              -1.500
     &  -5200,-5115,-5050,-4941,-4773,-4584,-4413,-4283,-4195,-4149,    -1.500
     &  -4132,-4117,-4070,-3943,-3724,-3446,-3137,-2813,-2484,-2156,    -1.500
     &  -1676,-1222, -829, -586, -508, -497, -505, -515, -520, -514,    -1.500
     &   -487, -474, -468, -477, -478, -469, -449, -436, -431, -428,    -1.500
     &   -431, -429, -426, -424, -426, -426, -425, -427, -432, -433,    -1.500
     &   -431, -426, -418, -405, -386, -360/                            -1.500
!
      data (ktab(it, 3, 3), it = 1, nnt) /                              -1.000
     &  -5242,-5131,-5057,-4981,-4852,-4677,-4494,-4337,-4220,-4145,    -1.000
     &  -4111,-4099,-4079,-4008,-3846,-3605,-3320,-3013,-2696,-2376,    -1.000
     &  -1901,-1440,-1013, -677, -510, -471, -475, -490, -503, -508,    -1.000
     &   -485, -469, -457, -465, -471, -464, -450, -432, -428, -423,    -1.000
     &   -426, -426, -423, -421, -423, -425, -423, -424, -429, -432,    -1.000
     &   -431, -426, -418, -405, -386, -360/                            -1.000
!
      data (ktab(it, 4, 3), it = 1, nnt) /                              -0.500
     &  -5285,-5154,-5057,-4986,-4894,-4747,-4569,-4392,-4247,-4144,    -0.500
     &  -4082,-4061,-4053,-4023,-3924,-3734,-3480,-3193,-2890,-2581,    -0.500
     &  -2113,-1651,-1206, -809, -536, -431, -421, -440, -464, -482,    -0.500
     &   -479, -461, -447, -448, -459, -458, -446, -427, -422, -418,    -0.500
     &   -420, -422, -421, -418, -419, -422, -421, -420, -425, -431,    -0.500
     &   -431, -426, -418, -405, -386, -360/                            -0.500
!
      data (ktab(it, 5, 3), it = 1, nnt) /                               0.000
     &  -5311,-5173,-5052,-4964,-4889,-4776,-4614,-4431,-4260,-4127,     0.000
     &  -4037,-3992,-3984,-3980,-3937,-3809,-3598,-3338,-3054,-2757,     0.000
     &  -2302,-1842,-1387, -954, -592, -385, -330, -346, -383, -420,     0.000
     &   -454, -443, -432, -425, -440, -445, -437, -422, -414, -411,     0.000
     &   -411, -417, -418, -415, -414, -418, -419, -416, -420, -428,     0.000
     &   -430, -426, -418, -405, -386, -360/                             0.000
!
      data (ktab(it, 6, 3), it = 1, nnt) /                               0.500
     &  -5313,-5165,-5027,-4916,-4833,-4743,-4606,-4429,-4241,-4076,     0.500
     &  -3954,-3879,-3853,-3859,-3857,-3798,-3644,-3420,-3159,-2880,     0.500
     &  -2443,-1992,-1533,-1079, -661, -346, -202, -191, -233, -292,     0.500
     &   -379, -400, -403, -397, -410, -422, -420, -410, -401, -401,     0.500
     &   -401, -408, -412, -411, -408, -412, -415, -413, -415, -424,     0.500
     &   -428, -425, -418, -405, -386, -360/                             0.500
!
      data (ktab(it, 7, 3), it = 1, nnt) /                               1.000
     &  -5305,-5120,-4963,-4827,-4724,-4637,-4525,-4366,-4175,-3986,     1.000
     &  -3827,-3715,-3654,-3646,-3664,-3659,-3577,-3403,-3171,-2912,     1.000
     &  -2500,-2067,-1617,-1157, -709, -316,  -57,   27,    0,  -72,     1.000
     &   -217, -303, -336, -349, -363, -386, -391, -384, -378, -384,     1.000
     &   -386, -395, -403, -406, -403, -405, -410, -409, -408, -417,     1.000
     &   -425, -424, -417, -405, -386, -360/                             1.000
!
      data (ktab(it, 8, 3), it = 1, nnt) /                               1.250
     &  -5300,-5088,-4915,-4765,-4647,-4555,-4451,-4303,-4118,-3922,     1.250
     &  -3748,-3616,-3536,-3511,-3526,-3538,-3491,-3348,-3135,-2887,     1.250
     &  -2489,-2068,-1628,-1172, -718, -300,   11,  150,  145,   73,     1.250
     &    -94, -221, -278, -310, -330, -359, -369, -366, -362, -370,     1.250
     &   -376, -385, -397, -402, -401, -400, -407, -407, -405, -413,     1.250
     &   -422, -423, -417, -405, -386, -360/                             1.250
!
      data (ktab(it, 9, 3), it = 1, nnt) /                               1.500
     &  -5294,-5057,-4859,-4691,-4555,-4454,-4355,-4219,-4042,-3844,     1.500
     &  -3658,-3507,-3406,-3360,-3365,-3386,-3369,-3261,-3070,-2836,     1.500
     &  -2451,-2043,-1615,-1168, -713, -280,   73,  273,  305,  242,     1.500
     &     57, -112, -200, -256, -285, -325, -344, -342, -339, -350,     1.500
     &   -362, -374, -389, -396, -397, -395, -403, -405, -402, -408,     1.500
     &   -420, -422, -416, -405, -386, -360/                             1.500
!
      data (ktab(it, 10, 3), it = 1, nnt) /                              1.750
     &  -5284,-5027,-4799,-4610,-4455,-4337,-4235,-4111,-3944,-3750,     1.750
     &  -3556,-3390,-3270,-3203,-3191,-3210,-3214,-3143,-2979,-2760,     1.750
     &  -2389,-1994,-1580,-1145, -696, -256,  128,  386,  470,  426,     1.750
     &    236,   28,  -97, -180, -228, -278, -309, -312, -310, -323,     1.750
     &   -343, -359, -378, -389, -393, -391, -397, -402, -399, -403,     1.750
     &   -416, -420, -416, -404, -386, -360/                             1.750
!
      data (ktab(it, 11, 3), it = 1, nnt) /                              2.000
     &  -5271,-4999,-4742,-4525,-4345,-4207,-4099,-3982,-3828,-3641,     2.000
     &  -3444,-3266,-3127,-3039,-3006,-3017,-3034,-2996,-2864,-2663,     2.000
     &  -2308,-1926,-1526,-1105, -666, -225,  179,  485,  629,  619,     2.000
     &    436,  197,   28,  -82, -153, -216, -265, -276, -273, -286,     2.000
     &   -317, -339, -364, -380, -387, -385, -391, -398, -395, -397,     2.000
     &   -411, -418, -415, -404, -385, -360/                             2.000
!
      data (ktab(it, 12, 3), it = 1, nnt) /                              2.250
     &  -5253,-4972,-4689,-4441,-4237,-4075,-3950,-3834,-3692,-3517,     2.250
     &  -3322,-3136,-2982,-2875,-2822,-2817,-2834,-2822,-2725,-2548,     2.250
     &  -2211,-1842,-1455,-1049, -623, -188,  226,  568,  771,  809,     2.250
     &    649,  393,  181,   41,  -60, -138, -206, -230, -230, -240,     2.250
     &   -279, -311, -343, -367, -379, -379, -383, -392, -391, -392,     2.250
     &   -405, -415, -413, -404, -385, -360/                             2.250
!
      data (ktab(it, 13, 3), it = 1, nnt) /                              2.500
     &  -5229,-4944,-4642,-4363,-4130,-3941,-3797,-3675,-3542,-3380,     2.500
     &  -3192,-3002,-2835,-2710,-2636,-2612,-2624,-2629,-2567,-2418,     2.500
     &  -2103,-1746,-1372, -981, -570, -144,  273,  640,  894,  989,     2.500
     &    869,  608,  358,  186,   54,  -40, -129, -173, -178, -184,     2.500
     &   -228, -274, -316, -350, -368, -372, -374, -385, -387, -386,     2.500
     &   -398, -411, -412, -403, -385, -360/                             2.500
!
      data (ktab(it, 14, 3), it = 1, nnt) /                              2.750
     &  -5199,-4914,-4599,-4293,-4031,-3816,-3647,-3512,-3381,-3230,     2.750
     &  -3053,-2864,-2689,-2549,-2456,-2412,-2411,-2421,-2389,-2272,     2.750
     &  -1984,-1641,-1279, -902, -506,  -92,  321,  701,  994, 1146,     2.750
     &   1087,  837,  559,  357,  194,   75,  -32,  -97, -114, -120,     2.750
     &   -160, -222, -277, -324, -351, -362, -364, -375, -381, -379,     2.750
     &   -390, -406, -409, -402, -385, -359/                             2.750
!
      data (ktab(it, 15, 3), it = 1, nnt) /                              3.000
     &  -5162,-4881,-4558,-4232,-3941,-3699,-3502,-3348,-3214,-3072,     3.000
     &  -2907,-2723,-2543,-2391,-2279,-2216,-2198,-2206,-2195,-2113,     3.000
     &  -1857,-1530,-1180, -815, -434,  -34,  371,  757, 1078, 1282,     3.000
     &   1299, 1072,  778,  544,  356,  212,   86,   -5,  -38,  -44,     3.000
     &    -78, -153, -223, -288, -329, -348, -351, -364, -373, -372,     3.000
     &   -381, -400, -406, -400, -384, -359/                             3.000
!
      data (ktab(it, 16, 3), it = 1, nnt) /                              3.200
     &  -5126,-4850,-4525,-4189,-3877,-3613,-3396,-3223,-3079,-2940,     3.200
     &  -2785,-2607,-2427,-2267,-2143,-2065,-2032,-2032,-2031,-1973,     3.200
     &  -1748,-1436,-1095, -741, -371,   16,  414,  799, 1134, 1371,     3.200
     &   1457, 1257,  964,  709,  506,  338,  198,   87,   36,   23,     3.200
     &     -2,  -83, -168, -247, -303, -332, -340, -351, -365, -366,     3.200
     &   -373, -393, -402, -399, -383, -359/                             3.200
!
      data (ktab(it, 17, 3), it = 1, nnt) /                              3.400
     &  -5084,-4816,-4492,-4148,-3821,-3535,-3297,-3103,-2947,-2806,     3.400
     &  -2659,-2489,-2310,-2145,-2010,-1917,-1869,-1859,-1861,-1826,     3.400
     &  -1634,-1339,-1008, -663, -304,   71,  459,  841, 1185, 1448,     3.400
     &   1603, 1442, 1156,  881,  665,  478,  324,  194,  124,  101,     3.400
     &     81,    0,  -99, -196, -269, -311, -325, -337, -354, -358,     3.400
     &   -363, -385, -398, -397, -383, -359/                             3.400
!
      data (ktab(it, 18, 3), it = 1, nnt) /                              3.600
     &  -5037,-4776,-4456,-4110,-3770,-3465,-3206,-2992,-2819,-2673,     3.600
     &  -2529,-2368,-2193,-2024,-1880,-1774,-1712,-1689,-1689,-1669,     3.600
     &  -1514,-1238, -918, -582, -233,  130,  507,  883, 1231, 1513,     3.600
     &   1734, 1621, 1350, 1063,  835,  633,  465,  318,  227,  192,     3.600
     &    173,   96,  -16, -130, -223, -282, -307, -320, -340, -349,     3.600
     &   -353, -375, -393, -394, -381, -358/                             3.600
!
      data (ktab(it, 19, 3), it = 1, nnt) /                              3.800
     &  -4985,-4731,-4416,-4070,-3724,-3402,-3124,-2890,-2698,-2541,     3.800
     &  -2398,-2245,-2074,-1904,-1754,-1636,-1559,-1523,-1516,-1506,     3.800
     &  -1386,-1133, -825, -498, -159,  192,  558,  926, 1274, 1569,     3.800
     &   1847, 1792, 1544, 1252, 1013,  800,  619,  459,  348,  296,     3.800
     &    273,  206,   82,  -49, -164, -243, -282, -299, -322, -336,     3.800
     &   -342, -363, -385, -390, -379, -357/                             3.800
!
      data (ktab(it, 20, 3), it = 1, nnt) /                              4.000
     &  -4926,-4680,-4372,-4031,-3681,-3347,-3050,-2796,-2585,-2414,     4.000
     &  -2267,-2119,-1955,-1786,-1630,-1502,-1411,-1362,-1345,-1339,     4.000
     &  -1252,-1024, -730, -412,  -83,  258,  612,  972, 1317, 1620,     4.000
     &   1945, 1953, 1734, 1445, 1195,  976,  784,  614,  483,  413,     4.000
     &    384,  327,  196,   48,  -91, -193, -251, -274, -301, -321,     4.000
     &   -328, -349, -377, -385, -377, -356/                             4.000
!
      data (ktab(it, 21, 3), it = 1, nnt) /                              4.200
     &  -4869,-4631,-4329,-3993,-3643,-3300,-2988,-2716,-2485,-2297,     4.200
     &  -2141,-1994,-1835,-1667,-1508,-1372,-1269,-1206,-1178,-1169,     4.200
     &  -1108, -910, -631, -323,   -3,  327,  670, 1019, 1359, 1665,     4.200
     &   2027, 2101, 1920, 1643, 1384, 1163,  962,  784,  639,  550,     4.200
     &    508,  459,  326,  163,    0, -126, -207, -243, -273, -300,     4.200
     &   -312, -332, -364, -379, -373, -354/                             4.200
!
      data (ktab(it, 22, 3), it = 1, nnt) /                              4.400
     &  -4805,-4578,-4281,-3952,-3605,-3258,-2932,-2643,-2395,-2187,     4.400
     &  -2019,-1869,-1714,-1549,-1388,-1245,-1132,-1055,-1014,-1000,     4.400
     &   -958, -790, -529, -232,   77,  398,  730, 1069, 1402, 1708,     4.400
     &   2097, 2233, 2101, 1841, 1575, 1353, 1146,  964,  806,  700,     4.400
     &    645,  599,  470,  294,  110,  -45, -151, -204, -239, -274,     4.400
     &   -292, -312, -349, -370, -368, -351/                             4.400
!
      data (ktab(it, 23, 3), it = 1, nnt) /                              4.600
     &  -4737,-4522,-4230,-3909,-3567,-3219,-2884,-2580,-2315,-2089,     4.600
     &  -1904,-1748,-1594,-1432,-1270,-1122, -999, -909, -856, -833,     4.600
     &   -801, -664, -423, -138,  162,  473,  793, 1121, 1447, 1750,     4.600
     &   2157, 2348, 2272, 2037, 1772, 1547, 1338, 1151,  986,  866,     4.600
     &    795,  747,  627,  440,  237,   54,  -79, -154, -198, -241,     4.600
     &   -267, -289, -330, -359, -362, -347/                             4.600
!
      data (ktab(it, 24, 3), it = 1, nnt) /                              4.800
     &  -4663,-4465,-4177,-3864,-3528,-3182,-2842,-2525,-2244,-2002,     4.800
     &  -1800,-1631,-1476,-1315,-1153,-1001, -870, -769, -703, -669,     4.800
     &   -639, -531, -313,  -41,  249,  550,  860, 1176, 1493, 1792,     4.800
     &   2210, 2447, 2432, 2230, 1972, 1744, 1536, 1346, 1177, 1044,     4.800
     &    959,  903,  795,  601,  382,  173,   10,  -91, -148, -199,     4.800
     &   -234, -261, -305, -343, -353, -342/                             4.800
!
      data (ktab(it, 25, 3), it = 1, nnt) /                              5.000
     &  -4581,-4404,-4119,-3814,-3486,-3144,-2802,-2477,-2182,-1923,     5.000
     &  -1704,-1521,-1360,-1200,-1037, -882, -745, -634, -555, -509,     5.000
     &   -476, -393, -199,   58,  339,  630,  929, 1234, 1542, 1834,     5.000
     &   2258, 2533, 2579, 2417, 2172, 1941, 1735, 1543, 1371, 1228,     5.000
     &   1130, 1066,  968,  771,  540,  308,  115,  -15,  -88, -148,     5.000
     &   -195, -226, -274, -323, -341, -335/                             5.000
!
      data (ktab(it, 26, 3), it = 1, nnt) /                              5.200
     &  -4490,-4342,-4059,-3760,-3441,-3105,-2764,-2433,-2127,-1854,     5.200
     &  -1619,-1421,-1251,-1087, -924, -766, -624, -504, -414, -355,     5.200
     &   -312, -246,  -78,  162,  432,  713, 1001, 1296, 1593, 1878,     5.200
     &   2303, 2604, 2708, 2597, 2375, 2144, 1941, 1748, 1575, 1427,     5.200
     &   1316, 1238, 1146,  954,  712,  463,  243,   82,  -14,  -86,     5.200
     &   -144, -184, -236, -295, -323, -323/                             5.200
!
      data (ktab(it, 27, 3), it = 1, nnt) /                              5.400
     &  -4388,-4277,-3997,-3703,-3392,-3063,-2726,-2393,-2078,-1793,     5.400
     &  -1543,-1330,-1147, -979, -814, -653, -506, -379, -277, -207,     5.400
     &   -151,  -96,   47,  270,  528,  799, 1076, 1360, 1647, 1924,     5.400
     &   2346, 2667, 2822, 2767, 2572, 2346, 2145, 1955, 1782, 1628,     5.400
     &   1506, 1415, 1325, 1141,  893,  630,  386,  196,   72,  -12,     5.400
     &    -83, -133, -188, -258, -300, -309/                             5.400
!
      data (ktab(it, 28, 3), it = 1, nnt) /                              5.600
     &  -4274,-4211,-3932,-3642,-3338,-3019,-2687,-2354,-2034,-1739,     5.600
     &  -1476,-1248,-1052, -877, -708, -545, -393, -259, -148,  -66,     5.600
     &      4,   56,  177,  381,  627,  887, 1155, 1428, 1705, 1972,     5.600
     &   2389, 2723, 2919, 2922, 2764, 2551, 2351, 2164, 1990, 1834,     5.600
     &   1703, 1601, 1505, 1333, 1080,  808,  545,  329,  177,   75,     5.600
     &     -9,  -72, -132, -212, -269, -289/                             5.600
!
      data (ktab(it, 29, 3), it = 1, nnt) /                              5.800
     &  -4144,-4143,-3865,-3579,-3282,-2970,-2645,-2316,-1993,-1691,     5.800
     &  -1417,-1176, -967, -782, -609, -441, -285, -144,  -26,   66,     5.800
     &    154,  209,  311,  496,  729,  979, 1236, 1499, 1765, 2024,     5.800
     &   2433, 2775, 3003, 3062, 2950, 2756, 2559, 2375, 2202, 2045,     5.800
     &   1908, 1793, 1688, 1526, 1272,  994,  719,  480,  300,  177,     5.800
     &     78,    1,  -66, -154, -228, -262/                             5.800
!
      data (ktab(it, 30, 3), it = 1, nnt) /                              6.000
     &  -3995,-4074,-3796,-3513,-3222,-2918,-2601,-2277,-1955,-1648,     6.000
     &  -1365,-1112, -890, -694, -515, -343, -182,  -35,   91,  194,     6.000
     &    299,  362,  449,  615,  833, 1073, 1321, 1573, 1829, 2078,     6.000
     &   2478, 2824, 3077, 3187, 3125, 2955, 2764, 2585, 2414, 2256,     6.000
     &   2113, 1988, 1873, 1718, 1468, 1186,  901,  643,  437,  294,     6.000
     &    180,   87,   11,  -83, -176, -226/                             6.000
!
      data (ktab(it, 31, 3), it = 1, nnt) /                              6.250
     &  -3776,-3986,-3709,-3429,-3143,-2848,-2541,-2225,-1907,-1598,     6.250
     &  -1307,-1043, -807, -597, -408, -231,  -63,   91,  227,  341,     6.250
     &    464,  540,  618,  761,  962, 1190, 1427, 1669, 1912, 2150,     6.250
     &   2536, 2882, 3156, 3317, 3322, 3197, 3026, 2853, 2685, 2526,     6.250
     &   2379, 2243, 2113, 1957, 1715, 1428, 1138,  868,  638,  466,     6.250
     &    331,  218,  125,   22,  -89, -165/                             6.250
!
      data (ktab(it, 32, 3), it = 1, nnt) /                              6.500
     &  -3511,-3896,-3621,-3342,-3061,-2773,-2475,-2169,-1858,-1550,     6.500
     &  -1256, -983, -735, -513, -313, -129,   45,  207,  353,  477,     6.500
     &    620,  713,  789,  911, 1094, 1310, 1537, 1768, 2000, 2226,     6.500
     &   2598, 2941, 3229, 3429, 3497, 3427, 3280, 3118, 2955, 2794,     6.500
     &   2643, 2500, 2356, 2196, 1964, 1673, 1379, 1102,  855,  657,     6.500
     &    505,  372,  260,  150,   20,  -83/                             6.500
!
      data (ktab(it, 33, 3), it = 1, nnt) /                              6.750
     &  -3189,-3806,-3531,-3255,-2976,-2694,-2405,-2108,-1806,-1503,     6.750
     &  -1208, -930, -673, -440, -230,  -38,  142,  311,  466,  602,     6.750
     &    763,  874,  956, 1064, 1228, 1430, 1646, 1867, 2089, 2305,     6.750
     &   2662, 3001, 3296, 3523, 3644, 3634, 3525, 3384, 3230, 3072,     6.750
     &   2916, 2761, 2605, 2437, 2211, 1921, 1621, 1336, 1077,  868,     6.750
     &    701,  553,  422,  298,  156,   25/                             6.750
!
      data (ktab(it, 34, 3), it = 1, nnt) /                              7.000
     &  -2985,-3716,-3441,-3166,-2890,-2613,-2331,-2043,-1750,-1454,     7.000
     &  -1162, -882, -619, -377, -157,   42,  228,  404,  567,  713,     7.000
     &    893, 1024, 1118, 1218, 1365, 1553, 1759, 1970, 2182, 2388,     7.000
     &   2731, 3063, 3362, 3608, 3772, 3820, 3758, 3640, 3502, 3349,     7.000
     &   3191, 3029, 2857, 2667, 2448, 2169, 1869, 1583, 1317, 1086,     7.000
     &    901,  747,  603,  467,  316,  160/                             7.000
!
      data (ktab(it, 35, 3), it = 1, nnt) /                              7.250
     &  -2876,-3625,-3351,-3077,-2804,-2530,-2254,-1974,-1690,-1402,     7.250
     &  -1116, -837, -571, -323,  -95,  112,  304,  484,  654,  809,     7.250
     &   1007, 1158, 1269, 1369, 1504, 1679, 1875, 2077, 2279, 2476,     7.250
     &   2806, 3130, 3429, 3687, 3880, 3979, 3969, 3885, 3770, 3629,     7.250
     &   3472, 3303, 3120, 2914, 2688, 2412, 2110, 1823, 1557, 1322,     7.250
     &   1125,  960,  800,  648,  495,  323/                             7.250
!
      data (ktab(it, 36, 3), it = 1, nnt) /                              7.500
     &  -2771,-3535,-3262,-2989,-2717,-2446,-2175,-1902,-1626,-1347,     7.500
     &  -1068, -793, -527, -276,  -41,  173,  370,  555,  729,  892,     7.500
     &   1106, 1277, 1407, 1516, 1643, 1807, 1993, 2187, 2381, 2570,     7.500
     &   2887, 3202, 3499, 3764, 3980, 4120, 4162, 4117, 4025, 3902,     7.500
     &   3751, 3578, 3383, 3164, 2929, 2659, 2356, 2068, 1802, 1562,     7.500
     &   1357, 1187, 1018,  854,  689,  503/                             7.500
!
      data (ktab(it, 37, 3), it = 1, nnt) /                              7.750
     &  -2670,-3445,-3172,-2900,-2630,-2362,-2095,-1828,-1559,-1289,     7.750
     &  -1018, -749, -486, -233,    4,  224,  427,  616,  793,  960,     7.750
     &   1186, 1377, 1528, 1650, 1777, 1935, 2113, 2300, 2486, 2668,     7.750
     &   2973, 3279, 3573, 3841, 4071, 4241, 4327, 4327, 4265, 4167,     7.750
     &   4027, 3854, 3650, 3417, 3169, 2898, 2597, 2306, 2041, 1802,     7.750
     &   1592, 1417, 1244, 1072,  897,  705/                             7.750
!
      data (ktab(it, 38, 3), it = 1, nnt) /                              8.000
     &  -2572,-3355,-3084,-2813,-2545,-2279,-2015,-1753,-1491,-1228,     8.000
     &   -965, -703, -445, -193,   45,  268,  476,  668,  848, 1018,     8.000
     &   1253, 1459, 1630, 1770, 1905, 2060, 2234, 2415, 2594, 2771,     8.000
     &   3066, 3362, 3651, 3921, 4160, 4353, 4478, 4520, 4492, 4415,     8.000
     &   4293, 4124, 3914, 3669, 3410, 3137, 2841, 2549, 2283, 2044,     8.000
     &   1831, 1651, 1477, 1297, 1116,  920/                             8.000
!
      data (ktab(it, 1, 4), it = 1, nnt) /                              -2.000
     &  -5141,-5078,-4993,-4840,-4650,-4468,-4323,-4223,-4165,-4145,    -2.000
     &  -4137,-4110,-4021,-3833,-3568,-3261,-2932,-2594,-2255,-1921,    -2.000
     &  -1439,-1000, -673, -535, -505, -506, -515, -522, -520, -507,    -2.000
     &   -484, -471, -473, -482, -478, -468, -444, -437, -432, -431,    -2.000
     &   -432, -429, -426, -426, -427, -426, -425, -428, -432, -433,    -2.000
     &   -431, -426, -417, -404, -385, -359/                            -2.000
!
      data (ktab(it, 2, 4), it = 1, nnt) /                              -1.500
     &  -5168,-5084,-5021,-4913,-4747,-4561,-4391,-4262,-4175,-4131,    -1.500
     &  -4118,-4108,-4064,-3938,-3719,-3441,-3131,-2807,-2477,-2149,    -1.500
     &  -1669,-1215, -821, -578, -498, -487, -496, -507, -514, -509,    -1.500
     &   -482, -467, -461, -473, -475, -465, -446, -434, -429, -426,    -1.500
     &   -429, -427, -424, -422, -425, -425, -423, -425, -431, -432,    -1.500
     &   -431, -426, -417, -405, -385, -359/                            -1.500
!
      data (ktab(it, 3, 4), it = 1, nnt) /                              -1.000
     &  -5209,-5097,-5024,-4951,-4823,-4650,-4469,-4313,-4197,-4123,    -1.000
     &  -4092,-4086,-4070,-4002,-3841,-3600,-3314,-3006,-2688,-2368,    -1.000
     &  -1893,-1433,-1006, -668, -501, -460, -463, -480, -494, -501,    -1.000
     &   -479, -461, -448, -459, -467, -460, -446, -429, -425, -421,    -1.000
     &   -424, -424, -422, -419, -421, -423, -421, -422, -428, -431,    -1.000
     &   -431, -426, -417, -404, -385, -359/                            -1.000
!
      data (ktab(it, 4, 4), it = 1, nnt) /                              -0.500
     &  -5252,-5119,-5021,-4952,-4862,-4716,-4540,-4365,-4220,-4118,    -0.500
     &  -4058,-4041,-4039,-4014,-3917,-3727,-3473,-3185,-2882,-2572,    -0.500
     &  -2104,-1643,-1198, -801, -527, -420, -408, -427, -453, -472,    -0.500
     &   -471, -452, -435, -438, -454, -453, -441, -423, -419, -415,    -0.500
     &   -417, -420, -419, -416, -417, -420, -419, -418, -424, -430,    -0.500
     &   -430, -426, -417, -405, -385, -359/                            -0.500
!
      data (ktab(it, 5, 4), it = 1, nnt) /                               0.000
     &  -5278,-5137,-5014,-4926,-4852,-4740,-4580,-4398,-4228,-4096,     0.000
     &  -4007,-3965,-3963,-3965,-3927,-3801,-3590,-3329,-3044,-2746,     0.000
     &  -2291,-1832,-1378, -945, -584, -375, -318, -331, -368, -407,     0.000
     &   -443, -433, -420, -413, -432, -439, -431, -416, -410, -407,     0.000
     &   -408, -414, -415, -413, -411, -416, -417, -414, -419, -427,     0.000
     &   -429, -425, -417, -404, -385, -359/                             0.000
!
      data (ktab(it, 6, 4), it = 1, nnt) /                               0.500
     &  -5280,-5129,-4988,-4876,-4793,-4704,-4567,-4391,-4205,-4041,     0.500
     &  -3919,-3845,-3824,-3838,-3842,-3785,-3633,-3409,-3147,-2867,     0.500
     &  -2429,-1979,-1522,-1070, -652, -337, -191, -177, -217, -276,     0.500
     &   -365, -388, -389, -381, -398, -415, -413, -402, -395, -397,     0.500
     &   -397, -405, -409, -409, -406, -410, -413, -411, -413, -422,     0.500
     &   -427, -425, -417, -404, -385, -359/                             0.500
!
      data (ktab(it, 7, 4), it = 1, nnt) /                               1.000
     &  -5272,-5083,-4924,-4785,-4681,-4595,-4483,-4323,-4134,-3945,     1.000
     &  -3787,-3676,-3619,-3617,-3642,-3642,-3563,-3389,-3157,-2897,     1.000
     &  -2484,-2052,-1603,-1146, -699, -307,  -47,   39,   14,  -55,     1.000
     &   -199, -287, -320, -331, -347, -376, -382, -375, -371, -378,     1.000
     &   -381, -390, -400, -403, -401, -402, -408, -407, -406, -415,     1.000
     &   -424, -423, -417, -404, -385, -359/                             1.000
!
      data (ktab(it, 8, 4), it = 1, nnt) /                               1.250
     &  -5267,-5053,-4876,-4724,-4605,-4512,-4407,-4259,-4074,-3880,     1.250
     &  -3707,-3577,-3499,-3479,-3500,-3518,-3475,-3333,-3119,-2871,     1.250
     &  -2472,-2051,-1613,-1159, -707, -291,   20,  161,  159,   90,     1.250
     &    -75, -203, -261, -291, -312, -347, -360, -356, -353, -363,     1.250
     &   -370, -381, -393, -399, -398, -398, -404, -405, -403, -411,     1.250
     &   -421, -422, -416, -404, -385, -359/                             1.250
!
      data (ktab(it, 9, 4), it = 1, nnt) /                               1.500
     &  -5261,-5022,-4820,-4650,-4513,-4411,-4311,-4175,-3998,-3802,     1.500
     &  -3616,-3467,-3367,-3325,-3337,-3364,-3351,-3245,-3054,-2819,     1.500
     &  -2432,-2025,-1599,-1154, -702, -270,   82,  284,  319,  259,     1.500
     &     77,  -92, -181, -235, -264, -310, -333, -331, -329, -343,     1.500
     &   -356, -368, -385, -393, -394, -392, -400, -402, -399, -406,     1.500
     &   -418, -421, -416, -404, -385, -359/                             1.500
!
      data (ktab(it, 10, 4), it = 1, nnt) /                              1.750
     &  -5253,-4993,-4762,-4571,-4413,-4295,-4192,-4067,-3901,-3708,     1.750
     &  -3515,-3350,-3232,-3168,-3161,-3187,-3195,-3125,-2962,-2742,     1.750
     &  -2370,-1975,-1563,-1130, -684, -246,  137,  396,  483,  442,     1.750
     &    256,   49,  -77, -158, -205, -260, -297, -301, -299, -314,     1.750
     &   -336, -353, -373, -386, -389, -387, -394, -399, -396, -401,     1.750
     &   -414, -419, -415, -404, -385, -359/                             1.750
!
      data (ktab(it, 11, 4), it = 1, nnt) /                              2.000
     &  -5241,-4966,-4706,-4486,-4305,-4167,-4058,-3940,-3786,-3601,     2.000
     &  -3405,-3227,-3090,-3004,-2975,-2992,-3013,-2977,-2846,-2644,     2.000
     &  -2288,-1906,-1508,-1089, -653, -214,  188,  495,  641,  635,     2.000
     &    457,  219,   49,  -59, -128, -195, -251, -264, -261, -276,     2.000
     &   -309, -332, -358, -376, -384, -381, -388, -395, -392, -395,     2.000
     &   -410, -417, -414, -404, -385, -359/                             2.000
!
      data (ktab(it, 12, 4), it = 1, nnt) /                              2.250
     &  -5224,-4940,-4654,-4404,-4199,-4036,-3912,-3795,-3653,-3478,     2.250
     &  -3285,-3100,-2947,-2842,-2791,-2791,-2812,-2802,-2706,-2529,     2.250
     &  -2191,-1822,-1437,-1033, -610, -177,  236,  578,  783,  823,     2.250
     &    669,  415,  203,   64,  -33, -114, -189, -216, -216, -229,     2.250
     &   -269, -303, -337, -362, -375, -375, -380, -389, -388, -389,     2.250
     &   -403, -414, -413, -403, -385, -359/                             2.250
!
      data (ktab(it, 13, 4), it = 1, nnt) /                              2.500
     &  -5203,-4914,-4609,-4328,-4094,-3905,-3762,-3640,-3507,-3345,     2.500
     &  -3158,-2969,-2803,-2679,-2606,-2586,-2602,-2609,-2548,-2398,     2.500
     &  -2082,-1726,-1353, -964, -555, -132,  283,  649,  905, 1002,     2.500
     &    888,  630,  379,  210,   81,  -13, -110, -158, -162, -171,     2.500
     &   -217, -265, -309, -344, -364, -368, -370, -382, -384, -382,     2.500
     &   -396, -410, -411, -402, -385, -359/                             2.500
!
      data (ktab(it, 14, 4), it = 1, nnt) /                              2.750
     &  -5174,-4886,-4567,-4260,-3998,-3783,-3615,-3480,-3349,-3199,     2.750
     &  -3022,-2834,-2660,-2521,-2429,-2387,-2388,-2400,-2369,-2252,     2.750
     &  -1963,-1621,-1260, -885, -491,  -80,  331,  710, 1005, 1159,     2.750
     &   1104,  858,  581,  380,  220,  102,  -10,  -81,  -98, -105,     2.750
     &   -148, -211, -269, -317, -347, -358, -360, -372, -378, -376,     2.750
     &   -388, -404, -408, -401, -384, -359/                             2.750
!
      data (ktab(it, 15, 4), it = 1, nnt) /                              3.000
     &  -5138,-4854,-4529,-4202,-3911,-3669,-3474,-3321,-3186,-3044,     3.000
     &  -2880,-2697,-2518,-2365,-2254,-2191,-2175,-2185,-2175,-2092,     3.000
     &  -1836,-1509,-1160, -798, -419,  -22,  382,  766, 1088, 1293,     3.000
     &   1315, 1091,  798,  565,  381,  240,  109,   12,  -21,  -27,     3.000
     &    -65, -142, -214, -280, -323, -344, -347, -360, -370, -369,     3.000
     &   -378, -398, -405, -400, -384, -359/                             3.000
!
      data (ktab(it, 16, 4), it = 1, nnt) /                              3.200
     &  -5103,-4825,-4498,-4160,-3849,-3586,-3371,-3198,-3055,-2916,     3.200
     &  -2761,-2584,-2404,-2244,-2120,-2042,-2010,-2011,-2010,-1953,     3.200
     &  -1727,-1415,-1076, -723, -355,   29,  424,  809, 1144, 1382,     3.200
     &   1471, 1276,  984,  729,  528,  365,  222,  106,   54,   41,     3.200
     &     12,  -70, -157, -239, -297, -327, -335, -347, -361, -362,     3.200
     &   -369, -391, -401, -398, -383, -358/                             3.200
!
      data (ktab(it, 17, 4), it = 1, nnt) /                              3.400
     &  -5063,-4792,-4466,-4122,-3795,-3511,-3274,-3081,-2925,-2785,     3.400
     &  -2637,-2468,-2290,-2124,-1989,-1896,-1848,-1839,-1841,-1805,     3.400
     &  -1613,-1318, -988, -645, -288,   84,  470,  851, 1195, 1459,     3.400
     &   1616, 1458, 1174,  900,  686,  503,  348,  213,  141,  120,     3.400
     &     97,   12,  -87, -187, -262, -306, -321, -333, -350, -354,     3.400
     &   -360, -382, -397, -396, -382, -358/                             3.400
!
      data (ktab(it, 18, 4), it = 1, nnt) /                              3.600
     &  -5018,-4753,-4431,-4085,-3747,-3443,-3186,-2973,-2801,-2654,     3.600
     &  -2511,-2350,-2175,-2005,-1861,-1754,-1692,-1669,-1669,-1649,     3.600
     &  -1492,-1217, -898, -564, -217,  144,  519,  893, 1241, 1524,     3.600
     &   1746, 1636, 1366, 1080,  854,  655,  488,  337,  245,  211,     3.600
     &    190,  110,   -3, -120, -216, -276, -302, -315, -336, -345,     3.600
     &   -349, -372, -391, -393, -381, -357/                             3.600
!
      data (ktab(it, 19, 4), it = 1, nnt) /                              3.800
     &  -4967,-4710,-4393,-4048,-3703,-3383,-3107,-2874,-2683,-2526,     3.800
     &  -2383,-2229,-2058,-1888,-1737,-1618,-1541,-1504,-1497,-1486,     3.800
     &  -1365,-1112, -805, -480, -143,  206,  570,  937, 1284, 1579,     3.800
     &   1858, 1805, 1559, 1268, 1030,  820,  640,  478,  365,  314,     3.800
     &    291,  220,   95,  -38, -156, -237, -277, -294, -318, -332,     3.800
     &   -338, -360, -384, -389, -379, -356/                             3.800
!
      data (ktab(it, 20, 4), it = 1, nnt) /                              4.000
     &  -4910,-4662,-4352,-4011,-3662,-3331,-3035,-2783,-2572,-2401,     4.000
     &  -2254,-2106,-1942,-1772,-1615,-1485,-1394,-1343,-1327,-1320,     4.000
     &  -1231,-1003, -710, -394,  -66,  272,  625,  983, 1327, 1630,     4.000
     &   1956, 1965, 1748, 1459, 1209,  993,  804,  631,  498,  430,     4.000
     &    402,  341,  209,   60,  -82, -186, -245, -269, -296, -317,     4.000
     &   -324, -345, -375, -384, -376, -355/                             4.000
!
      data (ktab(it, 21, 4), it = 1, nnt) /                              4.200
     &  -4853,-4613,-4310,-3974,-3625,-3285,-2974,-2704,-2474,-2286,     4.200
     &  -2130,-1983,-1823,-1655,-1495,-1357,-1254,-1189,-1160,-1151,     4.200
     &  -1088, -889, -610, -304,   13,  342,  683, 1030, 1369, 1675,     4.200
     &   2037, 2111, 1932, 1655, 1397, 1178,  979,  800,  653,  566,     4.200
     &    525,  473,  339,  175,   11, -119, -201, -237, -268, -296,     4.200
     &   -308, -329, -362, -377, -372, -353/                             4.200
!
      data (ktab(it, 22, 4), it = 1, nnt) /                              4.400
     &  -4790,-4561,-4263,-3935,-3588,-3244,-2920,-2633,-2385,-2178,     4.400
     &  -2009,-1860,-1704,-1539,-1377,-1232,-1118,-1039, -998, -982,     4.400
     &   -938, -770, -509, -213,   95,  413,  743, 1080, 1412, 1718,     4.400
     &   2106, 2242, 2111, 1852, 1587, 1365, 1160,  978,  819,  714,     4.400
     &    660,  613,  483,  307,  121,  -36, -144, -198, -234, -270,     4.400
     &   -287, -308, -347, -369, -368, -351/                             4.400
!
      data (ktab(it, 23, 4), it = 1, nnt) /                              4.600
     &  -4722,-4505,-4213,-3892,-3552,-3206,-2873,-2571,-2306,-2081,     4.600
     &  -1896,-1739,-1586,-1423,-1260,-1110, -986, -895, -841, -816,     4.600
     &   -783, -645, -403, -120,  179,  488,  807, 1133, 1457, 1759,     4.600
     &   2166, 2357, 2281, 2047, 1782, 1557, 1350, 1164,  997,  877,     4.600
     &    809,  760,  639,  453,  249,   63,  -72, -148, -193, -236,     4.600
     &   -262, -284, -327, -357, -361, -347/                             4.600
!
      data (ktab(it, 24, 4), it = 1, nnt) /                              4.800
     &  -4647,-4448,-4159,-3847,-3513,-3169,-2831,-2516,-2236,-1994,     4.800
     &  -1792,-1624,-1468,-1307,-1144, -991, -859, -757, -689, -654,     4.800
     &   -622, -513, -294,  -23,  266,  566,  873, 1188, 1503, 1801,     4.800
     &   2218, 2455, 2440, 2238, 1980, 1753, 1546, 1357, 1187, 1054,     4.800
     &    970,  915,  806,  612,  393,  183,   17,  -85, -142, -194,     4.800
     &   -230, -256, -302, -341, -352, -341/                             4.800
!
      data (ktab(it, 25, 4), it = 1, nnt) /                              5.000
     &  -4565,-4388,-4103,-3798,-3472,-3132,-2792,-2468,-2174,-1916,     5.000
     &  -1697,-1515,-1354,-1193,-1030, -874, -735, -623, -543, -495,     5.000
     &   -461, -376, -181,   75,  355,  645,  942, 1246, 1552, 1843,     5.000
     &   2266, 2540, 2586, 2424, 2179, 1948, 1743, 1552, 1380, 1236,     5.000
     &   1139, 1076,  978,  782,  552,  318,  124,   -8,  -82, -143,     5.000
     &   -190, -221, -270, -321, -340, -334/                             5.000
!
      data (ktab(it, 26, 4), it = 1, nnt) /                              5.200
     &  -4475,-4326,-4043,-3745,-3427,-3093,-2754,-2425,-2120,-1848,     5.200
     &  -1613,-1416,-1245,-1081, -917, -759, -615, -495, -403, -343,     5.200
     &   -299, -231,  -61,  179,  448,  728, 1015, 1308, 1604, 1887,     5.200
     &   2310, 2611, 2714, 2603, 2381, 2150, 1947, 1756, 1583, 1434,     5.200
     &   1324, 1247, 1154,  963,  723,  473,  252,   89,   -8,  -80,     5.200
     &   -139, -179, -232, -292, -322, -323/                             5.200
!
      data (ktab(it, 27, 4), it = 1, nnt) /                              5.400
     &  -4374,-4261,-3981,-3688,-3379,-3052,-2716,-2385,-2071,-1787,     5.400
     &  -1538,-1325,-1142, -974, -808, -647, -499, -370, -268, -196,     5.400
     &   -139,  -82,   62,  286,  543,  813, 1090, 1372, 1658, 1933,     5.400
     &   2353, 2673, 2827, 2772, 2577, 2351, 2151, 1961, 1788, 1634,     5.400
     &   1512, 1422, 1332, 1150,  903,  640,  395,  203,   79,   -6,     5.400
     &    -78, -128, -184, -256, -299, -308/                             5.400
!
      data (ktab(it, 28, 4), it = 1, nnt) /                              5.600
     &  -4260,-4195,-3917,-3628,-3326,-3008,-2678,-2346,-2028,-1733,     5.600
     &  -1471,-1243,-1048, -872, -703, -539, -386, -251, -140,  -57,     5.600
     &     15,   69,  190,  396,  641,  901, 1168, 1440, 1715, 1981,     5.600
     &   2396, 2729, 2924, 2926, 2768, 2555, 2355, 2169, 1995, 1839,     5.600
     &   1708, 1606, 1511, 1340, 1089,  818,  554,  336,  183,   81,     5.600
     &     -5,  -67, -128, -209, -268, -288/                             5.600
!
      data (ktab(it, 29, 4), it = 1, nnt) /                              5.800
     &  -4131,-4128,-3851,-3566,-3270,-2960,-2637,-2308,-1987,-1686,     5.800
     &  -1413,-1172, -963, -777, -604, -436, -279, -138,  -18,   75,     5.800
     &    164,  220,  324,  510,  742,  992, 1249, 1510, 1775, 2032,     5.800
     &   2439, 2779, 3007, 3065, 2953, 2759, 2562, 2379, 2206, 2049,     5.800
     &   1912, 1798, 1693, 1532, 1280, 1003,  727,  487,  306,  183,     5.800
     &     83,    6,  -61, -151, -226, -261/                             5.800
!
      data (ktab(it, 30, 4), it = 1, nnt) /                              6.000
     &  -3983,-4059,-3783,-3501,-3210,-2908,-2593,-2270,-1949,-1642,     6.000
     &  -1360,-1108, -886, -690, -510, -339, -176,  -29,   98,  201,     6.000
     &    307,  372,  460,  627,  845, 1085, 1332, 1584, 1839, 2086,     6.000
     &   2484, 2828, 3080, 3189, 3127, 2957, 2767, 2588, 2416, 2259,     6.000
     &   2116, 1991, 1877, 1723, 1475, 1194,  909,  650,  443,  300,     6.000
     &    184,   91,   15,  -80, -173, -225/                             6.000
!
      data (ktab(it, 31, 4), it = 1, nnt) /                              6.250
     &  -3766,-3972,-3696,-3417,-3132,-2838,-2533,-2218,-1901,-1593,     6.250
     &  -1303,-1039, -803, -593, -404, -226,  -57,   97,  233,  348,     6.250
     &    473,  550,  629,  772,  974, 1201, 1438, 1679, 1922, 2158,     6.250
     &   2542, 2886, 3159, 3319, 3324, 3199, 3028, 2855, 2687, 2528,     6.250
     &   2381, 2246, 2116, 1961, 1721, 1435, 1145,  874,  644,  471,     6.250
     &    336,  222,  129,   26,  -87, -164/                             6.250
!
      data (ktab(it, 32, 4), it = 1, nnt) /                              6.500
     &  -3502,-3883,-3608,-3331,-3051,-2764,-2468,-2163,-1852,-1545,     6.500
     &  -1251, -979, -731, -509, -309, -124,   50,  213,  359,  484,     6.500
     &    628,  723,  800,  922, 1106, 1321, 1548, 1778, 2009, 2234,     6.500
     &   2603, 2945, 3231, 3431, 3499, 3428, 3281, 3119, 2956, 2795,     6.500
     &   2645, 2501, 2358, 2199, 1968, 1680, 1386, 1108,  860,  662,     6.500
     &    510,  375,  264,  154,   23,  -81/                             6.500
!
      data (ktab(it, 33, 4), it = 1, nnt) /                              6.750
     &  -3182,-3793,-3520,-3244,-2967,-2686,-2398,-2103,-1800,-1498,     6.750
     &  -1204, -926, -669, -436, -226,  -33,  147,  317,  472,  609,     6.750
     &    771,  883,  967, 1075, 1240, 1441, 1657, 1877, 2098, 2312,     6.750
     &   2668, 3004, 3298, 3524, 3645, 3635, 3526, 3385, 3232, 3073,     6.750
     &   2917, 2763, 2607, 2439, 2214, 1926, 1627, 1342, 1082,  873,     6.750
     &    705,  557,  425,  301,  159,   27/                             6.750
!
      data (ktab(it, 34, 4), it = 1, nnt) /                              7.000
     &  -2979,-3704,-3430,-3157,-2882,-2605,-2325,-2038,-1745,-1449,     7.000
     &  -1158, -878, -615, -373, -153,   46,  233,  409,  573,  720,     7.000
     &    901, 1033, 1129, 1229, 1377, 1565, 1770, 1980, 2191, 2396,     7.000
     &   2737, 3067, 3364, 3610, 3773, 3821, 3759, 3641, 3502, 3350,     7.000
     &   3192, 3030, 2858, 2669, 2450, 2173, 1874, 1587, 1321, 1090,     7.000
     &    905,  750,  606,  471,  319,  162/                             7.000
!
      data (ktab(it, 35, 4), it = 1, nnt) /                              7.250
     &  -2870,-3614,-3341,-3068,-2796,-2523,-2248,-1969,-1685,-1398,     7.250
     &  -1112, -833, -568, -319,  -91,  116,  309,  490,  660,  816,     7.250
     &   1015, 1167, 1279, 1380, 1516, 1691, 1886, 2087, 2288, 2484,     7.250
     &   2811, 3134, 3432, 3689, 3881, 3979, 3970, 3885, 3770, 3629,     7.250
     &   3473, 3304, 3121, 2916, 2690, 2415, 2114, 1827, 1561, 1325,     7.250
     &   1128,  963,  802,  651,  497,  325/                             7.250
!
      data (ktab(it, 36, 4), it = 1, nnt) /                              7.500
     &  -2766,-3524,-3252,-2981,-2710,-2440,-2169,-1897,-1621,-1343,     7.500
     &  -1064, -789, -524, -272,  -37,  177,  375,  560,  735,  898,     7.500
     &   1113, 1287, 1418, 1527, 1655, 1818, 2004, 2197, 2389, 2577,     7.500
     &   2892, 3205, 3501, 3766, 3981, 4121, 4162, 4117, 4025, 3903,     7.500
     &   3751, 3579, 3384, 3164, 2930, 2661, 2359, 2071, 1805, 1565,     7.500
     &   1360, 1190, 1020,  856,  691,  505/                             7.500
!
      data (ktab(it, 37, 4), it = 1, nnt) /                              7.750
     &  -2666,-3434,-3164,-2893,-2624,-2357,-2090,-1823,-1555,-1285,     7.750
     &  -1014, -745, -482, -229,    8,  229,  432,  621,  799,  967,     7.750
     &   1194, 1386, 1538, 1660, 1788, 1945, 2124, 2309, 2494, 2675,     7.750
     &   2978, 3282, 3575, 3843, 4072, 4241, 4328, 4327, 4265, 4167,     7.750
     &   4027, 3854, 3650, 3417, 3170, 2900, 2600, 2309, 2044, 1804,     7.750
     &   1595, 1419, 1246, 1074,  899,  707/                             7.750
!
      data (ktab(it, 38, 4), it = 1, nnt) /                              8.000
     &  -2568,-3345,-3076,-2807,-2539,-2274,-2010,-1748,-1486,-1224,     8.000
     &   -961, -699, -440, -189,   49,  273,  481,  674,  855, 1025,     8.000
     &   1261, 1468, 1640, 1781, 1916, 2071, 2244, 2425, 2602, 2777,     8.000
     &   3071, 3365, 3653, 3922, 4161, 4354, 4478, 4521, 4492, 4415,     8.000
     &   4293, 4124, 3914, 3670, 3411, 3138, 2843, 2551, 2286, 2046,     8.000
     &   1833, 1652, 1479, 1299, 1118,  922/                             8.000
!
      data (ktab(it, 1, 5), it = 1, nnt) /                              -2.000
     &  -5105,-5043,-4959,-4808,-4620,-4441,-4297,-4198,-4143,-4128,    -2.000
     &  -4125,-4102,-4015,-3828,-3563,-3255,-2926,-2586,-2247,-1912,    -2.000
     &  -1430, -990, -662, -522, -491, -493, -503, -513, -513, -500,    -2.000
     &   -477, -460, -464, -477, -474, -462, -440, -434, -429, -428,    -2.000
     &   -429, -425, -422, -421, -422, -422, -421, -425, -430, -431,    -2.000
     &   -429, -424, -416, -403, -384, -358/                            -2.000
!
      data (ktab(it, 2, 5), it = 1, nnt) /                              -1.500
     &  -5132,-5046,-4984,-4877,-4713,-4529,-4362,-4234,-4149,-4108,    -1.500
     &  -4102,-4097,-4056,-3932,-3713,-3434,-3124,-2798,-2468,-2139,    -1.500
     &  -1659,-1204, -810, -565, -483, -471, -481, -494, -504, -501,    -1.500
     &   -473, -455, -448, -466, -469, -459, -440, -429, -425, -422,    -1.500
     &   -425, -423, -419, -418, -420, -420, -418, -421, -427, -430,    -1.500
     &   -429, -424, -416, -403, -384, -358/                            -1.500
!
      data (ktab(it, 3, 5), it = 1, nnt) /                              -1.000
     &  -5172,-5058,-4984,-4912,-4785,-4614,-4435,-4280,-4166,-4095,    -1.000
     &  -4068,-4069,-4058,-3993,-3833,-3592,-3305,-2997,-2678,-2357,    -1.000
     &  -1881,-1421, -994, -655, -485, -442, -445, -463, -480, -489,    -1.000
     &   -468, -448, -432, -448, -460, -453, -438, -423, -420, -416,    -1.000
     &   -420, -420, -417, -414, -416, -418, -416, -417, -424, -429,    -1.000
     &   -429, -424, -416, -403, -384, -358/                            -1.000
!
      data (ktab(it, 4, 5), it = 1, nnt) /                              -0.500
     &  -5215,-5079,-4978,-4909,-4820,-4676,-4501,-4327,-4184,-4083,    -0.500
     &  -4027,-4016,-4021,-4001,-3906,-3717,-3462,-3174,-2869,-2558,    -0.500
     &  -2091,-1629,-1185, -788, -512, -403, -388, -406, -434, -456,    -0.500
     &   -458, -438, -417, -423, -444, -444, -431, -415, -413, -409,    -0.500
     &   -412, -415, -414, -410, -412, -414, -413, -413, -420, -427,    -0.500
     &   -428, -424, -416, -403, -384, -358/                            -0.500
!
      data (ktab(it, 5, 5), it = 1, nnt) /                               0.000
     &  -5240,-5096,-4970,-4880,-4806,-4695,-4536,-4355,-4186,-4054,     0.000
     &  -3968,-3930,-3936,-3946,-3912,-3788,-3577,-3315,-3029,-2731,     0.000
     &  -2275,-1817,-1364, -932, -571, -360, -299, -309, -345, -385,     0.000
     &   -425, -415, -399, -391, -418, -428, -419, -406, -402, -400,     0.000
     &   -401, -408, -410, -407, -405, -410, -411, -408, -414, -423,     0.000
     &   -427, -424, -416, -403, -384, -358/                             0.000
!
      data (ktab(it, 6, 5), it = 1, nnt) /                               0.500
     &  -5243,-5087,-4943,-4828,-4743,-4653,-4517,-4341,-4155,-3992,     0.500
     &  -3872,-3801,-3786,-3808,-3820,-3768,-3617,-3393,-3130,-2848,     0.500
     &  -2410,-1961,-1505,-1054, -638, -322, -174, -155, -192, -250,     0.500
     &   -341, -366, -366, -354, -378, -401, -399, -388, -385, -387,     0.500
     &   -388, -398, -403, -403, -399, -404, -407, -404, -407, -418,     0.500
     &   -424, -423, -415, -403, -384, -358/                             0.500
!
      data (ktab(it, 7, 5), it = 1, nnt) /                               1.000
     &  -5235,-5042,-4878,-4737,-4630,-4540,-4427,-4267,-4078,-3891,     1.000
     &  -3734,-3625,-3571,-3577,-3611,-3619,-3543,-3370,-3137,-2876,     1.000
     &  -2461,-2030,-1583,-1128, -683, -292,  -31,   58,   38,  -27,     1.000
     &   -170, -260, -293, -299, -318, -357, -366, -358, -356, -366,     1.000
     &   -370, -381, -392, -396, -393, -395, -401, -400, -399, -410,     1.000
     &   -421, -421, -415, -403, -384, -358/                             1.000
!
      data (ktab(it, 8, 5), it = 1, nnt) /                               1.250
     &  -5230,-5012,-4831,-4676,-4553,-4457,-4350,-4201,-4016,-3823,     1.250
     &  -3651,-3523,-3449,-3435,-3466,-3492,-3452,-3312,-3098,-2848,     1.250
     &  -2447,-2027,-1590,-1139, -690, -276,   35,  179,  182,  118,     1.250
     &    -44, -173, -231, -257, -279, -325, -343, -337, -336, -350,     1.250
     &   -358, -371, -385, -391, -390, -390, -397, -398, -395, -405,     1.250
     &   -418, -420, -415, -403, -384, -358/                             1.250
!
      data (ktab(it, 9, 5), it = 1, nnt) /                               1.500
     &  -5225,-4982,-4776,-4603,-4462,-4357,-4253,-4115,-3939,-3744,     1.500
     &  -3559,-3412,-3315,-3279,-3299,-3335,-3326,-3222,-3031,-2794,     1.500
     &  -2407,-2000,-1575,-1133, -684, -255,   97,  301,  341,  286,     1.500
     &    109,  -60, -149, -198, -226, -284, -314, -311, -309, -327,     1.500
     &   -342, -357, -375, -385, -386, -384, -392, -395, -391, -400,     1.500
     &   -414, -419, -414, -403, -384, -358/                             1.500
!
      data (ktab(it, 10, 5), it = 1, nnt) /                              1.750
     &  -5217,-4954,-4720,-4525,-4364,-4242,-4136,-4008,-3842,-3651,     1.750
     &  -3459,-3296,-3180,-3121,-3121,-3154,-3168,-3100,-2937,-2716,     1.750
     &  -2343,-1949,-1538,-1108, -665, -229,  152,  412,  503,  468,     1.750
     &    289,   83,  -42, -119, -163, -229, -275, -279, -277, -296,     1.750
     &   -321, -339, -363, -377, -381, -378, -386, -391, -388, -394,     1.750
     &   -410, -416, -413, -402, -384, -358/                             1.750
!
      data (ktab(it, 11, 5), it = 1, nnt) /                              2.000
     &  -5206,-4928,-4665,-4442,-4258,-4116,-4005,-3883,-3729,-3545,     2.000
     &  -3351,-3175,-3039,-2956,-2934,-2958,-2984,-2951,-2821,-2618,     2.000
     &  -2260,-1879,-1482,-1066, -633, -197,  203,  511,  660,  659,     2.000
     &    490,  254,   85,  -19,  -83, -158, -226, -241, -236, -256,     2.000
     &   -292, -316, -346, -366, -375, -372, -379, -387, -383, -387,     2.000
     &   -404, -414, -412, -402, -384, -358/                             2.000
!
      data (ktab(it, 12, 5), it = 1, nnt) /                              2.250
     &  -5191,-4903,-4615,-4362,-4154,-3989,-3862,-3742,-3599,-3426,     2.250
     &  -3234,-3051,-2899,-2797,-2750,-2756,-2782,-2775,-2680,-2502,     2.250
     &  -2163,-1794,-1410,-1009, -589, -159,  251,  593,  800,  846,     2.250
     &    701,  451,  240,  105,   13,  -71, -159, -191, -190, -206,     2.250
     &   -250, -286, -323, -351, -365, -365, -370, -380, -379, -380,     2.250
     &   -397, -410, -410, -401, -384, -358/                             2.250
!
      data (ktab(it, 13, 5), it = 1, nnt) /                              2.500
     &  -5170,-4879,-4571,-4287,-4051,-3861,-3716,-3592,-3457,-3296,     2.500
     &  -3111,-2923,-2759,-2636,-2567,-2551,-2571,-2581,-2521,-2371,     2.500
     &  -2053,-1697,-1326, -940, -534, -114,  299,  664,  922, 1023,     2.500
     &    919,  666,  417,  250,  129,   33,  -76, -130, -134, -145,     2.500
     &   -196, -245, -292, -331, -353, -357, -360, -373, -374, -373,     2.500
     &   -389, -406, -409, -401, -383, -358/                             2.500
!
      data (ktab(it, 14, 5), it = 1, nnt) /                              2.750
     &  -5143,-4852,-4531,-4221,-3958,-3743,-3573,-3437,-3305,-3155,     2.750
     &  -2980,-2793,-2620,-2482,-2392,-2353,-2357,-2372,-2341,-2224,     2.750
     &  -1934,-1592,-1233, -860, -470,  -62,  347,  725, 1021, 1178,     2.750
     &   1133,  892,  617,  418,  267,  152,   28,  -50,  -67,  -76,     2.750
     &   -124, -189, -250, -303, -335, -346, -349, -362, -368, -366,     2.750
     &   -380, -400, -406, -399, -383, -358/                             2.750
!
      data (ktab(it, 15, 5), it = 1, nnt) /                              3.000
     &  -5109,-4822,-4494,-4165,-3873,-3632,-3436,-3283,-3148,-3006,     3.000
     &  -2843,-2660,-2482,-2330,-2219,-2159,-2145,-2156,-2147,-2064,     3.000
     &  -1806,-1480,-1132, -773, -397,   -3,  398,  782, 1104, 1312,     3.000
     &   1341, 1124,  833,  602,  424,  290,  152,   45,   10,    5,     3.000
     &    -38, -117, -193, -264, -310, -332, -335, -349, -360, -358,     3.000
     &   -369, -392, -402, -398, -382, -357/                             3.000
!
      data (ktab(it, 16, 5), it = 1, nnt) /                              3.200
     &  -5076,-4794,-4464,-4125,-3814,-3552,-3337,-3164,-3021,-2881,     3.200
     &  -2727,-2551,-2371,-2211,-2088,-2011,-1981,-1983,-1982,-1924,     3.200
     &  -1697,-1385,-1048, -698, -333,   48,  441,  824, 1160, 1400,     3.200
     &   1495, 1306, 1017,  764,  569,  412,  265,  140,   87,   76,     3.200
     &     40,  -44, -134, -221, -282, -315, -323, -336, -351, -351,     3.200
     &   -360, -385, -398, -396, -381, -357/                             3.200
!
      data (ktab(it, 17, 5), it = 1, nnt) /                              3.400
     &  -5037,-4762,-4434,-4088,-3762,-3479,-3244,-3051,-2896,-2755,     3.400
     &  -2608,-2439,-2261,-2095,-1959,-1866,-1820,-1811,-1813,-1777,     3.400
     &  -1583,-1288, -960, -619, -265,  104,  488,  867, 1210, 1476,     3.400
     &   1638, 1487, 1205,  932,  723,  547,  392,  249,  175,  156,     3.400
     &    128,   40,  -62, -166, -245, -292, -307, -321, -340, -343,     3.400
     &   -349, -376, -393, -393, -380, -357/                             3.400
!
      data (ktab(it, 18, 5), it = 1, nnt) /                              3.600
     &  -4993,-4726,-4401,-4054,-3716,-3414,-3159,-2947,-2775,-2629,     3.600
     &  -2485,-2324,-2149,-1979,-1834,-1727,-1664,-1642,-1642,-1621,     3.600
     &  -1463,-1187, -869, -537, -194,  164,  537,  910, 1257, 1540,     3.600
     &   1766, 1662, 1395, 1110,  887,  696,  531,  373,  277,  247,     3.600
     &    222,  138,   23,  -97, -197, -262, -288, -302, -325, -333,     3.600
     &   -338, -365, -387, -390, -379, -356/                             3.600
!
      data (ktab(it, 19, 5), it = 1, nnt) /                              3.800
     &  -4943,-4684,-4365,-4019,-3674,-3357,-3082,-2851,-2660,-2504,     3.800
     &  -2361,-2207,-2036,-1865,-1712,-1593,-1515,-1478,-1470,-1459,     3.800
     &  -1336,-1082, -776, -453, -119,  228,  589,  954, 1300, 1596,     3.800
     &   1877, 1829, 1585, 1295, 1060,  856,  680,  513,  396,  350,     3.800
     &    325,  249,  123,  -13, -136, -221, -263, -280, -306, -320,     3.800
     &   -326, -352, -379, -386, -377, -355/                             3.800
!
      data (ktab(it, 20, 5), it = 1, nnt) /                              4.000
     &  -4888,-4638,-4326,-3984,-3636,-3307,-3014,-2763,-2553,-2383,     4.000
     &  -2236,-2087,-1922,-1751,-1593,-1462,-1369,-1318,-1300,-1293,     4.000
     &  -1202, -973, -680, -367,  -42,  294,  644, 1000, 1343, 1646,     4.000
     &   1973, 1987, 1772, 1483, 1235, 1024,  840,  665,  528,  464,     4.000
     &    437,  370,  238,   87,  -59, -168, -230, -254, -284, -305,     4.000
     &   -311, -336, -369, -381, -374, -354/                             4.000
!
      data (ktab(it, 21, 5), it = 1, nnt) /                              4.200
     &  -4831,-4590,-4284,-3948,-3600,-3263,-2954,-2686,-2458,-2270,     4.200
     &  -2114,-1966,-1806,-1637,-1475,-1336,-1231,-1166,-1136,-1125,     4.200
     &  -1059, -859, -581, -277,   38,  364,  702, 1048, 1386, 1691,     4.200
     &   2053, 2130, 1954, 1678, 1420, 1204, 1011,  831,  681,  597,     4.200
     &    559,  502,  368,  203,   35,  -99, -184, -222, -255, -283,     4.200
     &   -295, -318, -356, -374, -370, -352/                             4.200
!
      data (ktab(it, 22, 5), it = 1, nnt) /                              4.400
     &  -4768,-4537,-4238,-3910,-3565,-3223,-2902,-2617,-2371,-2165,     4.400
     &  -1996,-1846,-1689,-1523,-1359,-1214,-1097,-1017, -975, -957,     4.400
     &   -911, -741, -480, -186,  120,  436,  764, 1098, 1429, 1734,     4.400
     &   2122, 2259, 2130, 1872, 1607, 1388, 1188, 1006,  844,  741,     4.400
     &    691,  641,  510,  335,  147,  -15, -127, -182, -220, -257,     4.400
     &   -274, -297, -339, -364, -365, -349/                             4.400
!
      data (ktab(it, 23, 5), it = 1, nnt) /                              4.600
     &  -4700,-4482,-4189,-3868,-3529,-3186,-2856,-2556,-2293,-2069,     4.600
     &  -1885,-1727,-1573,-1409,-1244,-1094, -968, -875, -819, -793,     4.600
     &   -756, -617, -375,  -92,  205,  512,  828, 1152, 1474, 1776,     4.600
     &   2181, 2372, 2298, 2064, 1799, 1576, 1373, 1189, 1020,  901,     4.600
     &    837,  787,  665,  481,  276,   87,  -53, -131, -178, -223,     4.600
     &   -248, -272, -319, -352, -358, -345/                             4.600
!
      data (ktab(it, 24, 5), it = 1, nnt) /                              4.800
     &  -4626,-4425,-4136,-3824,-3491,-3150,-2814,-2502,-2224,-1983,     4.800
     &  -1782,-1613,-1457,-1295,-1130, -976, -842, -739, -669, -632,     4.800
     &   -598, -486, -266,    4,  292,  589,  895, 1207, 1521, 1817,     4.800
     &   2233, 2468, 2455, 2253, 1996, 1769, 1565, 1378, 1207, 1074,     4.800
     &    994,  939,  830,  639,  421,  208,   37,  -67, -126, -181,     4.800
     &   -216, -243, -292, -335, -348, -339/                             4.800
!
      data (ktab(it, 25, 5), it = 1, nnt) /                              5.000
     &  -4544,-4365,-4079,-3775,-3451,-3113,-2776,-2455,-2163,-1906,     5.000
     &  -1688,-1505,-1343,-1182,-1018, -861, -721, -607, -524, -475,     5.000
     &   -437, -350, -154,  102,  381,  669,  964, 1266, 1570, 1859,     5.000
     &   2280, 2552, 2598, 2437, 2193, 1962, 1758, 1570, 1398, 1253,     5.000
     &   1158, 1097,  998,  807,  579,  344,  145,    9,  -64, -129,     5.000
     &   -177, -207, -260, -314, -336, -332/                             5.000
!
      data (ktab(it, 26, 5), it = 1, nnt) /                              5.200
     &  -4453,-4303,-4020,-3723,-3407,-3075,-2739,-2412,-2109,-1839,     5.200
     &  -1604,-1407,-1236,-1072, -907, -747, -602, -480, -387, -325,     5.200
     &   -278, -207,  -35,  205,  473,  752, 1037, 1328, 1621, 1903,     5.200
     &   2324, 2622, 2725, 2614, 2392, 2161, 1960, 1771, 1598, 1448,     5.200
     &   1340, 1265, 1173,  986,  749,  499,  274,  108,    9,  -65,     5.200
     &   -126, -165, -220, -285, -318, -320/                             5.200
!
      data (ktab(it, 27, 5), it = 1, nnt) /                              5.400
     &  -4353,-4239,-3959,-3667,-3359,-3035,-2702,-2373,-2061,-1778,     5.400
     &  -1530,-1317,-1134, -965, -799, -637, -487, -357, -254, -180,     5.400
     &   -120,  -60,   86,  310,  568,  836, 1112, 1392, 1676, 1949,     5.400
     &   2366, 2684, 2836, 2781, 2586, 2361, 2161, 1973, 1801, 1646,     5.400
     &   1524, 1438, 1348, 1169,  927,  666,  418,  222,   97,    9,     5.400
     &    -65, -114, -171, -247, -294, -305/                             5.400
!
      data (ktab(it, 28, 5), it = 1, nnt) /                              5.600
     &  -4239,-4173,-3895,-3607,-3307,-2991,-2664,-2335,-2018,-1725,     5.600
     &  -1463,-1236,-1040, -864, -695, -530, -376, -240, -127,  -42,     5.600
     &     32,   89,  213,  419,  665,  924, 1189, 1459, 1733, 1997,     5.600
     &   2408, 2738, 2932, 2934, 2776, 2563, 2363, 2178, 2006, 1849,     5.600
     &   1718, 1619, 1525, 1356, 1111,  842,  577,  355,  202,   97,     5.600
     &      8,  -53, -114, -200, -262, -285/                             5.600
!
      data (ktab(it, 29, 5), it = 1, nnt) /                              5.800
     &  -4111,-4106,-3829,-3546,-3252,-2944,-2623,-2297,-1978,-1677,     5.800
     &  -1405,-1165, -956, -770, -596, -428, -269, -127,   -6,   88,     5.800
     &    180,  239,  344,  531,  764, 1014, 1270, 1530, 1793, 2048,     5.800
     &   2451, 2789, 3014, 3072, 2959, 2765, 2568, 2386, 2215, 2058,     5.800
     &   1920, 1808, 1705, 1546, 1300, 1026,  750,  506,  324,  200,     5.800
     &     97,   19,  -48, -141, -220, -257/                             5.800
!
      data (ktab(it, 30, 5), it = 1, nnt) /                              6.000
     &  -3965,-4038,-3762,-3482,-3193,-2893,-2580,-2259,-1940,-1634,     6.000
     &  -1353,-1101, -879, -683, -503, -331, -168,  -19,  109,  214,     6.000
     &    322,  389,  479,  647,  866, 1105, 1352, 1603, 1856, 2101,     6.000
     &   2495, 2836, 3087, 3194, 3132, 2962, 2771, 2593, 2423, 2266,     6.000
     &   2122, 1999, 1887, 1734, 1491, 1215,  930,  669,  460,  317,     6.000
     &    198,  104,   29,  -69, -166, -221/                             6.000
!
      data (ktab(it, 31, 5), it = 1, nnt) /                              6.250
     &  -3750,-3951,-3677,-3399,-3116,-2824,-2521,-2208,-1892,-1585,     6.250
     &  -1296,-1032, -797, -587, -398, -219,  -49,  106,  244,  360,     6.250
     &    487,  567,  648,  792,  994, 1222, 1458, 1697, 1938, 2172,     6.250
     &   2552, 2894, 3165, 3323, 3328, 3203, 3032, 2859, 2692, 2534,     6.250
     &   2386, 2252, 2123, 1970, 1735, 1453, 1165,  892,  660,  487,     6.250
     &    350,  235,  142,   37,  -79, -159/                             6.250
!
      data (ktab(it, 32, 5), it = 1, nnt) /                              6.500
     &  -3489,-3863,-3590,-3315,-3036,-2751,-2457,-2153,-1844,-1538,     6.500
     &  -1245, -972, -725, -503, -303, -117,   58,  222,  369,  496,     6.500
     &    642,  739,  818,  942, 1126, 1341, 1567, 1795, 2025, 2248,     6.500
     &   2613, 2952, 3237, 3434, 3502, 3431, 3283, 3122, 2959, 2799,     6.500
     &   2649, 2506, 2363, 2206, 1979, 1695, 1404, 1125,  875,  676,     6.500
     &    523,  387,  276,  165,   31,  -76/                             6.500
!
      data (ktab(it, 33, 5), it = 1, nnt) /                              6.750
     &  -3172,-3775,-3502,-3229,-2953,-2674,-2388,-2093,-1792,-1491,     6.750
     &  -1198, -920, -663, -430, -220,  -26,  154,  325,  482,  620,     6.750
     &    784,  899,  984, 1094, 1260, 1461, 1676, 1895, 2114, 2326,     6.750
     &   2678, 3012, 3304, 3528, 3648, 3637, 3529, 3387, 3234, 3076,     6.750
     &   2920, 2766, 2611, 2444, 2223, 1938, 1643, 1357, 1096,  886,     6.750
     &    718,  567,  436,  313,  168,   33/                             6.750
!
      data (ktab(it, 34, 5), it = 1, nnt) /                              7.000
     &  -2971,-3686,-3414,-3142,-2869,-2594,-2315,-2029,-1737,-1442,     7.000
     &  -1152, -872, -610, -368, -147,   53,  240,  417,  582,  731,     7.000
     &    914, 1049, 1146, 1248, 1397, 1585, 1789, 1998, 2207, 2410,     7.000
     &   2747, 3074, 3370, 3613, 3775, 3823, 3761, 3642, 3504, 3352,     7.000
     &   3194, 3032, 2861, 2673, 2457, 2183, 1887, 1601, 1334, 1101,     7.000
     &    916,  760,  615,  481,  328,  168/                             7.000
!
      data (ktab(it, 35, 5), it = 1, nnt) /                              7.250
     &  -2863,-3597,-3326,-3055,-2785,-2513,-2239,-1961,-1678,-1391,     7.250
     &  -1106, -827, -562, -313,  -85,  123,  316,  498,  669,  827,     7.250
     &   1028, 1183, 1296, 1399, 1535, 1711, 1905, 2105, 2304, 2498,     7.250
     &   2821, 3141, 3436, 3692, 3883, 3981, 3971, 3887, 3772, 3631,     7.250
     &   3474, 3306, 3123, 2919, 2695, 2424, 2125, 1839, 1573, 1335,     7.250
     &   1138,  972,  810,  660,  506,  332/                             7.250
!
      data (ktab(it, 36, 5), it = 1, nnt) /                              7.500
     &  -2760,-3508,-3239,-2969,-2700,-2431,-2161,-1890,-1615,-1336,     7.500
     &  -1058, -783, -518, -266,  -31,  184,  383,  569,  745,  909,     7.500
     &   1126, 1301, 1435, 1545, 1674, 1838, 2023, 2214, 2405, 2591,     7.500
     &   2901, 3212, 3506, 3769, 3983, 4122, 4163, 4118, 4026, 3904,     7.500
     &   3753, 3580, 3386, 3167, 2934, 2667, 2368, 2081, 1815, 1574,     7.500
     &   1368, 1197, 1026,  864,  699,  512/                             7.500
!
      data (ktab(it, 37, 5), it = 1, nnt) /                              7.750
     &  -2660,-3420,-3151,-2882,-2615,-2348,-2083,-1816,-1548,-1279,     7.750
     &  -1008, -739, -476, -223,   14,  236,  440,  630,  809,  978,     7.750
     &   1207, 1401, 1555, 1679, 1807, 1964, 2142, 2326, 2510, 2688,     7.750
     &   2987, 3289, 3579, 3846, 4073, 4242, 4328, 4328, 4266, 4168,     7.750
     &   4028, 3855, 3652, 3419, 3173, 2905, 2607, 2318, 2052, 1812,     7.750
     &   1602, 1425, 1252, 1080,  906,  713/                             7.750
!
      data (ktab(it, 38, 5), it = 1, nnt) /                              8.000
     &  -2563,-3332,-3064,-2797,-2531,-2267,-2004,-1742,-1480,-1218,     8.000
     &   -955, -693, -434, -183,   56,  280,  489,  683,  864, 1036,     8.000
     &   1274, 1483, 1657, 1799, 1934, 2089, 2262, 2441, 2617, 2790,     8.000
     &   3079, 3371, 3657, 3925, 4163, 4355, 4479, 4521, 4492, 4415,     8.000
     &   4294, 4125, 3915, 3671, 3413, 3141, 2848, 2558, 2293, 2053,     8.000
     &   1839, 1658, 1483, 1303, 1123,  927/                             8.000
!
      data iftabk, tsave, vsave / .false., 0.0, -1.0 /
!
      data np / nnp /
!
      data (tabp(ip), ip = 1, nnp) /  
     & -2.000, -1.500, -1.000, -0.500,  0.000,  0.500,  1.000,  1.250, 
     &  1.500,  1.750,  2.000,  2.250,  2.500,  2.750,  3.000,  3.200,
     &  3.400,  3.600,  3.800,  4.000,  4.200,  4.400,  4.600,  4.800, 
     &  5.000,  5.200,  5.400,  5.600,  5.800,  6.000,  6.250,  6.500, 
     &  6.750,  7.000,  7.250,  7.500,  7.750,  8.000                 /
!
      data nt / nnt /
!
      data (tabt(it), it = 1, nnt) /
     &  3.320,  3.340,  3.360,  3.380,  3.400,  3.420,  3.440,  3.460, 
     &  3.480,  3.500,  3.520,  3.540,  3.560,  3.580,  3.600,  3.620, 
     &  3.640,  3.660,  3.680,  3.700,  3.730,  3.760,  3.790,  3.820, 
     &  3.850,  3.880,  3.910,  3.940,  3.970,  4.000,  4.050,  4.100, 
     &  4.150,  4.200,  4.250,  4.300,  4.350,  4.400,  4.450,  4.500, 
     &  4.550,  4.600,  4.650,  4.700,  4.750,  4.800,  4.850,  4.900, 
     &  4.950,  5.000,  5.050,  5.100,  5.150,  5.200,  5.250,  5.300 /
!
      data nv / nnv /
!
      data tabv / 0.0d5, 1.0d5, 2.0d5, 4.0d5, 8.0d5 /
!
!---------------------------- EXECUTION --------------------------------
!
      exp10(x) = exp(x * 2.30258509299405d0)
!
      if (.not. iftabk) then
         iftabk = .true.
!
         do iv = 1, nnv
!
            do ip = 1 , nnp
!
               do it = 1, nnt
                  tabkap(it, ip, iv) = dble(ktab(it, ip, iv)) * 0.001
               end do
!
            end do
!
         end do
!
      end if
!
      if ((t .ne. tsave) .or. (v .ne. vsave)) then
         tsave = t
         vsave = v
         tlog = log10(t)
!
         do ip = 1, nnp
!
            do iv = 1, nnv
               idum = map1(tabt, tabkap(1,ip,iv), nt, tlog,
     &                     pkapv(iv,ip), 1)
            end do
!
            idum = map1(tabv, pkapv(1,ip), nv, v, pkap(ip), 1)
         end do
!
      end if
!
      idum = map1(tabp, pkap, np, log10(p), ablog, 1)
      rosstab = exp10(ablog)
      end
!
!!!!!!!!!!! E N D   O F   F U N C T I O N   R O S S T A B !!!!!!!!!!!!!!
!
      subroutine ttaup (t, tau, abstd, ptotal, p, prad, pturb, vturb, 
     &                  grav, numtau)
c     implicit undefined (a-z)
      implicit none
!
!...  TAUS MUST HAVE LOG SPACINGS
!...  ASSUMES PTOTAL = PRAD = PTURB = 0 WHEN RHOX = 0
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision  abstd(*), ptotal(*), p(*), prad(*), pturb(*), 
     &                  t(*), tau(*), vturb(*), grav
      integer  numtau
!
!------------------------------- EXTERNALS -----------------------------
!
      external rosstab, introssk 
!
!-------------------------- LOCAL VARIABLES --------------------------
!
      double precision  dlgtau, dplog, dplog1, dplog2, dplog3, error, 
     &                  plog, plog1, plog2, plog3, plog4, pnew, rosstab
      double precision introssk

      integer  i, j, n, maxit 
!
!------------------------------- EXECUTION -----------------------------
!
      maxit = 10000
      if (numtau .gt. 1) then
!
!....    THIS INITIALIZES THE INTEGRATION AT THE SURFACE
!
         dlgtau = log(tau(2) / tau(1))
         plog1  = 0.0d0
         plog2  = 0.0d0
         plog3  = 0.0d0
         dplog1 = 0.0d0
         dplog2 = 0.0d0
!
!....    ASSUME CONSTANT OPACITY NEAR THE SURFACE.  FIRST GUESS = 0.1
!
         abstd(1) = 0.1
         if (prad(1) .gt. 0.0d0) abstd(1) = min(0.1d0, grav * tau(1) / 
     &        prad(1) / 2.0d0)
      end if
!
      do j = 1, numtau
         error = 1.0
         n = 0
!
         if (j .gt. 4) then
            plog = (3.0 * plog4 + 8.0 * dplog1 - 4.0 * dplog2 + 
     &               8.0 * dplog3) / 3.0
!
         else if (j .gt. 1) then
            plog = plog1 + dplog1
!
         else
            plog = log(grav / abstd(1) * tau(1))
         end if
!
         ptotal(j) = exp(plog)
!
         do while (n .lt. maxit  .and. error .gt. 0.00005)
            n = n + 1
            p(j) = ptotal(j) + (prad(1) - prad(j)) - pturb(j)
!
            if (p(j) .le. 0.0d0) then
               write (16, 9000) ' pressure le 0.0 at depth', j,
     &            'j', 'p(j)', 'ptotal(j)', 'prad(j)', 
     &                  'pturb(j)', 'abstd(j)',
     &            (i, p(i), ptotal(i), prad(i),
     &                pturb(i), abstd(i), i = 1, j)
9000           format (a, i5/
     &                 4x, a, 8x, a, 9x, a, 7x, a, 7x, a, 8x, a /
     &                 (i5, 1p5e15.4))
               stop
!
            else
!               abstd(j) = rosstab(t(j), p(j), vturb(j))
               abstd(j) = introssk(t(j), p(j))
!
               if (numtau .eq. 1) then
                  error = 0.0d0
!
               else
                  dplog = grav / abstd(j) * tau(j) / ptotal(j) * 
     &                    dlgtau
!
                  if (j .gt. 4) then
                     pnew = (126.0 * plog1 - 14.0 * plog3 + 9.0 *
     &                       plog4 + 42.0 * dplog + 108.0 * dplog1 -
     &                       54.0 * dplog2 + 24.0 * dplog3) / 121.0
!
                  else if (j .gt. 1) then
                     pnew = (plog + 2.0 * plog1 + dplog + dplog1) / 3.0
!
                  else
                     pnew = log (grav / abstd(1) * tau(1))
                  end if
!
                  error = abs(pnew - plog)
                  plog = 0.5 * (pnew + plog)
                  ptotal(j) = exp(plog)
               end if
!
            end if
!
         end do
!
         if (numtau .gt. 1) then
!
            if (n .lt. maxit ) then
               plog4 = plog3
               plog3 = plog2
               plog2 = plog1
               plog1 = plog
               dplog3 = dplog2
               dplog2 = dplog1
               dplog1 = dplog
!
            else
               print*, 'WARNING, ttaup did not converge!'

               write (16, 9001) ' n = 1000 at depth', j,
     &            'j', 'p(j)', 'ptotal(j)', 'pturb(j)', 'abstd(j)',
     &            (i, p(i), ptotal(i), prad(i), pturb(i), abstd(i), 
     &             i = 1, j)
9001           format (a, i5 /
     &                 4x, a, 8x, a, 9x, a,
     &                    7x, 'prad(j)', 7x, a, 8x, a /
     &                    (i5, 1p5e15.4))
               stop
            end if
!
         end if
!
      end do
!
      end
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   T T A U P !!!!!!!!!!!!!!

      subroutine init_freset(nsub, nfreq, nt, np,  ifstandart, ivar)
       use types
       use atlcomm

       implicit none

       integer, intent(in) ::   nfreq, nt, np, nsub, ivar  
       logical, intent(in) :: ifstandart  
       integer :: ifre



       if (ifstandart) then 

         call set_freq(nsub, nfreq, nt, np, ivar) 


!!.... VALUES ARE 1/60, 2/60, 3/60, 6/60
         do ifre =1, nfreq

         sbwt(1,ifre)=  1.66666667d-2
         sbwt(2,ifre)=  3.33333333d-2
         sbwt(3,ifre)=  5.00000000d-2
         sbwt(4,ifre)=  1.00000000d-1
        sbwt(5,ifre)=  1.00000000d-1
        sbwt(6,ifre)=  1.00000000d-1
        sbwt(7,ifre)=  1.00000000d-1
        sbwt(8,ifre)=  1.00000000d-1
        sbwt(9,ifre)=  1.00000000d-1
        sbwt(10,ifre)=  1.00000000d-1
        sbwt(11,ifre)=  1.00000000d-1
        sbwt(12,ifre)=  1.00000000d-1

         end do

         tabp(1:np)= (/ -4.0d0, -3.5d0, -3.0d0, -2.5d0, -2.0d0, -1.5d0,
     &             -1.0d0, -0.5d0, 0.0d0, 0.5d0, 1.0d0,
     &             1.5d0,  2.0d0,  2.5d0,  3.0d0, 3.5d0, 4.0d0, 4.5d0,
     &             5.0d0,  5.5d0,  6.0d0,  6.5d0, 7.0d0, 7.5d0, 8.0d0 /)
!
         tabt(1:np) =  (/
     &  3.30d0, 3.32d0, 3.34d0, 3.36d0, 3.38d0, 3.40d0, 3.42d0, 3.44d0,
     &  3.46d0,
     &  3.48d0, 3.50d0, 3.52d0, 3.54d0, 3.56d0, 3.58d0, 3.60d0, 3.62d0,
     &  3.64d0, 3.66d0, 3.68d0, 3.70d0, 3.73d0, 3.76d0, 3.79d0, 3.82d0,
     &  3.85d0, 3.88d0, 3.91d0, 3.94d0, 3.97d0, 4.00d0, 4.05d0, 4.10d0,
     &  4.15d0, 4.20d0, 4.25d0, 4.30d0, 4.35d0, 4.40d0, 4.45d0, 4.50d0,
     &  4.55d0, 4.60d0, 4.65d0, 4.70d0, 4.75d0, 4.80d0, 4.85d0, 4.90d0,
     &  4.95d0, 5.00d0, 5.05d0, 5.10d0, 5.15d0, 5.20d0, 5.25d0, 5.30d0/)

         inifreset(6:95) = (/
     &        9.09,       9.35,       9.61,       9.77,       9.96,
     &       10.20,      10.38,      10.56,      10.77,      11.04,
     &       11.40,      11.78,      12.13,      12.48,      12.71,
     &       12.84,      13.05,      13.24,      13.39,      13.66,
     &       13.98,      14.33,      14.72,      15.10,      15.52,
     &       15.88,      16.20,      16.60,      17.03,      17.34,
     &       17.68,      18.02,      18.17,      18.61,      19.10,
     &       19.39,      19.84,      20.18,      20.50,      21.05,
     &       21.62,      21.98,      22.30,      22.68,      23.00,
     &       23.40,      24.00,      24.65,      25.24,      25.68,
     &       26.00,      26.40,      26.85,      27.35,      27.85,
     &       28.40,      29.0,       29.6,       30.1,       30.8,
     &       31.8,       32.8,       33.8,       34.8,       35.7,
     &       36.6,       37.5,       38.5,       39.5,       40.5,
     &       41.4,       42.2,       43.0,       44.1,       45.1,
     &       46.0,       47.0,       48.0,       49.0,       50.0,
     &       50.6,       51.4,       53.0,       55.0,       56.7,
     &       58.5,       60.5,       62.5,       64.5,       66.3   /)
!
        inifreset( 96: 185) = (/
     &       68.0,       70.0,       71.6,       73.0,       75.0,
     &       77.0,       79.0,       81.0,       83.0,       85.0,
     &       87.0,       89.0,       90.6,       92.6,       96.0,
     &      100.0,      104.0,      108.0,      111.5,      114.5,
     &      118.0,      122.0,      126.0,      130.0,      134.0,
     &      138.0,      142.0,      146.0,      150.0,      154.0,
     &      160.0,      165.0,      169.0,      173.0,      177.5,
     &      182.0,      186.0,      190.5,      195.0,      200.0,
     &      204.5,      208.5,      212.5,      217.5,      222.5,
     &      227.5,      232.5,      237.5,      242.5,      248.0,
     &      253.0,      257.5,      262.5,      267.5,      272.5,
     &      277.5,      282.5,      287.5,      295.0,      305.0,
     &      315.0,      325.0,      335.0,      345.0,      355.0,
     &      362.0,      367.0,      375.0,      385.0,      395.0,
     &      405.0,      415.0,      425.0,      435.0,      455.0,
     &      465.0,      475.0,      485.0,      495.0,      505.0,
     &      515.0,      525.0,      535.0,      545.0,      555.0,
     &      565.0,      575.0,      585.0,      595.0,      605.0  /)
!
        inifreset( 186: 275) = (/ 
     &      615.0,      625.0,      635.0,      645.0,      655.0,
     &      665.0,      675.0,      685.0,      695.0,      705.0,
     &      715.0,      725.0,      735.0,      745.0,      755.0,
     &      765.0,      775.0,      785.0,      795.0,      805.0,
     &      815.0,      825.0,      835.0,      845.0,      855.0,
     &      865.0,      875.0,      885.0,      895.0,      905.0,
     &      915.0,      925.0,      935.0,      945.0,      955.0,
     &      965.0,      975.0,      985.0,      995.0,     1012.5,
     &     1037.5,     1062.5,     1087.5,     1112.5,     1137.5,
     &     1162.5,     1187.5,     1212.5,     1237.5,     1262.5,
     &     1287.5,     1312.5,     1337.5,     1362.5,     1387.5,
     &     1412.5,     1442.0,     1467.0,     1487.5,     1512.5,
     &     1537.5,     1562.5,     1587.5,     1620.0,     1660.0,
     &     1700.0,     1740.0,     1780.0,     1820.0,     1860.0,
     &     1900.0,     1940.0,     1980.0,     2025.0,     2075.0,
     &     2125.0,     2175.0,     2225.0,     2265.0,     2290.0,
     &     2325.0,     2375.0,     2425.0,     2475.0,     2525.0,
     &     2575.0,     2625.0,     2675.0,     2725.0,     2775.0  /)
!
      inifreset( 276: 342) = (/ 
     &     2825.0,     2875.0,     2925.0,     2975.0,     3025.0,
     &     3075.0,     3125.0,     3175.0,     3240.0,     3340.0,
     &     3450.0,     3550.0,     3650.0,     3750.0,     3850.0,
     &     3950.0,     4050.0,     4150.0,     4250.0,     4350.0,
     &     4450.0,     4550.0,     4650.0,     4750.0,     4850.0,
     &     4950.0,     5050.0,     5150.0,     5250.0,     5350.0,
     &     5450.0,     5550.0,     5650.0,     5750.0,     5850.0,
     &     5950.0,     6050.0,     6150.0,     6250.0,     6350.0,
     &     6500.0,     6700.0,     6900.0,     7100.0,     7300.0,
     &     7500.0,     7700.0,     7900.0,     8100.0,     8300.0,
     &     8500.0,     8700.0,     8900.0,     9100.0,     9300.0,
     &     9500.0,     9700.0,     9900.0,    10000.0,    20000.0,
     &    40000.0,    60000.0,    80000.0,   100000.0,   120000.0,
     &   140000.0,   160000.0                                      /)
!
      inifreset(343: 432) = (/ 
     &        9.09,       9.35,       9.61,       9.77,       9.96,
     &       10.20,      10.38,      10.56,      10.77,      11.04,
     &       11.40,      11.78,      12.13,      12.48,      12.71,
     &       12.84,      13.05,      13.24,      13.39,      13.66,
     &       13.98,      14.33,      14.72,      15.10,      15.52,
     &       15.88,      16.20,      16.60,      17.03,      17.34,
     &       17.68,      18.02,      18.17,      18.61,      19.10,
     &       19.39,      19.84,      20.18,      20.50,      21.05,
     &       21.62,      21.98,      22.30,      22.68,      23.00,
     &       23.40,      24.00,      24.65,      25.24,      25.68,
     &       26.00,      26.40,      26.85,      27.35,      27.85,
     &       28.40,      29.0,       29.6,       30.1,       30.8,
     &       31.8,       32.8,       33.8,       34.8,       35.7,
     &       36.6,       37.5,       38.5,       39.5,       40.5,
     &       41.4,       42.2,       43.0,       44.1,       45.1,
     &       46.0,       47.0,       48.0,       49.0,       50.0,
     &       50.6,       51.2,       52.0,       53.0,       54.0,
     &       55.0,       56.0,       57.0,       58.0,       59.0 /)
!
      inifreset(433: 522) =(/ 
     &       60.0,       61.0,       62.0,       63.0,       64.0,
     &       65.0,       65.8,       66.5,       67.5,       68.5,
     &       69.5,       70.5,       71.6,       72.6,       73.5,
     &       74.5,       75.5,       76.5,       77.5,       78.5,
     &       79.5,       80.5,       81.5,       82.5,       83.5,
     &       84.5,       85.5,       86.5,       87.5,       88.5,
     &       89.5,       90.5,       91.5,       92.5,       93.5,
     &       94.5,       95.5,       96.5,       97.5,       98.5,
     &       99.5,      100.5,      101.5,      102.5,      103.5,
     &      104.5,      105.5,      106.5,      107.5,      108.5,
     &      109.5,      110.5,      111.5,      112.5,      113.5,
     &      114.5,      115.5,      116.5,      117.5,      118.5,
     &      119.5,      120.5,      121.5,      122.5,      123.5,
     &      124.5,      125.5,      126.5,      127.5,      128.5,
     &      129.5,      130.5,      131.5,      132.5,      133.5,
     &      134.5,      135.5,      136.5,      137.5,      138.5,
     &      139.5,      140.5,      141.5,      142.5,      143.5,
     &      144.2,      144.7,      145.5,      146.5,      147.5  /)

 
      inifreset(523: 612) =(/ 
     &      148.5,      149.5,      150.5,      151.2,      151.7,
     &      152.5,      153.5,      154.5,      155.5,      156.5,
     &      157.5,      158.5,      159.5,      160.5,      161.5,
     &      162.5,      163.5,      164.5,      165.5,      166.5,
     &      167.2,      167.7,      168.5,      169.5,      170.5,
     &      171.5,      172.5,      173.5,      174.5,      175.5,
     &      176.5,      177.5,      178.5,      179.5,      180.5,
     &      181.5,      182.5,      183.5,      184.5,      185.5,
     &      186.5,      187.5,      188.5,      189.5,      190.5,
     &      191.5,      192.5,      193.5,      194.5,      195.5,
     &      196.7,      197.6,      198.4,      199.5,      200.5,
     &      201.5,      202.5,      203.5,      204.5,      205.5,
     &      206.5,      207.4,      207.8,      208.5,      209.5,
     &      210.5,      211.5,      212.5,      213.5,      214.5,
     &      215.5,      216.5,      217.5,      218.5,      219.5,
     &      220.5,      221.5,      222.5,      223.5,      224.5,
     &      225.5,      226.5,      227.5,      228.5,      229.5,
     &      230.5,      231.5,      232.5,      233.5,      234.5  /)
!
      inifreset(613: 702) = (/ 
     &      235.5,      236.5,      237.5,      238.5,      239.5,
     &      240.5,      241.5,      242.5,      243.5,      244.5,
     &      245.5,      246.5,      247.5,      248.5,      249.5,
     &      250.5,      251.3,      251.8,      252.5,      253.5,
     &      254.5,      255.5,      256.5,      257.5,      258.5,
     &      259.5,      260.5,      261.5,      262.5,      263.5,
     &      264.5,      265.5,      266.5,      267.5,      268.5,
     &      269.5,      270.5,      271.5,      272.5,      273.5,
     &      274.5,      275.5,      276.5,      277.5,      278.5,
     &      279.5,      280.5,      281.5,      282.5,      283.5,
     &      284.5,      285.5,      286.5,      287.5,      288.5,
     &      289.5,      291.0,      293.0,      295.0,      297.0,
     &      299.0,      301.0,      303.0,      305.0,      307.0,
     &      309.0,      311.0,      313.0,      315.0,      317.0,
     &      319.0,      321.0,      323.0,      325.0,      327.0,
     &      329.0,      331.0,      333.0,      335.0,      337.0,
     &      339.0,      341.0,      343.0,      345.0,      347.0,
     &      349.0,      351.0,      353.0,      355.0,      357.0  /)
!
      inifreset(703: 792) =(/ 
     &      359.0,      361.0,      363.0,      364.0,      365.0,
     &      367.0,      369.0,      371.0,      373.0,      375.0,
     &      377.0,      379.0,      381.0,      383.0,      385.0,
     &      387.0,      389.0,      391.0,      393.0,      395.0,
     &      397.0,      399.0,      401.0,      403.0,      405.0,
     &      407.0,      409.0,      411.0,      413.0,      415.0,
     &      417.0,      419.0,      421.0,      423.0,      425.0,
     &      427.0,      429.0,      431.0,      433.0,      435.0,
     &      437.0,      439.0,      441.0,      443.0,      445.0,
     &      447.0,      449.0,      451.0,      453.0,      455.0,
     &      457.0,      459.0,      461.0,      463.0,      465.0,
     &      467.0,      469.0,      471.0,      473.0,      475.0,
     &      477.0,      479.0,      481.0,      483.0,      485.0,
     &      487.0,      489.0,      491.0,      493.0,      495.0,
     &      497.0,      499.0,      501.0,      503.0,      505.0,
     &      507.0,      509.0,      511.0,      513.0,      515.0,
     &      517.0,      519.0,      521.0,      523.0,      525.0,
     &      527.0,      529.0,      531.0,      533.0,      535.0  /)
!
      inifreset(793: 882) = (/
     &      537.0,      539.0,      541.0,      543.0,      545.0,
     &      547.0,      549.0,      551.0,      553.0,      555.0,
     &      557.0,      559.0,      561.0,      563.0,      565.0,
     &      567.0,      569.0,      571.0,      573.0,      575.0,
     &      577.0,      579.0,      581.0,      583.0,      585.0,
     &      587.0,      589.0,      591.0,      593.0,      595.0,
     &      597.0,      599.0,      601.0,      603.0,      605.0,
     &      607.0,      609.0,      611.0,      613.0,      615.0,
     &      617.0,      619.0,      621.0,      623.0,      625.0,
     &      627.0,      629.0,      631.0,      633.0,      635.0,
     &      637.0,      639.0,      641.0,      643.0,      645.0,
     &      647.0,      649.0,      651.0,      653.0,      655.0,
     &      657.0,      659.0,      661.0,      663.0,      665.0,
     &      667.0,      669.0,      671.0,      673.0,      675.0,
     &      677.0,      679.0,      681.0,      683.0,      685.0,
     &      687.0,      689.0,      691.0,      693.0,      695.0,
     &      697.0,      699.0,      701.0,      703.0,      705.0,
     &      707.0,      709.0,      711.0,      713.0,      715.0  /)
!
      inifreset(883: 972) =(/ 
     &      717.0,      719.0,      721.0,      723.0,      725.0,
     &      727.0,      729.0,      731.0,      733.0,      735.0,
     &      737.0,      739.0,      741.0,      743.0,      745.0,
     &      747.0,      749.0,      751.0,      753.0,      755.0,
     &      757.0,      759.0,      761.0,      763.0,      765.0,
     &      767.0,      769.0,      771.0,      773.0,      775.0,
     &      777.0,      779.0,      781.0,      783.0,      785.0,
     &      787.0,      789.0,      791.0,      793.0,      795.0,
     &      797.0,      799.0,      801.0,      803.0,      805.0,
     &      807.0,      809.0,      811.0,      813.0,      815.0,
     &      817.0,      819.0,      821.0,      823.0,      825.0,
     &      827.0,      829.0,      831.0,      833.0,      835.0,
     &      837.0,      839.0,      841.0,      843.0,      845.0,
     &      847.0,      849.0,      851.0,      853.0,      855.0,
     &      857.0,      859.0,      861.0,      863.0,      865.0,
     &      867.0,      869.0,      871.0,      873.0,      875.0,
     &      877.0,      879.0,      881.0,      883.0,      885.0,
     &      887.0,      889.0,      891.0,      893.0,      895.0  /) 
!
       inifreset( 973: 1062)  = (/ 
     &      897.0,      899.0,      901.0,      903.0,      905.0,
     &      907.0,      909.0,      911.0,      913.0,      915.0,
     &      917.0,      919.0,      921.0,      923.0,      925.0,
     &      927.0,      929.0,      931.0,      933.0,      935.0,
     &      937.0,      939.0,      941.0,      943.0,      945.0,
     &      947.0,      949.0,      951.0,      953.0,      955.0,
     &      957.0,      959.0,      961.0,      963.0,      965.0,
     &      967.0,      969.0,      971.0,      973.0,      975.0,
     &      977.0,      979.0,      981.0,      983.0,      985.0,
     &      987.0,      989.0,      991.0,      993.0,      995.0,
     &      997.0,      999.0,     1002.5,     1007.5,     1012.5,
     &     1017.5,     1022.5,     1027.5,     1032.5,     1037.5,
     &     1042.5,     1047.5,     1052.5,     1057.5,     1062.5,
     &     1067.5,     1072.5,     1077.5,     1082.5,     1087.5,
     &     1092.5,     1097.5,     1102.5,     1107.5,     1112.5,
     &     1117.5,     1122.5,     1127.5,     1132.5,     1137.5,
     &     1142.5,     1147.5,     1152.5,     1157.5,     1162.5,
     &     1167.5,     1172.5,     1177.5,     1182.5,     1187.5  /)

!
       inifreset(1063: 1152) =(/ 
     &     1192.5,     1197.5,     1202.5,     1207.5,     1212.5,
     &     1217.5,     1222.5,     1227.5,     1232.5,     1237.5,
     &     1242.5,     1247.5,     1252.5,     1257.5,     1262.5,
     &     1267.5,     1272.5,     1277.5,     1282.5,     1287.5,
     &     1292.5,     1297.5,     1302.5,     1307.5,     1312.5,
     &     1317.5,     1322.5,     1327.5,     1332.5,     1337.5,
     &     1342.5,     1347.5,     1352.5,     1357.5,     1362.5,
     &     1367.5,     1372.5,     1377.5,     1382.5,     1387.5,
     &     1392.5,     1397.5,     1402.5,     1407.5,     1412.5,
     &     1417.5,     1422.5,     1427.5,     1432.5,     1437.5,
     &     1442.5,     1447.5,     1452.5,     1457.0,     1462.0,
     &     1467.5,     1472.5,     1477.5,     1482.5,     1487.5,
     &     1492.5,     1497.5,     1502.5,     1507.5,     1512.5,
     &     1517.5,     1522.5,     1527.5,     1532.5,     1537.5,
     &     1542.5,     1547.5,     1552.5,     1557.5,     1562.5,
     &     1567.5,     1572.5,     1577.5,     1582.5,     1587.5,
     &     1592.5,     1597.5,     1605.0,     1615.0,     1625.0,
     &     1635.0,     1645.0,     1655.0,     1665.0,     1675.0  /)
!
      inifreset(1153: 1242) = (/ 
     &     1685.0,     1695.0,     1705.0,     1715.0,     1725.0,
     &     1735.0,     1745.0,     1755.0,     1765.0,     1775.0,
     &     1785.0,     1795.0,     1805.0,     1815.0,     1825.0,
     &     1835.0,     1845.0,     1855.0,     1865.0,     1875.0,
     &     1885.0,     1895.0,     1905.0,     1915.0,     1925.0,
     &     1935.0,     1945.0,     1955.0,     1965.0,     1975.0,
     &     1985.0,     1995.0,     2005.0,     2015.0,     2025.0,
     &     2035.0,     2045.0,     2055.0,     2065.0,     2075.0,
     &     2085.0,     2095.0,     2105.0,     2115.0,     2125.0,
     &     2135.0,     2145.0,     2155.0,     2165.0,     2175.0,
     &     2185.0,     2195.0,     2205.0,     2215.0,     2225.0,
     &     2235.0,     2245.0,     2255.0,     2265.0,     2275.0,
     &     2285.0,     2295.0,     2305.0,     2315.0,     2325.0,
     &     2335.0,     2345.0,     2355.0,     2365.0,     2375.0,
     &     2385.0,     2395.0,     2405.0,     2415.0,     2425.0,
     &     2435.0,     2445.0,     2455.0,     2465.0,     2475.0,
     &     2485.0,     2495.0,     2505.0,     2515.0,     2525.0,
     &     2535.0,     2545.0,     2555.0,     2565.0,     2575.0  /) 
!
      inifreset( 1243: 1332) = (/ 
     &     2585.0,     2595.0,     2605.0,     2615.0,     2625.0,
     &     2635.0,     2645.0,     2655.0,     2665.0,     2675.0,
     &     2685.0,     2695.0,     2705.0,     2715.0,     2725.0,
     &     2735.0,     2745.0,     2755.0,     2765.0,     2775.0,
     &     2785.0,     2795.0,     2805.0,     2815.0,     2825.0,
     &     2835.0,     2845.0,     2855.0,     2865.0,     2875.0,
     &     2885.0,     2895.0,     2905.0,     2915.0,     2925.0,
     &     2935.0,     2945.0,     2955.0,     2965.0,     2975.0,
     &     2985.0,     2995.0,     3005.0,     3015.0,     3025.0,
     &     3035.0,     3045.0,     3055.0,     3065.0,     3075.0,
     &     3085.0,     3095.0,     3105.0,     3115.0,     3125.0,
     &     3135.0,     3145.0,     3155.0,     3165.0,     3175.0,
     &     3185.0,     3195.0,     3210.0,     3230.0,     3250.0,
     &     3270.0,     3290.0,     3310.0,     3330.0,     3350.0,
     &     3370.0,     3390.0,     3410.0,     3430.0,     3450.0,
     &     3470.0,     3490.0,     3510.0,     3530.0,     3550.0,
     &     3570.0,     3590.0,     3610.0,     3630.0,     3650.0,
     &     3670.0,     3690.0,     3710.0,     3730.0,     3750.0  /)



      inifreset(1333: 1422) = (/
     &     3770.0,     3790.0,     3810.0,     3830.0,     3850.0,
     &     3870.0,     3890.0,     3910.0,     3930.0,     3950.0,
     &     3970.0,     3990.0,     4010.0,     4030.0,     4050.0,
     &     4070.0,     4090.0,     4110.0,     4130.0,     4150.0,
     &     4170.0,     4190.0,     4210.0,     4230.0,     4250.0,
     &     4270.0,     4290.0,     4310.0,     4330.0,     4350.0,
     &     4370.0,     4390.0,     4410.0,     4430.0,     4450.0,
     &     4470.0,     4490.0,     4510.0,     4530.0,     4550.0,
     &     4570.0,     4590.0,     4610.0,     4630.0,     4650.0,
     &     4670.0,     4690.0,     4710.0,     4730.0,     4750.0,
     &     4770.0,     4790.0,     4810.0,     4830.0,     4850.0,
     &     4870.0,     4890.0,     4910.0,     4930.0,     4950.0,
     &     4970.0,     4990.0,     5010.0,     5030.0,     5050.0,
     &     5070.0,     5090.0,     5110.0,     5130.0,     5150.0,
     &     5170.0,     5190.0,     5210.0,     5230.0,     5250.0,
     &     5270.0,     5290.0,     5310.0,     5330.0,     5350.0,
     &     5370.0,     5390.0,     5410.0,     5430.0,     5450.0,
     &     5470.0,     5490.0,     5510.0,     5530.0,     5550.0  /)
!
      inifreset(1423: 1512) = (/ 
     &     5570.0,     5590.0,     5610.0,     5630.0,     5650.0,
     &     5670.0,     5690.0,     5710.0,     5730.0,     5750.0,
     &     5770.0,     5790.0,     5810.0,     5830.0,     5850.0,
     &     5870.0,     5890.0,     5910.0,     5930.0,     5950.0,
     &     5970.0,     5990.0,     6010.0,     6030.0,     6050.0,
     &     6070.0,     6090.0,     6110.0,     6130.0,     6150.0,
     &     6170.0,     6190.0,     6210.0,     6230.0,     6250.0,
     &     6270.0,     6290.0,     6310.0,     6330.0,     6350.0,
     &     6370.0,     6390.0,     6420.0,     6460.0,     6500.0,
     &     6540.0,     6580.0,     6620.0,     6660.0,     6700.0,
     &     6740.0,     6780.0,     6820.0,     6860.0,     6900.0,
     &     6940.0,     6980.0,     7020.0,     7060.0,     7100.0,
     &     7140.0,     7180.0,     7220.0,     7260.0,     7300.0,
     &     7340.0,     7380.0,     7420.0,     7460.0,     7500.0,
     &     7540.0,     7580.0,     7620.0,     7660.0,     7700.0,
     &     7740.0,     7780.0,     7820.0,     7860.0,     7900.0,
     &     7940.0,     7980.0,     8020.0,     8060.0,     8100.0,
     &     8140.0,     8180.0,     8220.0,     8260.0,     8300.0  /)
!
      inifreset(1513: 1563) = (/ 
     &     8340.0,     8380.0,     8420.0,     8460.0,     8500.0,
     &     8540.0,     8580.0,     8620.0,     8660.0,     8700.0,
     &     8740.0,     8780.0,     8820.0,     8860.0,     8900.0,
     &     8940.0,     8980.0,     9020.0,     9060.0,     9100.0,
     &     9140.0,     9180.0,     9220.0,     9260.0,     9300.0,
     &     9340.0,     9380.0,     9420.0,     9460.0,     9500.0,
     &     9540.0,     9580.0,     9620.0,     9660.0,     9700.0,
     &     9740.0,     9780.0,     9820.0,     9860.0,     9900.0,
     &     9940.0,     9980.0,    10020.0,    20000.0,    40000.0,
     &    60000.0,    80000.0,   100000.0,   120000.0,   140000.0,
     &   160000.0                                                  /)

      else

       call set_freq(nsub, nfreq, nt, np, ivar)

      end if 


      end subroutine 

      subroutine print_summary(mode1) 

      use types
      use atlcomm
!
      implicit none

      include 'common.teffbl'
      include 'common.junkbl'
      include 'common.convbl'
      include 'common.ifopbl'
      include 'common.steplg'
      include 'common.turbpr'
      include 'common.elembl'
      include 'common.xabund'


      integer, intent(in) :: mode1
      character(3) a, b
      real(kind=8) ::locteff, locglog  
      integer :: i , iz
      integer :: iop(20)



      if (mode1 .gt. 0) then 

       if (mode1 .eq. 1) print*, 'This run calculates Model'
       if (mode1 .eq. 2) print*, 'This run calculates surface fluxes'   
        
       iop = 0 

       do i = 1, 20
         if (ifop(i)) iop(i) = 1
       end do  

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
       if (mode1 .eq. 1) then 
         locteff = teff1
         locglog = log10(gnew) 
       else 
         locteff = teff
         locglog = glog
       endif 


!
         write (6, 552) 'teff ', locteff, 'gravity', locglog, wlte,
     &      ' opacity ifop', iop,
     &      ' convection ', a, mixlth, ' turbulence ', b,
     &      trbfdg, trbpow, trbsnd, trbcon,
     &      'abundance scale ', xscale, ' abundance change',
     &                                  (iz, abund(iz), iz = 1, 2),
     &      (' abundance change',
     &      (iz, abund(iz), iz = i, i+5), i = 3, 98, 6),
     &       ' abundance change', (iz, abund(iz), iz = 99, 99)
552      format (a, f7.0, 2x, a, f8.5, 1x, a /
     &         a, 20i2 /
     &         a, a, f6.2, a, a, 4f6.2 /
     &         a, f9.5, a, 2(i2, f8.5) /
     &        (a, 6(i3, f7.2)))

      if(mode1 .eq. 1) then 
       print*, 'The model is calculated on tau-grid:'
       print*, 'starting tau ', tau1lg, 'step = ', steplg
       print*, 'using nrhox ', krhox, 'finishing at ', taulalg
 

      endif  


      else

        print*, 'This run calculates ODF'

        if (ifkbin) then 
         print*, 'Standard Kurucz binning'
        else
         print*, 'User defined binning'
         if (ifilter .eq. 1 ) then 
           print*, 'Filters are used'
         endif  
        endif  
        
        if (ameanODF .eq. 1) print*, 'Arithmetic mean for ODF is used'
        if (ameanODF .eq. 0) print*, 'Geometric mean for ODF is used '
 
         write (6, 550) 'abundance scale ', xscale, ' abundance change',
     &                                  (iz, abund(iz), iz = 1, 2),
     &      (' abundance change',
     &      (iz, abund(iz), iz = i, i+5), i = 3, 98, 6),
     &       ' abundance change', (iz, abund(iz), iz = 99, 99)
550      format ( a, f9.5, a, 2(i2, f8.5) /
     &        (a, 6(i3, f7.2)))    

      endif 









      end subroutine print_summary
