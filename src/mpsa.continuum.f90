  subroutine calc_continuum 


    use types
    use atlcomm
    use kappa_cal
    use population
 
    implicit none

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




!
!--------------------------- LOCAL VARIABLES --------------------------
!
!--- for comparison with cotinuumn in RH code 
!      double precision fit(5), xo, yo, opaca(maxd) 

      double precision  contin, freq15, part(maxd, 6), rco, & 
     &                  rcowt, stepwt, sumwt, wave, x, rosstab, asixth
!     double precision exp10
      integer  i, j, mm,  mode1, n, nsteps, nu, nmod, imod
      logical  finish, more, readi0, stopfl
      logical otherselection

!
!-------------------------------  EXTERNALS ----------------------------
!
      external convec, high, josh,  & 
     &         output, outpt2, outpt3, outpt4, outpt5, & 
     &         radiap, radiap2, radiap3, readi0, ross, ross2, &  
     &         ross3, stateq, stateq2, stateq3, tcorr, tcorr2, tcorr3, & 
     &         turb, rosstab
      character lit*1



      integer max1, in, nt, nelem6, place   

!   implicit params
      integer  last1, nnnnu, nelem, ion, nu3, iT, iP 
      double precision save, edge,  eq, cutoff, frqlg, contmin 
      double precision freeff
      parameter (max1=maxmol+1)
      
      logical nosaha, checkinf, isinf 

      double precision xnfh2(maxd), xnfph2(maxd), xnfpco(maxd)
      real(kind=dp) xnfp(maxd,10,99),xnfpel(6,99),dopple(6,99)
      double precision dummyxn(maxd,6)
      double precision ablog(maxd)
      double precision a(377), continall(1131, maxd), frqedg(377) 
      double precision contabs(1131,maxd)
      double precision contscat(1131,maxd)

      real(kind=dp)  xnfpel4(6,99,maxd), dopple4(6,99, maxd), continall4(1131)
      real(kind=dp)  xnfdopmax(377,6,99), xnfdopcont 
      real(kind=dp), allocatable ::  contgrid(:,:,:,:)



!---- for netcdf files ----------------------------------------------!
       integer   ncid, ier, myrank, ntemp, nproc, nufreq, comm 
       character(80) :: name
!---------------------------------------------------------------------

      character card*80  
      DOUBLE PRECISION  WAVEBIG(329)
      DOUBLE PRECISION  BIGA(95),BIGB(125),BIGC(993)
      EQUIVALENCE (WAVEBIG(1),BIGA(1)),(WAVEBIG(96),BIGB(1))
      EQUIVALENCE (WAVEBIG(221),BIGC(1))

      BIGA=(/ & 
    &     8.97666,     9.2,         9.5,         9.71730,     9.81590, & 
    &    10.10901,    10.3,        10.46451,    10.65,       10.88545, & 
    &    11.2,        11.6,        11.95562,    12.3,        12.66565, & 
    &    12.74676,    12.93504,    13.16033,    13.31723,    13.46700, & 
    &    13.86253,    14.10462,    14.54782,    14.9,        15.3, & 
    &    15.74874,    16.00346,    16.4,        16.8,        17.25088, & 
    &    17.43258,    17.93660,    18.10114,    18.24438,    18.97796, & 
    &    19.22422,    19.55098,    20.11987,    20.23592,    20.8, & 
    &    21.3,        21.94783,    22.01983,    22.58070,    22.78377, & 
    &    23.2,        23.65221,    24.3,        25.00744,    25.46614, & 
    &    25.89239,    26.14316,    26.6,        27.1,        27.6, & 
    &    28.09263,    28.65871,    29.3,        29.95849,    30.26260, & 
    &    31.3,        32.3,        33.3,        34.3,        35.27931, &
    &    36.14100,    37.,         38.01600,    38.96642,    40., & 
    &    41.,         41.87992,    42.5,        43.5,        44.73260, & 
    &    45.5,        46.5,        47.5,        48.50178,    49.5, & 
    &    50.42590,    50.87630,    52.,         54.,         56., & 
    &    57.40614,    59.5,        61.5,        63.5,        65.53786, & 
    &    67.,         69.,         71.,         72.24011,    74. /) 

      BIGB=(/ & 
    &    76.,         78.,         80.,         82.,         84., & 
    &    86.,         88.,         90.,         91.17535,    94., & 
    &    98.,         102.,       106.,        110.03056,   113., & 
    &   116.,         120.,       123.92928,   128.,        132., & 
    &   136.,         140.,       144.43391,   148.,        151.43485, & 
    &   157.,         162.15072,  167.40282,   171.,        175., & 
    &   180.,         184.,       188.,        193.,        197.46990, & 
    &   202.,         207.13210,  210.,        215.,        220., & 
    &   225.,         230.,       235.,        240.,        245., & 
    &   251.12621,    255.,       260.,        265.,        270., & 
    &   275.,         280.,       285.,        290.,        300., & 
    &   310.,         320.,       330.,        340.,        350., & 
    &   360.,         364.70183,  370.,        380.,        390., & 
    &   400.,   410., 420., 430., 450., 460.,  470., 480.,  490., 500., & 
    &   510.,   520., 530., 540., 550., 560.,  570., 580.,  590., 600., & 
    &   610.,   620., 630., 640., 650., 660.,  670., 680.,  690., 700., & 
    &   710.,   720., 730., 740., 750., 760.,  770., 780.,  790., 800., & 
    &   810.,820.58271,830.,840., 850., 860.,  870., 880.,  890., 900., & 
    &   910.,   920., 930., 940., 950., 960.,  970., 980.,  990.,1000. /) 


      BIGC(1:109)=(/ & 
    & 1025.,1050.,1075.,1100.,1125.,1150.,1175.,1200.,1225.,1250., & 
    & 1275.,1300.,1325.,1350.,1375.,1400.,1425.,1458.81670,1475.,1500., & 
    & 1525.,1550.,1575.,1600.,1640.,1680.,1720.,1760.,1800.,1840., & 
    & 1880.,1920.,1960.,2000.,2050.,2100.,2150.,2200.,2250.,2279.40330, & 
    & 2300.,2350.,2400.,2450.,2500.,2550.,2600.,2650.,2700.,2750., & 
    & 2800.,2850.,2900.,2950.,3000.,3050.,3100.,3150.,3200.,3282.34320, & 
    & 3400.,3500.,3600.,3700.,3800.,3900.,4000.,4100.,4200.,4300., & 
    & 4400.,4500.,4600.,4700.,4800.,4900.,5000.,5100.,5200.,5300., &
    & 5400.,5500.,5600.,5700.,5800.,5900.,6000.,6100.,6200.,6300., & 
    & 6400.,6600.,6800.,7000.,7200.,7400.,7600.,7800.,8000.,8200., & 
    & 8400.,8600.,8800.,9000.,9200.,9400.,9600.,9800.,10000000. /) 
      BIGC(110) = 0.0d0 





     otherselection = .false.

     myrank = 0
     nufreq = 1131
 

       idmol=(/ 101.,  106.,  107.,  108.,  606.,  607.,  608., & 
     & 707.,  708., & 
     &   808.,  112.,  113.,  114.,  812.,  813.,  814.,  116.,  120., & 
     &   816.,  820.,  821.,  822.,  823.,  103.,  104.,  105.,  109., & 
     &   115.,  117.,  121.,  122.,  123.,  124.,  125.,  126.,106.01, &
     & 107.01,108.01,112.01,113.01,114.01,120.01,  408.,  508.,  815., &
     &   817.,  824.,  825.,  826.,10108.,60808.,10106.,60606.,  127., &
     &   128.,  129.,  827.,  828.,  829., 608.01 /) 
 
       momass= (/ 2.,   13.,   15.,   17.,   24.,   26.,   28.,   28., & 
     &   30., & 
     &    32.,   25.,   28.,   29.,   40.,   43.,   44.,   33.,   41., &
     &    48.,   56.,   61.,   64.,   67.,    8.,   10.,   12.,   20.,& 
     &    32.,   36.,   46.,   49.,   52.,   53.,   56.,   57.,   13., &
     &    15.,   17.,   25.,   28.,   29.,   41.,   25.,   27.,   47., &
     &    51.,   68.,   71.,   72.,   18.,   44.,   14.,   36.,   60., &
     &    59.,   64.,   75.,   74.,   79.,   28. /)

! READ in input header first, as later on you want to be able to allocate 
! the right size array for small and big wl bins (not subbins)!



!!
!!.... OPEN THE I/O FILES
!!
      open (unit = 15, file = fileODFinput, form = 'formatted', & 
     &      status = 'old', access = 'sequential')
      open (unit = 16, file = 'mpsa.print', form = 'formatted', &
     &      status = 'new', access = 'sequential')
      open (unit = 17, file = 'mpsa.punch', form = 'formatted', &
     &      status = 'new', access = 'sequential')
      open (unit = 18, file = 'mpsa.jnu', form = 'formatted', &
     &      status = 'new', access = 'sequential')
      open (unit = 33, file = filePTgrid, form = 'formatted', &
     &      status = 'old', access = 'sequential')

!---------------------------------------------------------------------
     
      open (unit = 77, file ='continuumall.dat', form = 'formatted', &
     &      status = 'new')

      open (unit = 88, file ='continuumabs.dat', form = 'formatted', &
     &      status = 'new')

      open (unit = 79, file ='continuumscat.dat', form = 'formatted', &
     &      status = 'new')


      open (unit = 44, file = './INPUT/continua.dat', form = 'formatted', &
     &      status = 'old', access = 'sequential') 

      more = .true.
      nosaha = .false.

!
      maxpow = 99
      last=80
      numcol=1
      in = 0

   
       do i = 1, 329

       wledge(i) = wavebig(i)
       cmedge(i) = 1.0d7/wavebig(i)
       frqedg(i) = c_nm/wavebig(i)

      end do

      in = 329



      numnu=0
      do i=1, in-1
       numnu= numnu +1
       freqset( numnu ) = abs(frqedg(i))/1.0000001
       numnu= numnu +1
       freqset( numnu ) = c_nm/(abs(wledge(i)) + abs(wledge(i+1)))*2.0d0
       numnu = numnu +1
       freqset( numnu ) = abs(frqedg(i+1))*1.0000001
 
      end do 
    
      nnnnu= numnu

!---------- only now read the model
      more = .true.
!-------------------------------------------------------------------
! call reading model (reads header and how many atmosphere models
      more = readi0 (mode1, 0, 1 )
         if(more) then
            read (33,*) nmod ! this will be the amount of T values 
            read (33,*) nrhox ! this will be the amount of p values 

            if (nrhox .gt. maxd ) then 
                write(*, *) 'amount of P values > maxd'
               stop
            end if  
 
            if(nmod .gt. maxd) then
               write( *,*) ' nmod = ', nmod, ', which is > maxd'
               stop
            end if
         end if


      allocate(contgrid(328,nrhox ,nmod ,3))


!     READ MODEL (first nmod T, then nrhox p values)
!     Because we need one atmosphere at a time keep T fixed
!     and calculate for all p values -> for not pretabulated 
!     calculations we need to adjust here!!! 
!-------------------------------------------------------!
            
               read (33, *) (tsave(j), j = 1, nmod)
!               read (33,*) (tsave(j), j=1, nrhox)
               read (33, *) (p(j), j=1, nrhox)
!-------------------------------------------------------!
! now we now dimensions of the calculations and can creat 
!  NetCDF file                                          !
      ntemp = nmod
      nproc = 1
      comm = 1


      numnu   = nnnnu
      cutoff  = 0.1d-03
      ifop(14)= .false.
      ifop(15)= .false. 
      ifop(16)= .false.
      ifop(17)= .false. 
      ifpres = .true. 
      iter   = 1  
      numits = 1
      glog = 0.0d0
! start loop over different p and T
! use a certain tsave paired with all nrhox values of p

     
      do i=1,  nmod



       nproc =i


       do j= 1, nrhox
         t(j) = tsave(i)
         tk(j)     = k * t(j)
         hckt(j)   = hc / tk(j)
         hkt(j)    = h / tk(j)
         tkev(j)   = k_ev * t(j)
         tlog(j)   = log(t(j))
         vturb(j)  = 0.0d0
         xnatom(j) = p(j) / tk(j) - xne(j)
         rho(j)    = xnatom(j) * wtmole * 1.660d-24
         pturb(j) = 0.5 * rho(j) * vturb(j) **2
 
       end do 

         nt = nrhox
         teff = t(1)

       itemp = itemp +1 


              call pops(0.0d0, 1 , xne)

!
               call pops( 1.01d0, 11, xnfph)
               call pops( 2.02d0, 11, xnfphe)
               call pops( 6.03d0, 11, xnfpc)
               call pops( 7.04d0, 11, xnfpn)
               call pops( 8.05d0, 11, xnfpo)
               call pops(10.05d0, 11, xnfpne)
               call pops(12.01d0, 11, xnfpmg)
               call pops(13.00d0, 11, xnfpal)
               call pops(14.01d0, 11, xnfpsi)
               call pops(20.01d0, 11, xnfpca)
               call pops(26.00d0, 11, xnfpfe)




!*
               if (.not. ifmol) then
                  call pops( 1.01d0, 12, xnfh)
                  call pops( 2.02d0, 12, xnfhe)
                  call pops( 6.05d0, 12, xnfc)
                  call pops( 7.05d0, 12, xnfn)
                  call pops( 8.05d0, 12, xnfo)
                  call pops(10.05d0, 12, xnfne)
                  call pops(12.05d0, 12, xnfmg)
                  call pops(14.05d0, 12, xnfsi)
                  call pops(16.05d0, 12, xnfs)
                  call pops(26.04d0, 12, xnffe)

!*
               else
                  call pops(106.00d0, 11, xnfpch)
                  call pops(108.00d0, 11, xnfpoh)
!*
!*...  THE POPS WILL NOT RETURN NUMBER DENSITIES WHEN MOLECULES ARE ON
!*...  SO WE COMPUTE NUMBER DENSITIES/PART FUNCTIONS  AND PART FUNCTIONS
!*
                  call pops( 6.05d0, 11, xnfc)
                  call pops( 7.05d0, 11, xnfn)
                  call pops( 8.05d0, 11, xnfo)
                  call pops(10.05d0, 11, xnfne)
                  call pops(12.05d0, 11, xnfmg)
                  call pops(14.05d0, 11, xnfsi)
                  call pops(16.05d0, 11, xnfs)
                  call pops(26.04d0, 11, xnffe)
!*
                  do j = 1,nrhox
!*
                     call pfsaha(j, 1, 1, 3, part)
                     xnfh(j, 1) = xnfph(j, 1) * part(j, 1)
                     xnfh(j, 2) = xnfph(j, 2)
!*
                     call pfsaha(j, 2, 2, 13, part)
                     xnfhe(j, 1) = xnfphe(j, 1) * part(j, 1)
                     xnfhe(j, 2) = xnfphe(j, 2) * part(j, 2)
                     xnfhe(j, 3) = xnfphe(j, 3)
!*
                     call pfsaha(j, 6, 6, 13, part)
                     xnfc(j, 1) = xnfc(j, 1) * part(j, 1)
                     xnfc(j, 2) = xnfc(j, 2) * part(j, 2)
                     xnfc(j, 3) = xnfc(j, 3) * part(j, 3)
                     xnfc(j, 4) = xnfc(j, 4) * part(j, 4)
                     xnfc(j, 5) = xnfc(j, 5) * part(j, 5)
                     xnfc(j, 6) = xnfc(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 7, 6, 13, part)
                     xnfn(j, 1) = xnfn(j, 1) * part(j, 1)
                     xnfn(j, 2) = xnfn(j, 2) * part(j, 2)
                     xnfn(j, 3) = xnfn(j, 3) * part(j, 3)
                     xnfn(j, 4) = xnfn(j, 4) * part(j, 4)
                     xnfn(j, 5) = xnfn(j, 5) * part(j, 5)
                     xnfn(j, 6) = xnfn(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 8, 6, 13, part)
                     xnfo(j, 1) = xnfo(j, 1) * part(j, 1)
                     xnfo(j, 2) = xnfo(j, 2) * part(j, 2)
                     xnfo(j, 3) = xnfo(j, 3) * part(j, 3)
                     xnfo(j, 4) = xnfo(j, 4) * part(j, 4)
                     xnfo(j, 5) = xnfo(j, 5) * part(j, 5)
                     xnfo(j, 6) = xnfo(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 10, 6, 13, part)
                     xnfne(j, 1) = xnfne(j, 1) * part(j, 1)
                     xnfne(j, 2) = xnfne(j, 2) * part(j, 2)
                     xnfne(j, 3) = xnfne(j, 3) * part(j, 3)
                     xnfne(j, 4) = xnfne(j, 4) * part(j, 4)
                     xnfne(j, 5) = xnfne(j, 5) * part(j, 5)
                     xnfne(j, 6) = xnfne(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 12, 6, 13, part)
                     xnfmg(j, 1) = xnfmg(j, 1) * part(j, 1)
                     xnfmg(j, 2) = xnfmg(j, 2) * part(j, 2)
                     xnfmg(j, 3) = xnfmg(j, 3) * part(j, 3)
                     xnfmg(j, 4) = xnfmg(j, 4) * part(j, 4)
                     xnfmg(j, 5) = xnfmg(j, 5) * part(j, 5)
                     xnfmg(j, 6) = xnfmg(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 14, 6, 13, part)
                     xnfsi(j, 1) = xnfsi(j, 1) * part(j, 1)
                     xnfsi(j, 2) = xnfsi(j, 2) * part(j, 2)
                     xnfsi(j, 3) = xnfsi(j, 3) * part(j, 3)
                     xnfsi(j, 4) = xnfsi(j, 4) * part(j, 4)
                     xnfsi(j, 5) = xnfsi(j, 5) * part(j, 5)
                     xnfsi(j, 6) = xnfsi(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 16, 6, 13, part)
                     xnfs(j, 1) = xnfs(j, 1) * part(j, 1)
                     xnfs(j, 2) = xnfs(j, 2) * part(j, 2)
                     xnfs(j, 3) = xnfs(j, 3) * part(j, 3)
                     xnfs(j, 4) = xnfs(j, 4) * part(j, 4)
                     xnfs(j, 5) = xnfs(j, 5) * part(j, 5)
                     xnfs(j, 6) = xnfs(j, 6) * part(j, 6)
!*
                     call pfsaha(j, 26, 5, 13, part)
                     xnffe(j, 1) = xnffe(j, 1) * part(j, 1)
                     xnffe(j, 2) = xnffe(j, 2) * part(j, 2)
                     xnffe(j, 3) = xnffe(j, 3) * part(j, 3)
                     xnffe(j, 4) = xnffe(j, 4) * part(j, 4)
                     xnffe(j, 5) = xnffe(j, 5) * part(j, 5)
                  end do
!*
               end if








       if (otherselection) then 

         call pops(1.01d0, 12, xnfh)   ! calculte 1.01 as ionisation fraction is needed for H f.f. transitions
         call pops(2.02d0, 12, xnfhe)

         call pops( 1.01d0, 11, xnfph)
         call pops( 2.02d0, 11, xnfphe)
         call pops( 5.00d0, 11, xnfpb) 
         call pops( 6.01d0, 11, xnfpc)    
         call pops( 7.00d0, 11, xnfpn)
         call pops( 8.00d0, 11, xnfpo)  
         call pops(11.00d0, 11, xnfpna)
         call pops(12.01d0, 11, xnfpmg)
         call pops(13.01d0, 11, xnfpal)
         call pops(14.01d0, 11, xnfpsi)
         call pops(19.00d0, 11, xnfpk) 
         call pops(20.01d0, 11, xnfpca)
         call pops(26.00d0, 11, xnfpfe)



        IF(IFMOL)CALL POPS(106.00D0,11,XNFPCH)
        IF(IFMOL)CALL POPS(108.00D0,11,XNFPOH)
 

! get the right number densities for the output file!
       if (ifmol) then
         do j = 1,nrhox

           call pfsaha(j, 1, 1, 3, part)
           xnfh(j, 1) = xnfph(j, 1) * part(j, 1)
           xnfh(j, 2) = xnfph(j, 2)

           call pfsaha(j, 2, 2, 13, part)
           xnfhe(j, 1) = xnfphe(j, 1) * part(j, 1)
           xnfhe(j, 2) = xnfphe(j, 2) * part(j, 2)
           xnfhe(j, 3) = xnfphe(j, 3)
         end do
       end  if

        do j = 1, nrhox

         xnfh2(j) = 0.0d0

         if (t(j) .le. 9000) then
          eq=exp(4.478/tkev(j)-4.64584d1+(1.63660d-3+(-4.93992d-7 &
     &       +(1.11822d-10+(-1.49567d-14+(1.06206d-18- &
     &        3.08720D-23*t(j))*t(j))*t(j))*t(j))*t(j))*t(j) &
     &       -1.5*tlog(j))

           xnfh2(j) = xnfh(j,1)**2*eq
           xnfph2(j) = xnfph(j,1)**2*eq
           xnfpco(j) = xnfpc(j,1)*xnfpo(j,1)*exp(11.091/tkev(j)-49.0414+ &
     &            14.0306d-4*t(j)-26.6341d-8*t(j)**2+35.382d-12*t(j)**3- &
     &            26.5424d-16*t(j)**4+8.32385d-20*t(j)**5-1.5*tlog(j))
         end if

        end do



      end if 



        do nu=1,  numnu

         freq=freqset(nu)
         freq15 = freq/1.0d15
         rco = 0.0d0
         frqlg = log(freq)
         freqln = frqlg
         freqlg = log10(freq)
 
         do j = 1, nrhox
           ehvkt(j) = exp(-freq*hkt(j))
           stim(j) = 1.0d0 - ehvkt(j)
           bnu(j) =1.47439d-02*freq15**3*ehvkt(j)/stim(j)  
           dbnudt(j) = bnu(j) * freq * hkt(j) / t(j) / &
     &                           stim(j)
          if(numnu .eq. 1) dbnudt(j) = 4.0d0 * sigma / pi * &
     &                                 t(j) ** 3




         end do 

         n=1
         nsteps = 1
         stepwt = 1.0d0 

         wave =c_nm/freq 
         waveno= 1.0d7/wave 
         call kapp 


         do j = 1, nrhox
           abtot(j) = (acont(j)+sigmac(j))
           continall(nu,j) =(abtot(j))
           contabs(nu,j) = (((acont(j))))
           contscat(nu, j) =( sigmac(j))
           ablog(j) = log10(abtot(j))

         end do


       ! use the background formula to calculate the other opacity
!        fit(1) =  0.129854
!        fit(2) =  2.46688
!        fit(3) =  1.48724
!        fit(4) = -0.191331
!        fit(5) =  0.00556552 
!
!        xo = 1.0d5 * wave * 1.0d-7
!        yo = fit(1) + xo*fit(3) +(xo*(fit(3) +xo*(fit(4) + xo*fit(5))))
!        yo = max(yo,0.0) 
!
!
!        xo = 1.0d8*wave*1.0d-7
!        do j = 1, nrhox
!         xnatom(j) = p(j) / tk(j) - xne(j)
!         rho(j)    = xnatom(j) * wtmole * 1.660d-24
!          opaca(j) = 4.58d-44 * xo * (1.0d0 +xo /941.)*xne(j)* xnfh(j,1) 
!          opaca(j) = opaca(j) + 1.04d-34 *xnfh(j,1) *xne(j) /t(j)/sqrt(t(j) )* exp(8761./t(j))*yo 
!          opaca(j) = opaca(j)
!          abtot(j) = acont(j) *rho(j)  
!        end do 


        end do !end of nu loop

        
       ! use the background formula to calculate the other opacity


       

!--------------------------------------------------------------------

           do  nu=1,1131
              continall4(nu)=continall(nu,j)
           end do
 
         
        do j=1, nrhox           
           do nu=1, 328
              contgrid (nu, j, i, 1) = continall((nu-1)*3 +2,j ) 
              contgrid (nu, j, i, 2) = contabs((nu-1)*3 +2, j ) 
             contgrid (nu , j, i, 3) = contscat((nu-1)*3 +2,j ) 
           end do 
        end do  
 


!      print*, 'this is nmod = ', i
      end do



      do nu = 1, 328
        do i = 1, nmod
          do j = 1, nrhox
            write(77,*) contgrid(nu,j,i,1) ! absolute
            write(88,*) contgrid(nu,j,i,2) ! absorption only 
            write(79,*) contgrid(nu,j,i,3) ! scattering only 
          end do 
        end do 
      end do 


      deallocate(contgrid)


     ! close all files
      close (unit =44) ! close conintinua.dat file
      close (unit =33) ! model file
      close (unit =15) ! INPUT model
      close (unit =16) !, status = 'delete')
      close (unit =17)!, status = 'delete')
      close (unit =18) !, status = 'delete')
      close (unit =77)
      close (unit =88)
      close (unit =79)
 


 end subroutine

