  subroutine calc_odf 


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

! Note, the common blockes that are commented are in one of the modules above


!
!--------------------------- LOCAL VARIABLES --------------------------
!
      real(kind=dp)  contin, freq15, part(maxd, 6), rco  
      real(kind=dp)  rcowt, stepwt, sumwt, wave, x, rosstab, asixth
      integer  i, j, mm,  mode1, n, nsteps, nu, nmod, imod
      logical  finish, more, readi0, stopfl
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

      integer  last1, nnnnu, nelem, ion, nu3, iT, iP 
      real(kind=dp) save, edge,  eq, cutoff, frqlg, contmin 
      real(kind=dp) freeff
      parameter (max1=maxmol+1)
      
      logical nosaha, checkinf, isinf 

      real(kind=dp) xnfh2(maxd), xnfph2(maxd), xnfpco(maxd)
      real(kind=dp) xnfp(maxd,10,99),xnfpel(6,99),dopple(6,99)
      real(kind=dp) dummyxn(maxd,6)
      real(kind=dp) ablog(maxd)
      real(kind=dp) a(377), continall(1131, maxd), frqedg(377) 
      real(kind=dp) contabs(1131,maxd)
      real(kind=dp) contscat(1131,maxd)
!
      real(kind=dp)  xnfh4(maxd), xnfhe4(maxd,2)
      real(kind=dp)  xnfpel4(6,99,maxd), dopple4(6,99, maxd), continall4(1131)
      real(kind=dp)  xnfdopmax(377,6,99), xnfdopcont 

!---- for netcdf files ----------------------------------------------!
       integer   ncid, ier, myrank, ntemp, nproc, nufreq, comm 
       character(80) :: name
!---------------------------------------------------------------------
!      common /xnmol/codemol(maxmol), xnfpmol(maxd,maxmol),nummol

      character card*80  
      real(kind=dp)  WAVEBIG(329)
      real(kind=dp)  BIGA(95),BIGB(125),BIGC(993)
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
     &      access = 'sequential')
      open (unit = 17, file = 'mpsa.punch', form = 'formatted', &
     &      access = 'sequential')
      open (unit = 18, file = 'mpsa.jnu', form = 'formatted', &
     &      access = 'sequential')
      open (unit = 33, file = filePTgrid, form = 'formatted', &
     &      status = 'old', access = 'sequential')

!---------------------------------------------------------------------
     
      open (unit = 44, file = './INPUT/continua.dat', form = 'formatted', &
     &      status = 'old', access = 'sequential') 

      more = .true.
      nosaha = .false.
      
!
      maxpow = 99
      last=80
      numcol=1
      in = 0


      place = 0  
      i = 0 
      do while (more)
! read file 44: continua.dat the final line has to have 'begin'! 
        i = i +1
        numcol=1
        if (place .eq. 0) then 
         read(44,100) card
100      format(a)
!         print*, card
        end if
        if (i .gt. 1000) then 
          more = .false.
        end if
       if (card .eq. 'begin' ) then 
        more = .false.
       else
        edge = freeff(card, place) 
        if (edge .eq. 0 )  then 
         place = 0  
        else 

        in= in+1
!         print *, in, edge
  
         if (abs(edge) .lt. 1.0d6) then
            wledge(in) = edge
            cmedge(in) = 1.0d7/edge
            frqedg(in) = c_nm/wledge(in)

         else if (abs(edge) .lt. 1.0d25) then 
            frqedg(in) = edge
            wledge(in) = c_nm/edge
            cmedge(in) = 1.0d7/wledge(in) 


         else
            cmedge(in) = edge/1.0d25
            wledge(in) = 1.0d7/cmedge(in)
            frqedg(in) = c_nm/wledge(in)

         end if
         a(in) = abs(wledge(in))
 
        end if 
       end if
  
      end do
 
 

      
!  sort a, in ascending order
!
      do last= 2, in
       last1= in -last +2

         do i= 2, last1
          if (a(i) .lt. a(i-1)) then
          
            save  = a(i-1)
            a(i-1)= a(i) 
            a(i)  = save 

            save  = frqedg(i-1)
            frqedg(i-1) = frqedg(i)
            frqedg(i)   = save

            save        = wledge(i-1)
            wledge(i-1) = wledge(i)
            wledge(i)   = save
           
            save        = cmedge(i-1)
            cmedge(i-1) = cmedge(i)
            cmedge(i)   = save
 
          end if
         end do

      end do

      numnu=0
      do i=1, in-1
       numnu= numnu +1
       freqset( numnu ) = abs(frqedg(i))/1.0000001
       numnu= numnu +1
       freqset( numnu ) = c_nm/(abs(wledge(i)) + abs(wledge(i+1)))*2.0d0
       numnu = numnu +1
       freqset( numnu ) = abs(frqedg(i+1))*1.0000001
 
      end do 
    
!      write(6,15) numnu
      nnnnu= numnu

!---------- only now read the model
      more = .true.
!-------------------------------------------------------------------
! call reading model (reads header and how many atmosphere models
      mode1 = 0
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
!     READ MODEL (first nmod T, then nrhox p values)
!     Because we need one atmosphere at a time keep T fixed
!     and calculate for all p values -> for not pretabulated 
!     calculations we need to adjust here!!! 
!-------------------------------------------------------!
            
               read (33, *) (tsave(j), j = 1, nmod)
               read (33, *) (p(j), j=1, nrhox)

         call print_summary(0) 

!-------------------------------------------------------!
! now we now dimensions of the calculations and can creat 
!  NetCDF file                                          !
      ntemp = nmod
      nproc = 1
      comm = 1

      call CreateContinum(myrank, ncid, 'Continum.nc', nrhox, ntemp, & 
     &                   nufreq, ier )
      name = 'Continum.nc'


      numnu   = nnnnu
      cutoff  = 1.0d-03
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


        if (nosaha) then 


! --- Calculate Saha for certain elements ----!
! carefull xnfh was a xnfh(nrohx) array but in the
! new atlas it is xnfh(nrhox,2) --> what is the 2nd index --> it is the first and second 
! and nth ion!!!

       else 
         call pops(0.0d0, 1 , xne)

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



        do nu=1, numnu
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
           continall(nu,j) =log10(abtot(j))
           contabs(nu,j) = (log10((acont(j))))
           contscat(nu, j) =log10( sigmac(j))
           ablog(j) = log10(abtot(j))

         end do


        end do !end of nu loop





        do j=1,nrhox
         xnfh4(j)=xnfh(j,1)
         xnfhe4(j,1)=xnfhe(j,1)
         xnfhe4(j,2)=xnfhe(j,2)
        end do


        do nelem = 1, 99
          do ion = 1, 10
            do  j=1, nrhox
              xnfp(j,ion,nelem)=0.0d0
            end do
          end do
        end do 


        call pops(1.01d0,11,xnfp(1,1,1))
        call pops(2.02d0,11,xnfp(1,1,2))
        call pops(3.03d0,11,xnfp(1,1,3))
        call pops(4.03d0,11,xnfp(1,1,4))
        call pops(5.03d0,11,xnfp(1,1,5))
        call pops(6.05d0,11,xnfp(1,1,6))
        call pops(7.05d0,11,xnfp(1,1,7))
        call pops(8.05d0,11,xnfp(1,1,8))
        call pops(9.05d0,11,xnfp(1,1,9))
        call pops(10.05d0,11,xnfp(1,1,10))
        call pops(11.05d0,11,xnfp(1,1,11))
        call pops(12.05d0,11,xnfp(1,1,12))
        call pops(13.05d0,11,xnfp(1,1,13))
        call pops(14.05d0,11,xnfp(1,1,14))
        call pops(15.05d0,11,xnfp(1,1,15))
        call pops(16.05d0,11,xnfp(1,1,16))
        call pops(17.04d0,11,xnfp(1,1,17))
        call pops(18.04d0,11,xnfp(1,1,18))
        call pops(19.04d0,11,xnfp(1,1,19))

        call pops(20.05d0,11,xnfp(1,1,20))

        call pops(21.05d0,11,xnfp(1,1,21))    !Ionsation stages to 5, before it was 9)

        call pops(22.05d0,11,xnfp(1,1,22))
        call pops(23.05d0,11,xnfp(1,1,23))
        call pops(24.05d0,11,xnfp(1,1,24))
        call pops(25.05d0,11,xnfp(1,1,25))
        call pops(26.05d0,11,xnfp(1,1,26))
        call pops(27.05d0,11,xnfp(1,1,27))
        call pops(28.05d0,11,xnfp(1,1,28))
 

        do  nelem=29,99
          call pops(float(nelem)+.02d0,11,xnfp(1,1,nelem))
        end do

         
        do  j=1,nrhox
          xnfp(j,10,40)=xnfph2(j)
          xnfp(j,10,46)=xnfpco(j)
        end do

        if(ifmol) then 
          do nelem=40,99
             call pops(idmol(nelem-39),1,xnfp(1,6,nelem))
          end do
        end if

        do  j=1,nrhox
          do  nelem=20,28
            xnfp(j,5,30+nelem)=xnfp(j,7,nelem)
            xnfp(j,5,40+nelem)=xnfp(j,8,nelem)
            xnfp(j,5,50+nelem)=xnfp(j,9,nelem)
            xnfp(j,5,60+nelem)=xnfp(j,10,nelem)
          end do
        end do



        do  nelem=1,99
          do  ion=1,6
           nu3=0
           do  nu=1,1131,3
            nu3=nu3+1
            xnfdopmax(nu3,ion,nelem)=0.
           end do
          end do
        end do 


!-----------------------------------------------      
!     large nrhox 300 loop
!----------------------------------------------
        do  j=1,nrhox

           do  nelem=1,99
             do  ion=1,6
               xnfpel(ion,nelem)=xnfp(j,ion,nelem)
             end do
           end do




           do  nelem=1,99
             dopple(1,nelem)=sqrt(2.*tk(j)/atmass(nelem)/1.660d-24+ &
     &                      vturb(j)**2)/2.99792458d10
             dopple(2,nelem)=dopple(1,nelem)
             dopple(3,nelem)=dopple(1,nelem)
             dopple(4,nelem)=dopple(1,nelem)
             dopple(5,nelem)=dopple(1,nelem)
             dopple(6,nelem)=dopple(1,nelem)
           end do


           do  nelem=20,28
             dopple(5,30+nelem)=dopple(1,nelem)
             dopple(5,40+nelem)=dopple(1,nelem)
             dopple(5,50+nelem)=dopple(1,nelem)
             dopple(5,60+nelem)=dopple(1,nelem)
           end do

           if(ifmol) then
             do  nelem=40,99
               dopple(6,nelem)=sqrt(2.*tk(j)/momass(nelem-39)/1.660d-24) &
     &                         /2.99792458d10
             end do 
           end if


           do  nelem=1,99
            do  ion=1,6
              xnfpel4(ion,nelem, j)=xnfpel(ion,nelem)
              dopple4(ion,nelem,j)=dopple(ion,nelem)

            end do
           end do 

           nu3=0
           do  nu=1,numnu,3
              nu3=nu3+1

           if (j .eq. 15) then
           end if

              contmin=min(continall(nu,j),continall(nu+1,j), &
     &                 continall(nu+2,j))
              contmin=10.**contmin
              do nelem=1,99
                do  ion=1,6
                   xnfdopcont=xnfpel(ion,nelem)/dopple(ion,nelem)/ &
     &                         freqset(nu+1)/rho(j)/contmin/cutoff
              xnfdopmax(nu3,ion,nelem)=max(xnfdopmax(nu3,ion,nelem) &
     &                                       ,xnfdopcont)
             end do
            end do
          if (j .eq. 15) then
          endif

           end do



          it=itemp
          ip=j
!--------------------------------------------------------------------

           do  nu=1,1131
              continall4(nu)=continall(nu,j)
           end do
 

        end do

!----- end of large nrhox 300 loop
!-------------------------------------------------------------------------
       call WriteContinum( ncid, myrank, nproc, comm, &
     &      wledge, frqedg, cmedge, idmol, momass, freqset, &
     &      tsave, p, t, rho, xne, xnatom, vturb, xnfh2, xnfhe4, &
     &      xnfh4, continall, contscat, xnfpel4, xnfdopmax, dopple4,&
     &      nrhox, ntemp, nufreq , ier)


      end do
! clean up and close stuff
     call CloseNetCDF (myrank, ncid, ier)
!--- make sure that this is deallocated again:
      call close_freq

     ! close all files
      close (unit =44) ! close conintinua.dat file
      close (unit =33) ! model file
      close (unit =15) ! INPUT model
      close (unit =16, status = 'delete')
      close (unit =17, status = 'delete')
      close (unit =18, status = 'delete')
      close (unit =77, status= 'delete')

 


 end subroutine



!---- the subroutine below was for checking the output produced by the old f70 version of the DFSYNTHE code (can be ignored)
  subroutine cont_rw 
  


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
      double precision  contin, freq15, part(maxd, 6), rco, & 
     &                  rcowt, stepwt, sumwt, wave, x, rosstab, asixth
!     double precision exp10
      integer  i, j, mm,  mode1, n, nsteps, nu, nmod, imod
      logical  finish, more, readi0, stopfl
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
      
      logical nosaha, checkinf, isinf, header, damien 

      double precision xnfh2(maxd), xnfph2(maxd), xnfpco(maxd)
      real(kind=dp) xnfp(maxd,10,99)
      real(kind=4) xnfpel(6,99),dopple(6,99)
      double precision dummyxn(maxd,6)
      double precision ablog(maxd)
! in common:  double precision  wledge(377), cmedge(377), freqedge(377),  freqset(1131) 

      double precision a(377), continall(1131, maxd), frqedg(377) 
      double precision CONTABS(1131,maxd)
      double precision CONTSCAT(1131,maxd)
      real(kind=4)  t4(25),tkev4(25), tk4(25),hkt4(25),tlog4(25)
      real(kind=4)  hckt4(25) 
      
      real(kind=4)  p4(25), xne4(25), xnatom4(25),rho4(25),rhox4(25)
      real(kind=4)  VTURB4(25),XNFH4(25),XNFHE4(25,2),XNFH24(25)


      real(kind=dp)  t8(25),tkev8(25), tk8(25),hkt8(25),tlog8(25)
      real(kind=dp)  hckt8(25)

      real(kind=dp)  p8(25), xne8(25), xnatom8(25),rho8(25),rhox8(25)
      real(kind=dp)  VTURB8(25),XNFH8(25),XNFHE8(25,2),XNFH28(25)


      real(kind=dp) XNFPEL4(6,99,maxd),DOPPLE4(6,99, maxd)
      real(kind=4) CONTINALL4(1131)
      real(kind=4) XNFDOPMAX(377,6,99),XNFDOPCONT
      real(kind=dp) XNFDOPMAX8(377,6,99) 
!---- for netcdf files ----------------------------------------------!
       integer   ncid, ier, myrank, ntemp, nproc, nufreq, comm 
       character(80) :: name
!---------------------------------------------------------------------
!      common /xnmol/codemol(maxmol), xnfpmol(maxd,maxmol),nummol
       character(len=80) card  
!      dimension card(81)
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







      myrank = 0
      nufreq = 1131
 

      idmol= (/ 101.,  106.,  107.,  108.,  606.,  607.,  608., & 
     & 707.,  708., & 
     &   808.,  112.,  113.,  114.,  812.,  813.,  814.,  116.,  120., & 
     &   816.,  820.,  821.,  822.,  823.,  103.,  104.,  105.,  109., & 
     &   115.,  117.,  121.,  122.,  123.,  124.,  125.,  126.,106.01, &
     & 107.01,108.01,112.01,113.01,114.01,120.01,  408.,  508.,  815., &
     &   817.,  824.,  825.,  826.,10108.,60808.,10106.,60606.,  127., &
     &   128.,  129.,  827.,  828.,  829., 608.01 /) 
 
      momass=(/ 2.,   13.,   15.,   17.,   24.,   26.,   28.,   28., & 
     &   30., & 
     &    32.,   25.,   28.,   29.,   40.,   43.,   44.,   33.,   41., &
     &    48.,   56.,   61.,   64.,   67.,    8.,   10.,   12.,   20.,& 
     &    32.,   36.,   46.,   49.,   52.,   53.,   56.,   57.,   13., &
     &    15.,   17.,   25.,   28.,   29.,   41.,   25.,   27.,   47., &
     &    51.,   68.,   71.,   72.,   18.,   44.,   14.,   36.,   60., &
     &    59.,   64.,   75.,   74.,   79.,   28. /)

! READ in input header first, as later on you want to be able to allocate 
! the right size array for small and big wl bins (not subbins)!

 nproc = 1
 nmod = 57
 nrhox = 25





!-------------------------------------------------------!
! now we now dimensions of the calculations and can creat 
!  NetCDF file                                          !
      ntemp = nmod
      nproc = 1
      comm = 1
      call CreateContinum( myrank, ncid, 'Continum_orig.nc', nrhox, ntemp, & 
     &                   nufreq, ier )
      name = 'Continum_orig.nc'

! --- open the files!

     
      do i=1,  nmod


       nproc =i 

       read(10)nt,teff,glog,title
       read(11)nt,teff,glog,title
       read(10)in,(frqedg(j),wledge(j),cmedge(j),j=1,in),idmol,momass
       read(11)in,(frqedg(j),wledge(j),cmedge(j),j=1,in)
       read(10)numnu,(freqset(nu),nu=1,numnu)
       read(11)numnu,(freqset(nu),nu=1,numnu)
 
       write(20,*)nt,teff,glog,title
       write(21,*)nt,teff,glog,title
       write(20,*)in,(frqedg(j),wledge(j),cmedge(j),j=1,in),idmol,momass
       write(21,*)in,(frqedg(j),wledge(j),cmedge(j),j=1,in)
       write(20,*)numnu,(freqset(nu),nu=1,numnu)
       write(21,*)numnu,(freqset(nu),nu=1,numnu)



        read(10)t4,tkev4,tk4,hkt4,tlog4,hckt4,p4,xne4,xnatom4,rho4, & 
     &    rhox4,vturb4,xnfh4,xnfhe4,xnfh24

        
!        write(20,*)t4,tkev4,tk4,hkt4,tlog4,hckt4,p4,xne4,xnatom4,rho4, & 
!     &    rhox4,vturb4,xnfh4,xnfhe4,xnfh24



        do j = 1,nrhox


          t8(j)       = dble(t4(j))
          p8(j)       = dble(p4(j))
          xne8(j)     = dble(xne4(j))
          xnatom8(j)  = dble(xnatom4(j))
          rho8(j)     = dble(rho4(j))
          VTURB8(j)   = dble(vturb4(j)) 
          XNFHE8(j,1) = dble(xnfhe4(j,1))
          XNFHE8(j,2) = dble(xnfhe4(j,2))
          XNFH28(j)   = dble(xnfh24(j)) 
          xnfh8(j)    = dble(xnfh4(j))

!--
!        read(10) string
        read(10) xnfpel, dopple, continall4

!--------------------------------------------------------------------
          do  nelem=1,99
            do  ion=1,6
              xnfpel4(ion,nelem, j)=dble(xnfpel(ion,nelem))
              dopple4(ion,nelem,j)=dble(dopple(ion,nelem))

            end do
           end do

           do  nu=1,1131
              continall(nu,j) = dble(continall4(nu)) 
           end do
 
        WRITE(20,*)'XNFPEL,dopple, CONTINALL4'
        WRITE(20,*) xnfpel4(:,:,j), dopple4(:,:,j), CONTINALL(:,j) 




        end do

!----- end of large nrhox 300 loop
!-------------------------------------------------------------------------


        write(20,*)t8,tkev4,tk4,hkt4,tlog4,hckt4,p8,xne8,xnatom8,rho8, &
     &    rhox4,vturb8,xnfh8,xnfhe8,xnfh28




         read(11) xnfdopmax

         xnfdopmax8 = dble(xnfdopmax)
         WRITE(21,*)XNFDOPMAX8 
       call WriteContinum( ncid, myrank, nproc, comm, &
     &      wledge, frqedg, cmedge, idmol, momass, freqset, &
     &      tsave, p8, t8, rho8, xne8, xnatom8, vturb8, xnfh28, xnfhe8, &
     &      xnfh8, continall, contscat, xnfpel4, xnfdopmax8, dopple4,&
     &      nrhox, ntemp, nufreq , ier)



      end do

     call CloseNetCDF (myrank, ncid, ier)



 end subroutine

