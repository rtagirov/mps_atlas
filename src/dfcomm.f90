MODULE DFCOMM

  use types
  implicit none
!      This module allocated arrays for the opacity pre-tabulation calcualtions!
!
!.... maxd   = maximum number of depths in the atmosphere
!.... maxmeq = maximum number of molecular equations
!.... maxloc = maximum number of molecular components
!.... maxmol = maximum number of molecules
!.... maxmu  = maximum number of angles for the radiation field
!.... maxnu  = the maximum number of frequencies


! --------------------------------------------

  integer, parameter :: maxd = 360 
  integer, parameter :: lenbuff = 3510000
  integer, parameter :: maxprof = 10000
!
!  for ODF's T and P gird

  real(kind=dp) :: profile(10000), b10000(10000)
  real(kind=dp) :: extab(10001), e1tab(2000)
  real(kind=dp) :: h0tab(2001), h1tab(2001), h2tab(2001), atab(281,1001)

  real(kind=dp) :: xnfpel(594),dopple(594),xnfdop(594) !6*99=594 -> 6 levels and 99 atoms 
  real(kind=dp) :: buffer(lenbuff), continuum(lenbuff)

 ! for  Kuruzc  bin-sizes
  real(kind=dp) :: biga(95),bigb(125),bigc(109)
  real(kind=dp) :: lita(95),litb(190),litc(190),litd(190),lite(190)
  real(kind=dp) :: litf(190),litg(168)

  real(kind=dp) :: wlbeg,wlend, wlbeg0, resolu,ratio,ratiolg,wbegin, wend
  integer :: ixwlbeg, length, mlines 
  integer :: nsizebig, nsizelit
  integer :: nufreq, nufreq3

! for dynamical binning
  integer :: nsubbin

!--- all allocatable arrays

  integer,  allocatable :: nbeg(:), nend(:)
  integer,  allocatable ::  nbigbeg(:), nbigend(:), nlitbeg(:), nlitend(:)
  integer,  allocatable :: nn(:), nsteps(:,:)
  integer(kind=2), allocatable :: insteps(:), iodfstep(:,:,:,:,:)
  integer(kind=2), allocatable :: iodfsendl(:,:,:,:), iodfrecvl(:,:,:,:)
  integer(kind=2), allocatable :: iodfsendb(:,:,:,:), iodfrecvb(:,:,:,:)
 
  integer, allocatable :: nvt(:)

  real(kind=dp), allocatable :: wavebig(:), wavelit(:)
  real(kind=dp), allocatable ::  emerge(:), txnxn(:), bstim(:), xnfh2(:)

  real(kind=dp), allocatable :: dopmax(:,:,:)

  real(kind=dp), allocatable :: continall(:,:), contscat(:,:) 
  real(kind=dp), allocatable :: dopple4(: ,: ,:), xnfpel4(:,:,:) 
  real(kind=dp), allocatable  :: vtprof(:,:),vtcenter(:)

! allocatable frequency grid
  real(kind=dp), allocatable  ::  inifreset(:)

! for dynamical binning
  real(kind=dp), allocatable :: binbound(:,:),sbwith(:,:)
!-----------------------------------------------------------------------------------
!clsa
!related to filters
      integer :: nnsize, Nlambda_phi, numssbin, nntot
      integer :: iphi, ic
      integer, allocatable :: nnsteps(:,:),indexx(:)
      real (kind=dp),allocatable :: wavbigtilda(:),wavlittilda(:)
      real (kind=dp), allocatable :: dlambda(:,:), dlamtilda(:,:), phi_hr(:,:)
      real (kind=dp), allocatable :: xlambda_phi(:),filter(:)
      real (kind=dp) :: phi_lp,phi_avg,dlam,lambda0,lambda,lambdatilda,lambdatildaus,totdlam
      real (kind=dp) :: xlambda_min, xlambda_max
!-------------------------------------------------------------------------------------------
      real (kind=dp), allocatable :: sb_weight_send(:,:), sb_weight_recv(:,:)
!
!-------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

CONTAINS

SUBROUTINE set_arrays(numsbin, numssbin,ntemp) 

  implicit none
  include 'common.rhoxbl'
  include 'common.turbpr'
  include 'common.ifopbl'
 
  integer, intent(in) ::ntemp, numsbin
  integer, intent(in) :: numssbin
  integer :: ierr(10), i

 do i = 1, 10 
  ierr(i)=-1
 end do 
 
 if (numsbin .eq. 0) then 
  nsizebig = 328 
  nsizelit = 1212 
 else
  nsizebig = numsbin-1 
  nsizelit = 1
 endif


 allocate(iodfsendb(numssbin,nsizebig,nrhox,ntemp), stat = ierr(9))
 allocate(iodfrecvb(numssbin,nsizebig,nrhox,ntemp), stat = ierr(10))

 allocate(insteps(numssbin), stat=ierr(3))
 allocate(iodfstep(numssbin,nsizebig+nsizelit,nrhox,ntemp,numvt), stat = ierr(4) )
 
 allocate(nbeg(nsizebig+nsizelit), stat = ierr(5))
 allocate(nend(nsizebig+nsizelit), stat = ierr(6))
 allocate(nn(nsizebig+nsizelit), stat = ierr(7))
 allocate(nsteps(numssbin+2,nsizebig+nsizelit), stat= ierr(8))

 allocate(nbigbeg(nsizebig), stat = ierr(9))
 allocate(nbigend(nsizebig), stat = ierr(10))

 if (maxval(ierr) .gt. 0) then
   print*, 'first ten arrays could not be allocated'
   stop
 end if


!clsa
 nnsize=nsizebig+nsizelit
 if (ifilter .eq. 1) then 
   allocate (indexx(lenbuff), stat=ierr(1))
   allocate (dlambda(nnsize,lenbuff), stat=ierr(2))
   allocate (dlamtilda(nnsize,lenbuff), stat= ierr(3))
   allocate (phi_hr(nnsize,lenbuff), stat=ierr(4))
   allocate (wavbigtilda(nsizebig+1), stat=ierr(5))
   allocate (wavlittilda(nsizelit+1), stat=ierr(6))
   allocate (nnsteps(numssbin+1,nnsize), stat=ierr(7))

   allocate(sb_weight_send(numssbin,nsizebig+nsizelit), stat = ierr(10))
   allocate(sb_weight_recv(numssbin,nsizebig+nsizelit), stat = ierr(2))

  if (maxval(ierr) .gt. 0) then
   print*, 'some of the arrays for filters could not be allocated'
   stop
 end if


 end if 

!clsa

 

 allocate(nlitbeg(nsizelit), stat = ierr(1))
 allocate(nlitend(nsizelit), stat = ierr(2))

 allocate(wavebig(nsizebig+1), stat = ierr(3))
 allocate(wavelit(nsizelit+1), stat = ierr(4))

 allocate(binbound(nsizebig+nsizelit, numssbin),stat = ierr(5)) 
 allocate(sbwith(numssbin, nsizebig+nsizelit), stat = ierr(6))

 allocate(emerge(nrhox), stat=ierr(7))
 allocate(txnxn(nrhox), stat=ierr(8))
 allocate(bstim(nrhox), stat=ierr(9))

 allocate(xnfh2(nrhox), stat=ierr(10))

 if (maxval(ierr) .gt. 0) then
   print*, 'second ten arrays could not be allocated'
   stop
 end if 

! allocate(tsave(nrhox))
 allocate(dopmax(nufreq3, 6, 99), stat=ierr(1))
 allocate(continall(nufreq, nrhox), stat=ierr(2))
 allocate(contscat(nufreq, nrhox), stat=ierr(3))
 allocate(dopple4(6,99, nrhox), stat=ierr(4))
 allocate(xnfpel4(6,99, nrhox), stat=ierr(5)) 


 allocate(vtprof(120,numvt), stat=ierr(6))
 allocate(vtcenter(numvt), stat=ierr(7))
 allocate(nvt(numvt),stat=ierr(8))


 allocate(iodfsendl(numssbin,nsizelit,nrhox,ntemp), stat = ierr(9))
 allocate(iodfrecvl(numssbin,nsizelit,nrhox,ntemp), stat = ierr(10))

 if (maxval(ierr) .gt. 0) then
   print*, 'second ten arrays could not be allocated'
   stop
 end if

 nbeg=0
 nend=0
 nn =0
 nsteps=0
 nbigbeg=0
 nbigend =0
 nlitbeg= 0
 nlitend= 0
 
 wavebig=0.0d0
 wavelit=0.0d0 

 emerge=0.0d0
 txnxn = 0.0d0
 bstim = 0.0d0
 xnfh2 = 0.0d0
 dopmax = 0.0d0
 continall = 0.0d0
 contscat = 0.0d0
 dopple4 = 0.0d0
 xnfpel4 = 0.0d0

 iodfstep = 0
 iodfsendl = 0
 iodfrecvl = 0 
 iodfsendb = 0
 iodfrecvb = 0 

!lsa
 if (ifilter .eq. 1 ) then 
  sb_weight_send=0.
  sb_weight_recv=0.
 end if 

END SUBROUTINE set_arrays

SUBROUTINE close_arrays

  implicit none
  include 'common.ifopbl'

 deallocate(nbeg)
 deallocate(nend)
 deallocate(nn)
 deallocate(nsteps)

 deallocate(nbigbeg)
 deallocate(nbigend)
 deallocate(nlitbeg)
 deallocate(nlitend)

 deallocate(insteps)
 deallocate(iodfstep)
 deallocate(iodfsendb)
 deallocate(iodfrecvb)

 deallocate(iodfsendl)
 deallocate(iodfrecvl)

 deallocate(wavebig)
 deallocate(wavelit)

 deallocate(binbound)


 deallocate(emerge)
 deallocate(txnxn)
 deallocate(bstim)
 deallocate(xnfh2)
! deallocate(tsave) 
 deallocate(dopmax)
 deallocate(continall)
 deallocate(contscat)
 deallocate(dopple4)
 deallocate(xnfpel4)

 deallocate(vtprof)
 deallocate(vtcenter)
 deallocate(nvt) 

!lsa
 if (ifilter .eq. 1) then 
   deallocate (indexx)
   deallocate (dlambda)
   deallocate (dlamtilda)
   deallocate (phi_hr)
   deallocate (wavbigtilda)
   deallocate (wavlittilda)
   deallocate (nnsteps)

   deallocate (sb_weight_send)
   deallocate (sb_weight_recv)
!lsa
 endif 
 

END SUBROUTINE close_arrays

subroutine alloc_wave(numsbin, numssbin)

  implicit none
  integer, intent(in) :: numsbin, numssbin

  include 'common.rhoxbl'
  include 'common.turbpr'

 if (numsbin .eq. 0) then
  nsizebig = 328
  nsizelit = 1212

 else
  nsizebig = numsbin-1
  nsizelit = 1
 endif

 allocate(wavebig(nsizebig+1))
 allocate(wavelit(nsizelit+1))

end subroutine alloc_wave

subroutine close_wave

  implicit none
  deallocate(wavebig)
  deallocate(wavelit)


end subroutine close_wave



SUBROUTINE def_binsize

   implicit none 
   include 'common.ifopbl'

   integer i

   if (ifkbin) then 

     biga=(/8.97666,     9.2,         9.5,         9.71730,     9.81590, &
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
     &    67.,         69.,         71.,         72.24011,    74./)

    bigb=(/ 76.,         78.,         80.,         82.,         84., &
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
     &   910.,   920., 930., 940., 950., 960.,  970., 980.,  990.,1000./)
     bigc=(/1025.,1050.,1075.,1100.,1125.,1150.,1175.,1200.,1225.,1250., &
     & 1275.,1300.,1325.,1350.,1375.,1400.,1425.,1458.81670,1475.,1500., &
     & 1525.,1550.,1575.,1600.,1640.,1680.,1720.,1760.,1800.,1840., &
     & 1880.,1920.,1960.,2000.,2050.,2100.,2150.,2200.,2250.,2279.40330, &
     & 2300.,2350.,2400.,2450.,2500.,2550.,2600.,2650.,2700.,2750., &
     & 2800.,2850.,2900.,2950.,3000.,3050.,3100.,3150.,3200.,3282.34320, &
     & 3400.,3500.,3600.,3700.,3800.,3900.,4000.,4100.,4200.,4300., &
     & 4400.,4500.,4600.,4700.,4800.,4900.,5000.,5100.,5200.,5300., &
     & 5400.,5500.,5600.,5700.,5800.,5900.,6000.,6100.,6200.,6300., &
     & 6400.,6600.,6800.,7000.,7200.,7400.,7600.,7800.,8000.,8200., &
     & 8400.,8600.,8800.,9000.,9200.,9400.,9600.,9800.,10000./)

     lita=(/8.97666,     9.2,         9.5,         9.71730,     9.81590, &
     &    10.10901,    10.3,        10.46451,    10.65000,    10.88545, &
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
     &    50.42590,    50.87630,    51.5,        52.5,        53.5, &
     &    54.5,        55.5,        56.5,        57.40614,    58.5, &
     &    59.5,        60.5,        61.5,        62.5,        63.5/)
    litb=(/ 64.5, 65.53786, 66., 67., 68., 69., 70., 71., 72.24011,  73., &
     &  74., 75., 76., 77., 78., 79., 80., 81., 82., 83., &
     &  84., 85., 86., 87., 88., 89., 90., 91.17535, 92., 93., &
     &  94., 95., 96., 97., 98., 99., 100., 101., 102., 103., &
     & 104., 105., 106., 107., 108., 109., 110.03056, 111., 112., 113., &
     & 114., 115., 116., 117., 118., 119., 120., 121., 122., 123., &
     & 123.92928, 125., 126., 127., 128., 129., 130., 131., 132., 133., &
     & 134., 135., 136., 137., 138., 139., 140., 141., 142., 143., &
     & 144.,144.43391,145.,146.,147.,148.,149.,150.,151.,151.43485, &
     & 152., 153., 154., 155., 156., 157., 158., 159., 160., 161., &
     & 162.15072,163.,164.,165.,166.,167.,167.40282,168.,169.,170., &
     & 171., 172., 173., 174., 175., 176., 177., 178., 179., 180., &
     & 181., 182., 183., 184., 185., 186., 187., 188., 189., 190., &
     & 191.,192.,193.,194.,195.,196.,197.46990,197.81149,199.,200., &
     & 201.,202.,203.,204.,205.,206.,207.13210,207.61400,208.,209., &
     & 210., 211., 212., 213., 214., 215., 216., 217., 218., 219., &
     & 220., 221., 222., 223., 224., 225., 226., 227., 228., 229., &
     & 230., 231., 232., 233., 234., 235., 236., 237., 238., 239., &
     & 240., 241., 242., 243., 244., 245., 246., 247., 248., 249./)
    litc=(/250.,251.12621,251.51005,252.,253.,254.,255.,256.,257.,258., &
     & 259., 260., 261., 262., 263., 264., 265., 266., 267., 268., &
     & 269., 270., 271., 272., 273., 274., 275., 276., 277., 278., &
     & 279., 280., 281., 282., 283., 284., 285., 286., 287., 288., &
     & 289., 290., 292., 294., 296., 298., 300., 302., 304., 306., &
     & 308., 310., 312., 314., 316., 318., 320., 322., 324., 326., &
     & 328., 330., 332., 334., 336., 338., 340., 342., 344., 346., &
     & 348., 350., 352., 354., 356., 358., 360., 362., 364., 364.70183, &
     & 366., 368., 370., 372., 374., 376., 378., 380., 382., 384., &
     & 386., 388., 390., 392., 394., 396., 398., 400., 402., 404., &
     & 406., 408., 410., 412., 414., 416., 418., 420., 422., 424., &
     & 426., 428., 430., 432., 434., 436., 438., 440., 442., 444., &
     & 446., 448., 450., 452., 454., 456., 458., 460., 462., 464., &
     & 466., 468., 470., 472., 474., 476., 478., 480., 482., 484., &
     & 486., 488., 490., 492., 494., 496., 498., 500., 502., 504., &
     & 506., 508., 510., 512., 514., 516., 518., 520., 522., 524., &
     & 526., 528., 530., 532., 534., 536., 538., 540., 542., 544., &
     & 546., 548., 550., 552., 554., 556., 558., 560., 562., 564., &
     & 566., 568., 570., 572., 574., 576., 578., 580., 582., 584./)
    litd=(/586., 588., 590., 592., 594., 596., 598., 600., 602., 604., &
     & 606., 608., 610., 612., 614., 616., 618., 620., 622., 624., &
     & 626., 628., 630., 632., 634., 636., 638., 640., 642., 644., &
     & 646., 648., 650., 652., 654., 656., 658., 660., 662., 664., &
     & 666., 668., 670., 672., 674., 676., 678., 680., 682., 684., &
     & 686., 688., 690., 692., 694., 696., 698., 700., 702., 704., &
     & 706., 708., 710., 712., 714., 716., 718., 720., 722., 724., &
     & 726., 728., 730., 732., 734., 736., 738., 740., 742., 744., &
     & 746., 748., 750., 752., 754., 756., 758., 760., 762., 764., &
     & 766., 768., 770., 772., 774., 776., 778., 780., 782., 784., &
     & 786., 788., 790., 792., 794., 796., 798., 800., 802., 804., &
     & 806., 808., 810., 812., 814., 816., 818., 820.58271, 822., 824., &
     & 826., 828., 830., 832., 834., 836., 838., 840., 842., 844., &
     & 846., 848., 850., 852., 854., 856., 858., 860., 862., 864., &
     & 866., 868., 870., 872., 874., 876., 878., 880., 882., 884., &
     & 886., 888., 890., 892., 894., 896., 898., 900., 902., 904., &
     & 906., 908., 910., 912., 914., 916., 918., 920., 922., 924., &
     & 926., 928., 930., 932., 934., 936., 938., 940., 942., 944., &
     & 946., 948., 950., 952., 954., 956., 958., 960., 962., 964./)
    lite=(/966., 968., 970., 972., 974., 976., 978., 980., 982., 984., &
     & 986., 988., 990., 992., 994., 996., 998.,1000.,1005.,1010., &
     & 1015.,1020.,1025.,1030.,1035.,1040.,1045.,1050.,1055.,1060., &
     & 1065.,1070.,1075.,1080.,1085.,1090.,1095.,1100.,1105.,1110., &
     & 1115.,1120.,1125.,1130.,1135.,1140.,1145.,1150.,1155.,1160., &
     & 1165.,1170.,1175.,1180.,1185.,1190.,1195.,1200.,1205.,1210., &
     & 1215.,1220.,1225.,1230.,1235.,1240.,1245.,1250.,1255.,1260., &
     & 1265.,1270.,1275.,1280.,1285.,1290.,1295.,1300.,1305.,1310., &
     & 1315.,1320.,1325.,1330.,1335.,1340.,1345.,1350.,1355.,1360., &
     & 1365.,1370.,1375.,1380.,1385.,1390.,1395.,1400.,1405.,1410., &
     & 1415.,1420.,1425.,1430.,1435.,1440.,1445.,1450.,1455.,1458.81670, &
     & 1465.,1470.,1475.,1480.,1485.,1490.,1495.,1500.,1505.,1510., &
     & 1515.,1520.,1525.,1530.,1535.,1540.,1545.,1550.,1555.,1560., &
     & 1565.,1570.,1575.,1580.,1585.,1590.,1595.,1600.,1610.,1620., &
     & 1630.,1640.,1650.,1660.,1670.,1680.,1690.,1700.,1710.,1720., &
     & 1730.,1740.,1750.,1760.,1770.,1780.,1790.,1800.,1810.,1820., &
     & 1830.,1840.,1850.,1860.,1870.,1880.,1890.,1900.,1910.,1920., &
     & 1930.,1940.,1950.,1960.,1970.,1980.,1990.,2000.,2010.,2020., &
     & 2030.,2040.,2050.,2060.,2070.,2080.,2090.,2100.,2110.,2120./)

    litf=(/ 2130.,2140.,2150.,2160.,2170.,2180.,2190.,2200.,2210.,2220., &
     & 2230.,2240.,2250.,2260.,2270.,2279.40330,2290.,2300.,2310.,2320., &
     & 2330.,2340.,2350.,2360.,2370.,2380.,2390.,2400.,2410.,2420., &
     & 2430.,2440.,2450.,2460.,2470.,2480.,2490.,2500.,2510.,2520., &
     & 2530.,2540.,2550.,2560.,2570.,2580.,2590.,2600.,2610.,2620., &
     & 2630.,2640.,2650.,2660.,2670.,2680.,2690.,2700.,2710.,2720., &
     & 2730.,2740.,2750.,2760.,2770.,2780.,2790.,2800.,2810.,2820., &
     & 2830.,2840.,2850.,2860.,2870.,2880.,2890.,2900.,2910.,2920., &
     & 2930.,2940.,2950.,2960.,2970.,2980.,2990.,3000.,3010.,3020., &
     & 3030.,3040.,3050.,3060.,3070.,3080.,3090.,3100.,3110.,3120., &
     & 3130.,3140.,3150.,3160.,3170.,3180.,3190.,3200.,3220.,3240., &
     & 3260.,3282.34820,3300.,3320.,3340.,3360.,3380.,3400.,3420.,3440., &
     & 3460.,3480.,3500.,3520.,3540.,3560.,3580.,3600.,3620.,3640., &
     & 3660.,3680.,3700.,3720.,3740.,3760.,3780.,3800.,3820.,3840., &
     & 3860.,3880.,3900.,3920.,3940.,3960.,3980.,4000.,4020.,4040., &
     & 4060.,4080.,4100.,4120.,4140.,4160.,4180.,4200.,4220.,4240., &
     & 4260.,4280.,4300.,4320.,4340.,4360.,4380.,4400.,4420.,4440., &
     & 4460.,4480.,4500.,4520.,4540.,4560.,4580.,4600.,4620.,4640., &
     & 4660.,4680.,4700.,4720.,4740.,4760.,4780.,4800.,4820.,4840./)
    litg=(/4860.,4880.,4900.,4920.,4940.,4960.,4980.,5000.,5020.,5040., &
     & 5060.,5080.,5100.,5120.,5140.,5160.,5180.,5200.,5220.,5240., &
     & 5260.,5280.,5300.,5320.,5340.,5360.,5380.,5400.,5420.,5440., &
     & 5460.,5480.,5500.,5520.,5540.,5560.,5580.,5600.,5620.,5640., &
     & 5660.,5680.,5700.,5720.,5740.,5760.,5780.,5800.,5820.,5840., &
     & 5860.,5880.,5900.,5920.,5940.,5960.,5980.,6000.,6020.,6040., &
     & 6060.,6080.,6100.,6120.,6140.,6160.,6180.,6200.,6220.,6240., &
     & 6260.,6280.,6300.,6320.,6340.,6360.,6380.,6400.,6440.,6480., &
     & 6520.,6560.,6600.,6640.,6680.,6720.,6760.,6800.,6840.,6880., &
     & 6920.,6960.,7000.,7040.,7080.,7120.,7160.,7200.,7240.,7280., &
     & 7320.,7360.,7400.,7440.,7480.,7520.,7560.,7600.,7640.,7680., &
     & 7720.,7760.,7800.,7840.,7880.,7920.,7960.,8000.,8040.,8080., &
     & 8120.,8160.,8200.,8240.,8280.,8320.,8360.,8400.,8440.,8480., &
     & 8520.,8560.,8600.,8640.,8680.,8720.,8760.,8800.,8840.,8880., &
     & 8920.,8960.,9000.,9040.,9080.,9120.,9160.,9200.,9240.,9280., &
     & 9320.,9360.,9400.,9440.,9480.,9520.,9560.,9600.,9640.,9680., &
     & 9720.,9760.,9800.,9840.,9880.,9920.,9960.,10000./)

      wavebig(1:95)=biga
      wavebig(96:220)=bigb
      wavebig(221:nsizebig+1) = bigc
      wavelit(1:95)=lita
      wavelit(96:285)=litb
      wavelit(286:475)=litc
      wavelit(476:665)=litd
      wavelit(666:855)=lite
      wavelit(856:1045)=litf
      wavelit(1046:nsizelit+1)=litg


   else

      open(unit = 95, file = './INPUT/bin_grid_sizes.dat', form= 'formatted', &
     &     status = 'old', access = 'sequential' )
!
      read(95,*) nsizebig
      do i = 1, nsizebig 
        read(95,*) wavebig(i) 
      end do 
      wavelit(1) = wavebig(1)  
      wavelit(2) = wavebig(2)
      nsizebig = nsizebig -1 
!
      close(unit = 95)
 
 
   end if



END SUBROUTINE def_binsize


END MODULE DFCOMM

