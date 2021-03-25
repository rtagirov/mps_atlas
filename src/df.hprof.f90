MODULE hprof
    use types
    use dfcomm
  
    implicit none

CONTAINS 

double precision function hprof4(n,m,j,wl0,delw,dopph)

!!     faster and stripped from synthe version
!!     from deane peterson
   use types
   use dfcomm
 
   implicit none
    include 'common.rhoxbl'
    include 'common.tempbl'
    include 'common.turbpr'
    include 'common.stateb'
    include 'common.depart'
    include 'common.opsblk'
    include 'common.ionsbl'
!-- careful k is defined as a constant in common.constb ' but here it is dummy integer


!--- Internal---------------!
    integer(kind=4) k, n,  m
    integer j, n1, m1, itemp1, mmn 
    real(kind=4)  xn, xn2, xm, xm2, xmn2, xm2mn2
    integer ifins, i, ipos, nwid, ifcore, icut 

    real(kind=dp) rydh, gnm, xne16, t4, t43, xknm, y1num, y1wht, freqnm, dbeta, wavenm
    real(kind=dp) c1con, c2con, radamp, resont, vdw, hwvdw, hwrad, stark
    real(kind=dp) hwstk, hwres, hwlor, hwdop, hfwid, dop, d, hhw, hprofres
    real(kind=dp) cutoff, spacing, freq22000,  cutfreq, hprofrad, hprofvdw
    real(kind=dp) top, wty1, y1scal, c1, c2, g1, gnot, beta, y1, y2, hproflor, gam
    real(kind=dp) prqs, freq15000, beta4000, prqsp4000, f, p1, fnx, cutoff4000, fns

   
    real(kind=dp)  delw,wl0, wl, freq, freq0, del
    real(kind=dp)  pp(maxd),fo(maxd),gcon1(maxd),gcon2(maxd),y1b(maxd),y1s(maxd)
    real(kind=dp)  c1d(maxd),c2d(maxd)
    real(kind=dp)  t3nhe(maxd),t3nh2(maxd)
    real(kind=dp), dimension(2,2) :: y1wtm
    real(kind=dp), dimension(4,3) :: xknmtb
    real(kind=dp), dimension(5,4) ::  stcomp, stcpwt

    real(kind=dp)  stalph(34),istal(4),lnghal(4),stwtal(34)
    real(kind=dp)  lncomp(4),finest(14),finswt(14)

    real(kind=dp) lorwing,asum(100),asumlyman(100)
    real(kind=dp)  cutoffh2plus(111),cutoffh2(91)
    real(kind=dp)  dopph(maxd)


!-----------------------------------------------------
!   external



    itemp1=0
    n1= 0 
    m1= 0
    rydh= 3.2880515e15

!     fine structure components for alpha lines in freq*10**-7
    stalph=(/-730.0,370.0,188.0,515.0,327.0,619.0,-772.0,-473.0,-369.0,120.0,  & 
     &   256.0,162.0,285.0,-161.0,-38.3,6.82,-174.0,-147.0,-101.0,-77.5,55.0,126.0, & 
     &   75.0,139.0,-60.0,3.7,27.0,-69.0,-42.0,-18.0,-5.5,-9.1,-33.0,-24./) 

!     alpha component weights
    stwtal=(/1.0,2.0,1.0,2.0,1.0,2.0,1.0,2.0,3.0,1.0,2.0,1.0,2.0,1.0,4.0,6.0,1.0,2.0, & 
     &      3.0,4.0,1.0,2.0,1.0,2.0,1.0,4.0,6.0,1.0,7.0,6.0,4.0,4.0,4.0,5./) 

    istal=(/1,3,10,21/)
    lnghal=(/2,7,11,14/) 
!
!     fine structure for m.eq.infinity in freq*10**-7

    stcomp(1:5,1) = (/0.0,0.0,0.0,0.0,0.0/)
    stcomp(1:5,2) = (/468.0,576.0,-522.0,0.0,0.0/)
    stcomp(1:5,3) = (/260.0,290.0,-33.0,-140.0,0.0/)
    stcomp(1:5,4) = (/140.0,150.0,18.0,-27.0,-51.0/) 

!     weights
    stcpwt(1:5,1) =(/1.0,0.0,0.0,0.0,0.0/)
    stcpwt(1:5,2) =(/1.0,1.0,2.0,0.0,0.0/)
    stcpwt(1:5,3) =(/1.0,1.0,4.0,3.0,0.0/) 
    stcpwt(1:5,4) =(/1.0,1.0,14.0,6.0,4.0/) 


    xknmtb(1:4,1) =(/.0001716,.009019,.1001,.5820/) 
    xknmtb(1:4,2) =(/.0005235,.01772,.171,.866/)  
    xknmtb(1:4,3) =(/.0008912,.02507,.223,1.02/)

    y1wtm(1:2,1)=(/1.e18,1.e17/)
    y1wtm(1:2,2)=(/1.e16,1.e14/)

!
    lncomp=(/1,3,4,5/) 
!     lyman alpha quasi h2+ cutoff
!     delta waveno =  -15000+100*(n-1) n=1,111   up to -4000
!
     cutoffh2plus=(/& 
     &-15.14,-15.06,-14.97,-14.88,-14.80,-14.71,-14.62,-14.53, & 
     &-14.44,-14.36,-14.27,-14.18,-14.09,-14.01,-13.92,-13.83, & 
     &-13.74,-13.65,-13.57,-13.48,-13.39,-13.30,-13.21,-13.13, & 
     &-13.04,-12.95,-12.86,-12.77,-12.69,-12.60,-12.51,-12.40, & 
     &-12.29,-12.15,-12.02,-11.90,-11.76,-11.63,-11.53,-11.41, & 
     &-11.30,-11.22,-11.15,-11.09,-11.07,-11.06,-11.07,-11.09, & 
     &-11.12,-11.16,-11.19,-11.21,-11.24,-11.27,-11.30,-11.33, & 
     &-11.36,-11.39,-11.42,-11.45,-11.48,-11.48,-11.48,-11.48, & 
     &-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48, & 
     &-11.48,-11.48,-11.48,-11.48,-11.41,-11.40,-11.39,-11.38, & 
     &-11.37,-11.36,-11.35,-11.34,-11.33,-11.32,-11.30,-11.29, & 
     &-11.28,-11.27,-11.27,-11.27,-11.26,-11.25,-11.24,-11.23, & 
     &-11.22,-11.21,-11.20,-11.19,-11.18,-11.17,-11.15,-11.14, & 
     &-11.13,-11.12,-11.11,-11.10,-11.09,-11.08,-11.07/) 
!     lyman alpha quasi h2 cutoff
!     delta waveno = -22000+200*(n-1)  n=1,91  -4000

    cutoffh2=(/& 
     &-13.64,-13.52,-13.39,-13.27,-13.14,-13.01,-12.87,-12.74, & 
     &-12.63,-12.56,-12.51,-12.48,-12.47,-12.49,-12.52,-12.55, & 
     &-12.57,-12.61,-12.65,-12.69,-12.72,-12.76,-12.79,-12.82, & 
     &-12.84,-12.85,-12.87,-12.90,-12.93,-12.94,-12.93,-12.95, & 
     &-12.95,-12.96,-12.97,-12.96,-12.96,-12.95,-12.95,-12.96, & 
     &-12.98,-12.99,-12.95,-12.96,-13.00,-13.00,-12.98,-12.97, & 
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, & 
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, & 
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-12.89,-12.88, & 
     &-12.87,-12.86,-12.85,-12.84,-12.83,-12.81,-12.80,-12.79, & 
     &-12.78,-12.76,-12.74,-12.72,-12.70,-12.68,-12.65,-12.62, & 
     &-12.59,-12.56,-12.53/) 
!
     asumlyman=(/& 
     & 0.000e+00, 6.265e+08, 1.897e+08, 8.126e+07, 4.203e+07, 2.450e+07, & 
     & 1.236e+07, 8.249e+06, 5.782e+06, 4.208e+06, 3.158e+06, 2.430e+06, & 
     & 1.910e+06, 1.567e+06, 1.274e+06, 1.050e+06, 8.752e+05, 7.373e+05, & 
     & 6.269e+05, 5.375e+05, 4.643e+05, 4.038e+05, 3.534e+05, 3.111e+05, & 
     & 2.752e+05, 2.447e+05, 2.185e+05, 1.959e+05, 1.763e+05, 1.593e+05, & 
     & 1.443e+05, 1.312e+05, 1.197e+05, 1.094e+05, 1.003e+05, 9.216e+04, & 
     & 8.489e+04, 7.836e+04, 7.249e+04, 6.719e+04, 6.239e+04, 5.804e+04, & 
     & 5.408e+04, 5.048e+04, 4.719e+04, 4.418e+04, 4.142e+04, 3.888e+04, &
     & 3.655e+04, 3.440e+04, 3.242e+04, 3.058e+04, 2.888e+04, 2.731e+04, & 
     & 2.585e+04, 2.449e+04, 2.322e+04, 2.204e+04, 2.094e+04, 1.991e+04, & 
     & 1.894e+04, 1.804e+04, 1.720e+04, 1.640e+04, 1.566e+04, 1.496e+04, & 
     & 1.430e+04, 1.368e+04, 1.309e+04, 1.254e+04, 1.201e+04, 1.152e+04, & 
     & 1.105e+04, 1.061e+04, 1.019e+04, 9.796e+03, 9.419e+03, 9.061e+03, & 
     & 8.721e+03, 8.398e+03, 8.091e+03, 7.799e+03, 7.520e+03, 7.255e+03, & 
     & 7.002e+03, 6.760e+03, 6.530e+03, 6.310e+03, 6.100e+03, 5.898e+03, & 
     & 5.706e+03, 5.522e+03, 5.346e+03, 5.177e+03, 5.015e+03, 4.860e+03, & 
     & 4.711e+03, 4.569e+03, 4.432e+03, 4.300e+03/) 

     asum=(/& 
     & 0.000e+00, 4.696e+08, 9.980e+07, 3.017e+07, 1.155e+07, 5.189e+06, & 
     & 2.616e+06, 1.437e+06, 8.444e+05, 5.234e+05, 3.389e+05, 2.275e+05, & 
     & 1.575e+05, 1.120e+05, 8.142e+04, 6.040e+04, 4.560e+04, 3.496e+04, & 
     & 2.719e+04, 2.141e+04, 1.711e+04, 1.377e+04, 1.119e+04, 9.166e+03, & 
     & 7.572e+03, 6.341e+03, 5.338e+03, 4.523e+03, 3.854e+03, 3.302e+03, & 
     & 2.844e+03, 2.460e+03, 2.138e+03, 1.866e+03, 1.635e+03, 1.438e+03, & 
     & 1.269e+03, 1.124e+03, 9.983e+02, 8.894e+02, 7.947e+02, 7.120e+02, & 
     & 6.396e+02, 5.759e+02, 5.198e+02, 4.703e+02, 4.263e+02, 3.873e+02, & 
     & 3.526e+02, 3.215e+02, 2.938e+02, 2.689e+02, 2.465e+02, 2.264e+02, & 
     & 2.082e+02, 1.918e+02, 1.769e+02, 1.634e+02, 1.512e+02, 1.400e+02, & 
     & 1.298e+02, 1.206e+02, 1.121e+02, 1.043e+02, 9.720e+01, 9.066e+01, & 
     & 8.465e+01, 7.912e+01, 7.403e+01, 6.933e+01, 6.498e+01, 6.097e+01, & 
     & 5.725e+01, 5.381e+01, 5.061e+01, 4.765e+01, 4.489e+01, 4.232e+01, & 
     & 3.994e+01, 3.771e+01, 3.563e+01, 3.369e+01, 3.188e+01, 3.019e+01, & 
     & 2.860e+01, 2.712e+01, 2.572e+01, 2.442e+01, 2.319e+01, 2.204e+01, & 
     & 2.096e+01, 1.994e+01, 1.898e+01, 1.808e+01, 1.722e+01, 1.642e+01, & 
     & 1.566e+01, 1.495e+01, 1.427e+01, 1.363e+01/) 

!------------------------------------------------------------------------!
!------------------------------------------------------------------------!
!------------------------------------------------------------------------!

  if(itemp .ne. itemp1) then 
!     set up depth vectors
    itemp1=itemp

    do k=1,nrhox

      xne16 = xne(k)**0.1666667
      pp(k) = xne16*0.08989/sqrt(t(k))
      fo(k) = xne16**4*1.25e-9
      y1b(k)= 2.0d0/(1.0d0+0.012/t(k)*sqrt(xne(k)/t(k)))
      t4 = t(k)/10000.0d0
      t43 =t4**0.3d0 
      y1s(k)= t43/xne16
      t3nhe(k)=t43*xnfhe(k,1)
      t3nh2(k)=t43*xnfh2(k)
      c1d(k)=fo(k)*78940.0d0/t(k)
      c2d(k)=fo(k)**2/5.96e-23/xne(k)
      gcon1(k)=0.2d0+0.09d0*sqrt(t4)/(1.0d0+xne(k)/1.0e13)
      gcon2(k)=0.2d0/(1.0d0+xne(k)/1.0e15)
    end do 

  end if 

!     set up for this line
  if((n.ne.n1) .or. (m.ne.m1)) then 

    n1=n
    m1=m
    mmn=m-n
    xn=dble(n) 
    xn2=xn*xn
    xm=dble(m) 
    xm2=xm*xm
    xmn2=xm2*xn2
    xm2mn2=xm2-xn2
    gnm= dble(xm2mn2)/dble(xmn2) 

    if(mmn.le.3.and.n.le.4) xknm=xknmtb(n,mmn)
    if(mmn.gt.3.or.n.gt.4) xknm=5.5e-5/gnm*xmn2/(1.+.13/float(mmn))

    y1num=320.

    if(m.eq.2) y1num=550.
    if(m.eq.3) y1num=380.

    y1wht=1.e13

    if(mmn.le.3)y1wht=1.e14

    if(mmn.le.2.and.n.le.2)y1wht=y1wtm(n,mmn)


    freqnm=rydh*gnm
    dbeta=2.99792458e18/freqnm**2/xknm
    wavenm=2.99792458e18/freqnm
    c1con=xknm/wavenm*gnm*dble(xm2mn2) 
    c2con=(xknm/wavenm)**2
    radamp=asum(n)+asum(m)
 
    if(n.eq.1)radamp=asumlyman(m)

    radamp=radamp/12.5664
    radamp=radamp/freqnm
    resont=hfnm(1,m)/(1.-1./dble(xm2))

    if(n.ne.1)resont=resont+hfnm(1,n)/(1.-1./xn2)
    resont=resont*3.579e-24*0.5773/gnm
    vdw=4.45e-26/gnm*(dble(xm2)*(7.*dble(xm2)+5.))**.4

!   guess that h2 is twice as strong as he in txnxn
    hwvdw=vdw*t3nhe(j)+2.*vdw*t3nh2(j)
    hwrad=radamp
    stark=1.6678e-18*freqnm*xknm

!!     fine structure components
!
!!      if(n.gt.4)then
    if(n.gt.4.or.m.gt.10)then

      ifins=1
      finest(1)=0.
      finswt(1)=1.

    else 
!
      if( mmn .ne. 1) then
!     use m.eq.inf structure
        ifins=lncomp(n)

        do  i=1,ifins
         finest(i)=stcomp(i,n)*1.e7
         finswt(i)=stcpwt(i,n)/xn2
        end do 
     
      else 

!     for alpha lines
        ifins=lnghal(n)
        ipos=istal(n)
     
         do i=1,ifins
           k=ipos-1+i
           finest(i)=stalph(k)*1.e7
           finswt(i)=stwtal(k)/xn2/3.
         end do 
     
      end if 
    end if
!
!
  end if


  wl=wl0*10.+delw*10.
  freq=2.99792458e18/wl
  freq0=2.99792458e17/wl0

!      del=abs(freq-freqnm)
  del=abs(freq-freq0)

!     wl in nm
  wl=wl/10.

!     these half-widths are really dnu/nu
  hwstk=stark*fo(j)
  hwvdw=vdw*t3nhe(j)+2.*vdw*t3nh2(j)
  hwrad=radamp

!     xnfph(j,1)*2 is the number in the ground state
  hwres=resont*xnfph(j,1)*2.
  hwlor=hwres+hwvdw+hwrad
  hwdop=dopph(j)


!     specify largest half width in case of core calc
!     nwid=1, doppler  =2, lorentz  =3, stark
   nwid=1



   if((hwdop.lt.hwstk) .or.  (hwdop.lt.hwlor)) then 
     nwid=2

     if(hwlor .lt. hwstk) then 
      nwid=3
     end if
   end if 

   hfwid=freqnm*max(hwdop,hwlor,hwstk)


!     sets flag if in a line core
!     hprofl=0.

   hprof4=0.
   ifcore=0


   if(abs(del).le.hfwid)ifcore=1

!      dop=freqnm*hwdop
   dop=freq0*hwdop

!
   if(ifcore.eq.1)go to (32,40,50),nwid
!-----------------------------------------------!
!     do doppler
!     put fine structure in doppler core

   32 do  i=1,ifins
        d=abs(freq-freq0-finest(i))/dop
        if(d .le. 7.0d0)hprof4=hprof4+h0tab(int(200.*d+1.5))*finswt(i)
      end do 
      
      if(ifcore.eq.1)return
 
   
!     do lorentz
   40 if(n.ne.1)go to 48
      if(m.ne.2)go to 48

!     lyman alpha
!     near center
!     modify old resonance broadening to match at 4000 cm-1

      hwres=hwres*4.
      hwlor=hwres+hwvdw+hwrad
      hhw=freqnm*hwlor

!     if (lambda<1277a)
      if(freq.gt.(82259.105-4000.)*2.99792458e10)then
        hprofres=hwres*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
        go to 44
      endif

!     only far red wing
!     data from n.f. allard, march 1997
!     insert lyman alpha cutoff a la n.f. allard and d. koester, a&a, 258,
!     464-468. 1992.

      cutoff=0.
      if(freq.lt.50000.*2.99792458e10)go to 43

!     tabulated at 200 cm-1 spacing
      spacing=200.*2.99792458e10
      freq22000=(82259.105-22000.)*2.99792458e10

      if(freq.lt.freq22000)then
        cutoff = (cutoffh2(2)-cutoffh2(1))/spacing*(freq-freq22000)+ & 
    &            cutoffh2(1)
      else

        icut=(freq-freq22000)/spacing
        icut=min(icut,89)
        cutfreq=icut*spacing+freq22000
        cutoff=(cutoffh2(icut+2)-cutoffh2(icut+1))/spacing*  & 
    &          (freq-cutfreq)+cutoffh2(icut+1)
      endif


      cutoff=(10.**(cutoff-14.))*xnfph(j,1)*2./2.99792458e10

   43 hprofres=cutoff*1.77245*dop
!
   44 hprofrad=0.


!     rayleighscattering except near doppler core 
!     crossover from absorption to rayleigh scattering in hrayop
!     lambda>/=1217a
!
      if(freq .gt. 2.4190611e15 .and. freq.lt..77*3.28805e15) then 
        hprofrad=hwrad*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
      end if


!     correction to lorentz profile   aller p.164   not used
!     hprofrad=hprofrad*4*freq**2/(freq**2+freqnm**2)

      hprofvdw=hwvdw*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop


!     van der waal cutoff for he and for h2
!     guess both 60000 cm-1=1.8e15 hz

      if(freq.lt.1.8e15)hprofvdw=0.

      hproflor=hprofres+hprofrad+hprofvdw
      hprof4=hprof4+hproflor

      if(ifcore.eq.1)return
      go to 50    
!------------------------------------------------
!     not lyman alpha

   48 hhw=freqnm*hwlor
      top=hhw

      if(n.ne.1.or.m.gt.5)go to 49
!     lyman beta
      if(m.eq.3.and.freq.gt..885*3.288051e15.and. & 
     &   freq.lt..890*3.288051e15)top=hhw-freqnm*hwrad

!     lyman gamma
      if(m.eq.4.and.freq.gt..936*3.288051e15.and. & 
     & freq.lt..938*3.288051e15)top=hhw-freqnm*hwrad

!     lyman delta
      if(m.eq.5.and.freq.gt..959*3.288051e15.and. &
     & freq.lt..961*3.288051e15)top=hhw-freqnm*hwrad

   49 hproflor=top/3.14159/(del**2+hhw**2)*1.77245*dop
      hprof4=hprof4+hproflor

      if(ifcore.eq.1)return


!  Do Stark
!---------------------------------------
   50 wty1=1./(1.+xne(j)/y1wht)
      y1scal=y1num*y1s(j)*wty1+y1b(j)*(1.-wty1)
      c1=c1d(j)*c1con*y1scal

      c2=c2d(j)*c2con
      g1=6.77*sqrt(c1)
      gnot=g1*max(0.0d0,.2114d0+log(sqrt(c2)/c1))*(1.-gcon1(j)-gcon2(j))
      beta=abs(del)/fo(j)*dbeta
      y1=c1*beta
      y2=c2*beta**2
      gam=gnot

      if(y2.le.1.e-4.and.y1.le.1.e-5)go to 51

!
!     change for dfsynthe for speed
!

! There was a bug -> y1 became in some cases smaller than 0.001!
!------------------------------ > make sure extab does not have
!-------------------------------> an index eq 0
!      gam = g1*(.5*extab(max(1,(nint(100*min(80.0d0,y1)))))+  & 
      gam = g1*(.5*extab(max(1,(nint(100*min(80.0d0,y1)))))+  &
     &    vcse1f(y1)-.5*vcse1f(y2))*  & 
     &    (1.-gcon1(j)/(1.+(90.*y1)**3)-gcon2(j)/(1.+2000.*y1))

      if(gam.le.1.e-20)gam=0.

   51 prqs=sofbet(beta,pp(j),n,m)

      if(m.gt.2)go to 53

!     assume quasistati! profile is half protons, half electrons
      prqs=prqs*.5
      cutoff=0.

!     lyman alpha quasi h2+ cutoff
!     data from n.f. allard, march 1997
      if(freq.lt.(82259.105-20000.)*2.99792458e10)go to 53
      if(freq.gt.(82259.105-4000.)*2.99792458e10)go to 52
!     tabulated at 100 cm-1 spacing
      freq15000=(82259.105-15000.)*2.99792458e10
      spacing=100.*2.99792458e10

      if(freq.lt.freq15000)then
        cutoff=(cutoffh2plus(2)-cutoffh2plus(1))/spacing*  & 
     &         (freq-freq15000)+cutoffh2plus(1)
      else
        icut=(freq-freq15000)/spacing
        icut=min(icut,109)
        cutfreq=icut*spacing+freq15000
        cutoff=(cutoffh2plus(icut+2)-cutoffh2plus(icut+1))/ & 
     & spacing*(freq-cutfreq)+cutoffh2plus(icut+1)
      endif

      cutoff=(10.**(cutoff-14.))/2.99792458e10*xnfph(j,2)
      hprof4=hprof4+cutoff*1.77245*dop
      go to 53

   52 beta4000=4000.*2.99792458e10/fo(j)*dbeta
      prqsp4000=sofbet(beta4000,pp(j),n,m)*.5/fo(j)*dbeta
      cutoff4000=(10.**(-11.07-14.))/2.99792458e10*xnfph(j,2)
      hprof4=hprof4+cutoff4000/prqsp4000*prqs/fo(j)*dbeta*1.77245*dop

   53 f=0.
      if(gam.gt.0.)f=gam/3.14159/(gam**2+beta**2)
      p1=(.9*y1)**2
      fns=(p1+.03*sqrt(y1))/(p1+1.)


!     same normalization as voigt function
      hprof4=hprof4+(prqs*(1.+fns)+f)/fo(j)*dbeta*1.77245*dop
      return

end function 



!---------------------------------------------------------------
!
double precision function hprof4m(n,m,j,wl0,delw,dopph)


!!     faster and stripped from synthe version
!!     from deane peterson
   use types
   use dfcomm

   implicit none
    include 'common.rhoxbl'
    include 'common.tempbl'
    include 'common.turbpr'
    include 'common.stateb'
    include 'common.depart'
    include 'common.opsblk'
    include 'common.ionsbl'
!-- careful k is defined as a constant in common.constb ' but here it is dummy integer


!--- Internal---------------!

!--- Internal---------------!
    integer k, n, j, m
    integer n1, m1, itemp1, mmn 
    real(kind=4)  xn, xn2, xm, xm2, xmn2, xm2mn2
    integer ifins, i, ipos, nwid, icut, ifcore

    real(kind=dp) rydh, gnm,  xne16, t4, t43, xknm, y1num, y1wht, freqnm, dbeta, wavenm
    real(kind=dp) c1con, c2con, radamp, resont, vdw, hwvdw, hwrad, stark
    real(kind=dp) hwstk, hwres, hwlor, hwdop, hfwid, dop, d, hhw, hprofres
    real(kind=dp) cutoff, spacing, freq22000,  cutfreq, hprofrad, hprofvdw
    real(kind=dp) top, wty1, y1scal, c1, c2, g1, gnot, beta, y1, y2, hproflor, gam
    real(kind=dp) prqs, freq15000, beta4000, prqsp4000, f, p1, fnx, cutoff4000, fns 


    real(kind=dp)  delw,wl0, wl, freq, freq0, del
    real(kind=dp)  pp(maxd),fo(maxd),gcon1(maxd),gcon2(maxd),y1b(maxd),y1s(maxd)
    real(kind=dp)  c1d(maxd),c2d(maxd)
    real(kind=dp)  t3nhe(maxd),t3nh2(maxd)
    real(kind=dp), dimension(2,2) :: y1wtm
    real(kind=dp), dimension(4,3) :: xknmtb
    real(kind=dp), dimension(5,4) ::  stcomp, stcpwt

    real(kind=dp)  stalph(34),istal(4),lnghal(4),stwtal(34)
    real(kind=dp)  lncomp(4),finest(14),finswt(14)

    real(kind=dp)  lorwing,asum(100),asumlyman(100)
    real(kind=dp)  cutoffh2plus(111),cutoffh2(91)
    real(kind=dp)  dopph(maxd)


!-----------------------------------------------------
!   external


!-----------------------------------------------------

    itemp1=0
    n1= 0
    m1= 0
    rydh= 3.2880515e15

!     fine structure components for alpha lines in freq*10**-7
    stalph=(/-730.0,370.0,188.0,515.0,327.0,619.0,-772.0,-473.0,-369.0,120.0,  &
     &   256.0,162.0,285.0,-161.0,-38.3,6.82,-174.0,-147.0,-101.0,-77.5,55.0,126.0, &
     &   75.0,139.0,-60.0,3.7,27.0,-69.0,-42.0,-18.0,-5.5,-9.1,-33.0,-24./)

!     alpha component weights
    stwtal=(/1.0,2.0,1.0,2.0,1.0,2.0,1.0,2.0,3.0,1.0,2.0,1.0,2.0,1.0,4.0,6.0,1.0,2.0, &
     &      3.0,4.0,1.0,2.0,1.0,2.0,1.0,4.0,6.0,1.0,7.0,6.0,4.0,4.0,4.0,5./)

    istal=(/1,3,10,21/)
    lnghal=(/2,7,11,14/)
!
!     fine structure for m.eq.infinity in freq*10**-7

    stcomp(1:5,1) = (/0.0,0.0,0.0,0.0,0.0/)
    stcomp(1:5,2) = (/468.0,576.0,-522.0,0.0,0.0/)
    stcomp(1:5,3) = (/260.0,290.0,-33.0,-140.0,0.0/)
    stcomp(1:5,4) = (/140.0,150.0,18.0,-27.0,-51.0/)

!     weights
    stcpwt(1:5,1) =(/1.0,0.0,0.0,0.0,0.0/)
    stcpwt(1:5,2) =(/1.0,1.0,2.0,0.0,0.0/)
    stcpwt(1:5,3) =(/1.0,1.0,4.0,3.0,0.0/) 
    stcpwt(1:5,4) =(/1.0,1.0,14.0,6.0,4.0/)   


    xknmtb(1:4,1) =(/.0001716,.009019,.1001,.5820/) 
    xknmtb(1:4,2) =(/.0005235,.01772,.171,.866/)  
    xknmtb(1:4,3) =(/.0008912,.02507,.223,1.02/)

    y1wtm(1:2,1)=(/1.e18,1.e17/)
    y1wtm(1:2,2)=(/1.e16,1.e14/)

!
    lncomp=(/1,3,4,5/)


!     lyman alpha quasi h2+ cutoff
!     delta waveno =  -15000+100*(n-1) n=1,111   up to -4000
!
     cutoffh2plus=(/&
     &-15.14,-15.06,-14.97,-14.88,-14.80,-14.71,-14.62,-14.53, &
     &-14.44,-14.36,-14.27,-14.18,-14.09,-14.01,-13.92,-13.83, &
     &-13.74,-13.65,-13.57,-13.48,-13.39,-13.30,-13.21,-13.13, &
     &-13.04,-12.95,-12.86,-12.77,-12.69,-12.60,-12.51,-12.40, &
     &-12.29,-12.15,-12.02,-11.90,-11.76,-11.63,-11.53,-11.41, &
     &-11.30,-11.22,-11.15,-11.09,-11.07,-11.06,-11.07,-11.09, &
     &-11.12,-11.16,-11.19,-11.21,-11.24,-11.27,-11.30,-11.33, &
     &-11.36,-11.39,-11.42,-11.45,-11.48,-11.48,-11.48,-11.48, &
     &-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48, &
     &-11.48,-11.48,-11.48,-11.48,-11.41,-11.40,-11.39,-11.38, &
     &-11.37,-11.36,-11.35,-11.34,-11.33,-11.32,-11.30,-11.29, &
     &-11.28,-11.27,-11.27,-11.27,-11.26,-11.25,-11.24,-11.23, &
     &-11.22,-11.21,-11.20,-11.19,-11.18,-11.17,-11.15,-11.14, &
     &-11.13,-11.12,-11.11,-11.10,-11.09,-11.08,-11.07/)
!     lyman alpha quasi h2 cutoff
!     delta waveno = -22000+200*(n-1)  n=1,91  -4000

    cutoffh2=(/&
     &-13.64,-13.52,-13.39,-13.27,-13.14,-13.01,-12.87,-12.74, &
     &-12.63,-12.56,-12.51,-12.48,-12.47,-12.49,-12.52,-12.55, &
     &-12.57,-12.61,-12.65,-12.69,-12.72,-12.76,-12.79,-12.82, &
     &-12.84,-12.85,-12.87,-12.90,-12.93,-12.94,-12.93,-12.95, &
     &-12.95,-12.96,-12.97,-12.96,-12.96,-12.95,-12.95,-12.96, &
     &-12.98,-12.99,-12.95,-12.96,-13.00,-13.00,-12.98,-12.97, &
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, &
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, &
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-12.89,-12.88, &
     &-12.87,-12.86,-12.85,-12.84,-12.83,-12.81,-12.80,-12.79, &
     &-12.78,-12.76,-12.74,-12.72,-12.70,-12.68,-12.65,-12.62, &
     &-12.59,-12.56,-12.53/)
!
     asumlyman=(/&
     & 0.000e+00, 6.265e+08, 1.897e+08, 8.126e+07, 4.203e+07, 2.450e+07, &
     & 1.236e+07, 8.249e+06, 5.782e+06, 4.208e+06, 3.158e+06, 2.430e+06, &
     & 1.910e+06, 1.567e+06, 1.274e+06, 1.050e+06, 8.752e+05, 7.373e+05, &
     & 6.269e+05, 5.375e+05, 4.643e+05, 4.038e+05, 3.534e+05, 3.111e+05, &
     & 2.752e+05, 2.447e+05, 2.185e+05, 1.959e+05, 1.763e+05, 1.593e+05, &
     & 1.443e+05, 1.312e+05, 1.197e+05, 1.094e+05, 1.003e+05, 9.216e+04, &
     & 8.489e+04, 7.836e+04, 7.249e+04, 6.719e+04, 6.239e+04, 5.804e+04, &
     & 5.408e+04, 5.048e+04, 4.719e+04, 4.418e+04, 4.142e+04, 3.888e+04, &
     & 3.655e+04, 3.440e+04, 3.242e+04, 3.058e+04, 2.888e+04, 2.731e+04, &
     & 2.585e+04, 2.449e+04, 2.322e+04, 2.204e+04, 2.094e+04, 1.991e+04, &
     & 1.894e+04, 1.804e+04, 1.720e+04, 1.640e+04, 1.566e+04, 1.496e+04, &
     & 1.430e+04, 1.368e+04, 1.309e+04, 1.254e+04, 1.201e+04, 1.152e+04, &
     & 1.105e+04, 1.061e+04, 1.019e+04, 9.796e+03, 9.419e+03, 9.061e+03, &
     & 8.721e+03, 8.398e+03, 8.091e+03, 7.799e+03, 7.520e+03, 7.255e+03, &
     & 7.002e+03, 6.760e+03, 6.530e+03, 6.310e+03, 6.100e+03, 5.898e+03, &
     & 5.706e+03, 5.522e+03, 5.346e+03, 5.177e+03, 5.015e+03, 4.860e+03, &
     & 4.711e+03, 4.569e+03, 4.432e+03, 4.300e+03/)

     asum=(/&
     & 0.000e+00, 4.696e+08, 9.980e+07, 3.017e+07, 1.155e+07, 5.189e+06, & 
     & 2.616e+06, 1.437e+06, 8.444e+05, 5.234e+05, 3.389e+05, 2.275e+05, & 
     & 1.575e+05, 1.120e+05, 8.142e+04, 6.040e+04, 4.560e+04, 3.496e+04, & 
     & 2.719e+04, 2.141e+04, 1.711e+04, 1.377e+04, 1.119e+04, 9.166e+03, & 
     & 7.572e+03, 6.341e+03, 5.338e+03, 4.523e+03, 3.854e+03, 3.302e+03, & 
     & 2.844e+03, 2.460e+03, 2.138e+03, 1.866e+03, 1.635e+03, 1.438e+03, & 
     & 1.269e+03, 1.124e+03, 9.983e+02, 8.894e+02, 7.947e+02, 7.120e+02, & 
     & 6.396e+02, 5.759e+02, 5.198e+02, 4.703e+02, 4.263e+02, 3.873e+02, & 
     & 3.526e+02, 3.215e+02, 2.938e+02, 2.689e+02, 2.465e+02, 2.264e+02, & 
     & 2.082e+02, 1.918e+02, 1.769e+02, 1.634e+02, 1.512e+02, 1.400e+02, & 
     & 1.298e+02, 1.206e+02, 1.121e+02, 1.043e+02, 9.720e+01, 9.066e+01, & 
     & 8.465e+01, 7.912e+01, 7.403e+01, 6.933e+01, 6.498e+01, 6.097e+01, & 
     & 5.725e+01, 5.381e+01, 5.061e+01, 4.765e+01, 4.489e+01, 4.232e+01, & 
     & 3.994e+01, 3.771e+01, 3.563e+01, 3.369e+01, 3.188e+01, 3.019e+01, & 
     & 2.860e+01, 2.712e+01, 2.572e+01, 2.442e+01, 2.319e+01, 2.204e+01, & 
     & 2.096e+01, 1.994e+01, 1.898e+01, 1.808e+01, 1.722e+01, 1.642e+01, & 
     & 1.566e+01, 1.495e+01, 1.427e+01, 1.363e+01/)



!------------------------------------------------------------------------!
!------------------------------------------------------------------------!
!------------------------------------------------------------------------!



    if(itemp.eq.itemp1)go to 20
!     set up depth vectors
      itemp1=itemp

    do  k=1,nrhox
      xne16=xne(k)**.1666667
      pp(k)=xne16*.08989/sqrt(t(k))
      fo(k)=xne16**4*1.25e-9
      y1b(k)=2./(1.+.012/t(k)*sqrt(xne(k)/t(k)))
      t4=t(k)/10000.
      t43=t4**.3
      y1s(k)=t43/xne16
!      t3nhe(k)=t43*xnfphe(k,1)
      t3nhe(k)=t43*xnfhe(k,1)
      t3nh2(k)=t43*xnfh2(k)
      c1d(k)=fo(k)*78940./t(k)
      c2d(k)=fo(k)**2/5.96e-23/xne(k)
      gcon1(k)=.2+.09*sqrt(t4)/(1.+xne(k)/1.e13)
      gcon2(k)=.2/(1.+xne(k)/1.e15)
    end do 
!     set up for this line
   20 if(n.eq.n1.and.m.eq.m1)go to 30
      n1=n
      m1=m
      mmn=m-n
      xn=dble(n) 
      xn2=xn*xn
      xm=dble(m) 
      xm2=xm*xm
      xmn2=xm2*xn2
      xm2mn2=xm2-xn2
      gnm= dble(xm2mn2)/dble(xmn2) 
      if(mmn.le.3.and.n.le.4)xknm=xknmtb(n,mmn)
      if(mmn.gt.3.or.n.gt.4)xknm=5.5e-5/gnm*xmn2/(1.+.13/float(mmn))
      y1num=320.
      if(m.eq.2)y1num=550.
      if(m.eq.3)y1num=380.
      y1wht=1.e13
      if(mmn.le.3)y1wht=1.e14
      if(mmn.le.2.and.n.le.2)y1wht=y1wtm(n,mmn)
      freqnm=rydh*gnm
      dbeta=2.99792458e18/freqnm**2/xknm
      wavenm=2.99792458e18/freqnm


      c1con=xknm/wavenm*gnm*dble(xm2mn2) 
      c2con=(xknm/wavenm)**2
      radamp=asum(n)+asum(m)
      if(n.eq.1)radamp=asumlyman(m)
      radamp=radamp/12.5664
      radamp=radamp/freqnm
      resont=hfnm(1,m)/(1.-1./dble(xm2))
      if(n.ne.1)resont=resont+hfnm(1,n)/(1.-1./xn2)
!      fudge to baschek*2
      resont=resont*3.579e-24*0.5773/gnm
      vdw=4.45e-26/gnm*(dble(xm2)*(7.*dble(xm2)+5.))**.4
!     guess that h2 is twice as strong as he in txnxn
      hwvdw=vdw*t3nhe(j)+2.*vdw*t3nh2(j)
      hwrad=radamp
      stark=1.6678e-18*freqnm*xknm
!     fine structure components
      if(n.gt.4.or.m.gt.10)then
      ifins=1
      finest(1)=0.
      finswt(1)=1.
      go to 30
      endif
      
      if(mmn .eq. 1) go to 22
!     use m.eq.inf structure
      ifins=lncomp(n)
      do 21 i=1,ifins
      finest(i)=stcomp(i,n)*1.e7
      finswt(i)=stcpwt(i,n)/xn2
   21 continue
      go to 30
!     for alpha lines
   22 ifins=lnghal(n)
      ipos=istal(n)
      do 23 i=1,ifins
      k=ipos-1+i
      finest(i)=stalph(k)*1.e7
      finswt(i)=stwtal(k)/xn2/3.
   23 continue
!     now do this depth
   30 wl=wl0*10.+delw*10.
      freq=2.99792458e18/wl
      freq0=2.99792458e17/wl0
!      del=abs(freq-freqnm)
      del=abs(freq-freq0)
!     wl in nm
      wl=wl/10.
!     these half-widths are really dnu/nu
      hwstk=stark*fo(j)
      hwvdw=vdw*t3nhe(j)+2.*vdw*t3nh2(j)
      hwrad=radamp
!     xnfph(j,1)*2 is the number in the ground state
      hwres=resont*xnfph(j,1)*2.
      hwlor=hwres+hwvdw+hwrad
      hwdop=dopph(j)
!     specify largest half width in case of core calc
!     nwid=1, doppler  =2, lorentz  =3, stark
      nwid=1
      if(hwdop.ge.hwstk.and.hwdop.ge.hwlor)go to 31
      nwid=2
      if(hwlor.ge.hwstk)go to 31
      nwid=3
   31 hfwid=freqnm*max(hwdop,hwlor,hwstk)
!     sets flag if in a line core
      hprof4m=0.
      ifcore=0
      if(abs(del).le.hfwid)ifcore=1
!      dop=freqnm*hwdop
      dop=freq0*hwdop
      if(ifcore.eq.1)go to (32,40,50),nwid
!     do doppler
!     put fine structure in doppler core
   32 do 33 i=1,ifins
      d=abs(freq-freq0-finest(i))/dop
      if(d.le.7.)hprof4m=hprof4m+h0tab(int(200.*d+1.5))*finswt(i)
   33 continue
      if(ifcore.eq.1)return
!     do lorentz
   40 if(n.ne.1)go to 48
      if(m.ne.2)go to 48
!     lyman alpha
!     near center
!     modify old resonance broadening to match at 4000 cm-1
      hwres=hwres*4.
      hwlor=hwres+hwvdw+hwrad
      hhw=freqnm*hwlor
!     if (lambda<1277a)
      if(freq.gt.(82259.105-4000.)*2.99792458e10)then
      hprofres=hwres*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
      go to 44
      endif
!     only far red wing
!     data from n.f. allard, march 1997
!     insert lyman alpha cutoff a la n.f. allard and d. koester, a&a, 258,
!     464-468. 1992.
      cutoff=0.
      if(freq.lt.50000.*2.99792458e10)go to 43
!     tabulated at 200 cm-1 spacing
      spacing=200.*2.99792458e10
      freq22000=(82259.105-22000.)*2.99792458e10
      if(freq.lt.freq22000)then
      cutoff=(cutoffh2(2)-cutoffh2(1))/spacing*(freq-freq22000)+ & 
     &       cutoffh2(1)
      else
      icut=(freq-freq22000)/spacing
      icut=min(icut,89)
      cutfreq=icut*spacing+freq22000
      cutoff=(cutoffh2(icut+2)-cutoffh2(icut+1))/spacing* & 
     &   (freq-cutfreq)+cutoffh2(icut+1)
      endif
      cutoff=(10.**(cutoff-14.))*xnfph(j,1)*2./2.99792458e10
   43 hprofres=cutoff*1.77245*dop
!
   44 hprofrad=0.
!     rayleighscattering except near doppler core 
!     crossover from absorption to rayleigh scattering in hrayop
!     lambda>/=1217a
      if(freq.gt.2.4190611e15.and.freq.lt..77*3.28805e15) & 
     &  hprofrad=hwrad*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
!     correction to lorentz profile   aller p.164   not used
!
      hprofvdw=hwvdw*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
!     van der waal cutoff for he and for h2
!     guess both 60000 cm-1=1.8e15 hz
      if(freq.lt.1.8e15)hprofvdw=0.
      hproflor=hprofres+hprofrad+hprofvdw
      hprof4m=hprof4m+hproflor
      if(ifcore.eq.1)return
      go to 50
!
!     not lyman alpha
   48 hhw=freqnm*hwlor
      top=hhw
      if(n.ne.1.or.m.gt.5)go to 49
!     lyman beta
      if(m.eq.3.and.freq.gt..885*3.288051e15.and. & 
     &  freq.lt..890*3.288051e15)top=hhw-freqnm*hwrad

!     lyman gamma
      if(m.eq.4.and.freq.gt..936*3.288051e15.and. & 
     & freq.lt..938*3.288051e15)top=hhw-freqnm*hwrad

!     lyman delta
      if(m.eq.5.and.freq.gt..959*3.288051e15.and. &
     & freq.lt..961*3.288051e15)top=hhw-freqnm*hwrad

   49 hproflor=top/3.14159/(del**2+hhw**2)*1.77245*dop
      hprof4m=hprof4m+hproflor
      if(ifcore.eq.1)return
!     do stark
   50 wty1=1./(1.+xne(j)/y1wht)
      y1scal=y1num*y1s(j)*wty1+y1b(j)*(1.-wty1)

      c1=c1d(j)*c1con*y1scal
      c2=c2d(j)*c2con
      g1=6.77*sqrt(c1)

      gnot=g1*max(0.0d0,.2114d0+log(sqrt(c2)/c1))*(1.-gcon1(j)-gcon2(j))
      beta=abs(del)/fo(j)*dbeta
      y1=c1*beta
      y2=c2*beta**2
      gam=gnot

      if(y2.le.1.e-4.and.y1.le.1.e-5)go to 51

!     change for dfsynthe for speed

      gam=g1*(.5*extab(nint(100.*min(80.0d0,y1)))+  & 
     &    vcse1f(y1)-.5*vcse1f(y2))* & 
     &   (1.-gcon1(j)/(1.+(90.*y1)**3)-gcon2(j)/(1.+2000.*y1))

      if(gam.le.1.e-20)gam=0.
   51 prqs=sofbet(beta,pp(j),n,m)
      if(m.gt.2)go to 53
!     assume quasistatic profile is half protons, half electrons
      prqs=prqs*.5
      cutoff=0.
!     lyman alpha quasi h2+ cutoff
!     data from n.f. allard, march 1997
      if(freq.lt.(82259.105-20000.)*2.99792458e10)go to 53
      if(freq.gt.(82259.105-4000.)*2.99792458e10)go to 52
!     tabulated at 100 cm-1 spacing
      freq15000=(82259.105-15000.)*2.99792458e10
      spacing=100.*2.99792458e10
      if(freq.lt.freq15000)then
      cutoff=(cutoffh2plus(2)-cutoffh2plus(1))/spacing*  & 
     & (freq-freq15000)+cutoffh2plus(1)
      else
      icut=(freq-freq15000)/spacing
      icut=min(icut,109)
      cutfreq=icut*spacing+freq15000
      cutoff=(cutoffh2plus(icut+2)-cutoffh2plus(icut+1))/ & 
     & spacing*(freq-cutfreq)+cutoffh2plus(icut+1)
      endif
!     xnfph(j,2)=xnfh(j,2)
      cutoff=(10.**(cutoff-14.))/2.99792458e10*xnfph(j,2)
      hprof4m=hprof4m+cutoff*1.77245*dop
      go to 53
   52 beta4000=4000.*2.99792458e10/fo(j)*dbeta
      prqsp4000=sofbet(beta4000,pp(j),n,m)*.5/fo(j)*dbeta
      cutoff4000=(10.**(-11.07-14.))/2.99792458e10*xnfph(j,2)
      hprof4m=hprof4m+cutoff4000/prqsp4000*prqs/fo(j)*dbeta*1.77245*dop
   53 f=0.
      if(gam.gt.0.)f=gam/3.14159/(gam**2+beta**2)
      p1=(.9*y1)**2
      fns=(p1+.03*sqrt(y1))/(p1+1.)
!     same normalization as voigt function
      hprof4m=hprof4m+(prqs*(1.+fns)+f)/fo(j)*dbeta*1.77245*dop
      return
  end function  

double precision  function hprof4p(n,m,j,wl0,delw,dopph)

   use types
   use dfcomm

   implicit none
    include 'common.rhoxbl'
    include 'common.tempbl'
    include 'common.turbpr'
    include 'common.stateb'
    include 'common.depart'
    include 'common.opsblk'
    include 'common.ionsbl'
!-- careful k is defined as a constant in common.constb ' but here it is dummy integer


!--- Internal---------------!
!--- Internal---------------!
    integer k, n, j, m
    integer n1, m1, itemp1, mmn
    real(kind=4)  xn, xn2, xm, xm2, xmn2, xm2mn2
    integer ifins, i, ipos, nwid, icut,  ifcore

    real(kind=dp) rydh, gnm, xne16, t4, t43, xknm, y1num, y1wht, freqnm, dbeta, wavenm
    real(kind=dp) c1con, c2con, radamp, resont, vdw, hwvdw, hwrad, stark
    real(kind=dp) hwstk, hwres, hwlor, hwdop, hfwid, dop, d, hhw, hprofres
    real(kind=dp) cutoff, spacing, freq22000,  cutfreq, hprofrad, hprofvdw
    real(kind=dp) top, wty1, y1scal, c1, c2, g1, gnot, beta, y1, y2, hproflor, gam
    real(kind=dp) prqs, freq15000, beta4000, prqsp4000, f, p1, fnx, cutoff4000, fns 


    real(kind=dp)  delw,wl0, wl, freq, freq0, del
    real(kind=dp)  pp(maxd),fo(maxd),gcon1(maxd),gcon2(maxd),y1b(maxd),y1s(maxd)
    real(kind=dp)  c1d(maxd),c2d(maxd)
    real(kind=dp)  t3nhe(maxd),t3nh2(maxd)
    real(kind=dp), dimension(2,2) :: y1wtm
    real(kind=dp), dimension(4,3) :: xknmtb
    real(kind=dp), dimension(5,4) ::  stcomp, stcpwt

    real(kind=dp)  stalph(34),istal(4),lnghal(4),stwtal(34)
    real(kind=dp)  lncomp(4),finest(14),finswt(14)

    real(kind=dp) lorwing,asum(100),asumlyman(100)
    real(kind=dp) cutoffh2plus(111),cutoffh2(91)
    real(kind=dp) dopph(maxd)


!-----------------------------------------------------
!   external



!-----------------------------------------------------

    itemp1=0
    n1= 0
    m1= 0
    rydh= 3.2880515e15

!     fine structure components for alpha lines in freq*10**-7
    stalph=(/-730.0,370.0,188.0,515.0,327.0,619.0,-772.0,-473.0,-369.0,120.0,  &
     &   256.0,162.0,285.0,-161.0,-38.3,6.82,-174.0,-147.0,-101.0,-77.5,55.0,126.0, &
     &   75.0,139.0,-60.0,3.7,27.0,-69.0,-42.0,-18.0,-5.5,-9.1,-33.0,-24./)

!     alpha component weights
    stwtal=(/1.0,2.0,1.0,2.0,1.0,2.0,1.0,2.0,3.0,1.0,2.0,1.0,2.0,1.0,4.0,6.0,1.0,2.0, &
     &      3.0,4.0,1.0,2.0,1.0,2.0,1.0,4.0,6.0,1.0,7.0,6.0,4.0,4.0,4.0,5./)

    istal=(/1,3,10,21/)
    lnghal=(/2,7,11,14/)
!
!     fine structure for m.eq.infinity in freq*10**-7

    stcomp(1:5,1) = (/0.0,0.0,0.0,0.0,0.0/)
    stcomp(1:5,2) = (/468.0,576.0,-522.0,0.0,0.0/)
    stcomp(1:5,3) = (/260.0,290.0,-33.0,-140.0,0.0/)
    stcomp(1:5,4) = (/140.0,150.0,18.0,-27.0,-51.0/)

!     weights
    stcpwt(1:5,1) =(/1.0,0.0,0.0,0.0,0.0/)
    stcpwt(1:5,2) =(/1.0,1.0,2.0,0.0,0.0/)
    stcpwt(1:5,3) =(/1.0,1.0,4.0,3.0,0.0/) 
    stcpwt(1:5,4) =(/1.0,1.0,14.0,6.0,4.0/)   


    xknmtb(1:4,1) =(/.0001716,.009019,.1001,.5820/) 
    xknmtb(1:4,2) =(/.0005235,.01772,.171,.866/)  
    xknmtb(1:4,3) =(/.0008912,.02507,.223,1.02/)

    y1wtm(1:2,1)=(/1.e18,1.e17/)
    y1wtm(1:2,2)=(/1.e16,1.e14/)

!
    lncomp=(/1,3,4,5/)


!     lyman alpha quasi h2+ cutoff
!     delta waveno =  -15000+100*(n-1) n=1,111   up to -4000
!
     cutoffh2plus=(/&
     &-15.14,-15.06,-14.97,-14.88,-14.80,-14.71,-14.62,-14.53, &
     &-14.44,-14.36,-14.27,-14.18,-14.09,-14.01,-13.92,-13.83, &
     &-13.74,-13.65,-13.57,-13.48,-13.39,-13.30,-13.21,-13.13, &
     &-13.04,-12.95,-12.86,-12.77,-12.69,-12.60,-12.51,-12.40, &
     &-12.29,-12.15,-12.02,-11.90,-11.76,-11.63,-11.53,-11.41, &
     &-11.30,-11.22,-11.15,-11.09,-11.07,-11.06,-11.07,-11.09, &
     &-11.12,-11.16,-11.19,-11.21,-11.24,-11.27,-11.30,-11.33, &
     &-11.36,-11.39,-11.42,-11.45,-11.48,-11.48,-11.48,-11.48, &
     &-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48,-11.48, &
     &-11.48,-11.48,-11.48,-11.48,-11.41,-11.40,-11.39,-11.38, &
     &-11.37,-11.36,-11.35,-11.34,-11.33,-11.32,-11.30,-11.29, &
     &-11.28,-11.27,-11.27,-11.27,-11.26,-11.25,-11.24,-11.23, &
     &-11.22,-11.21,-11.20,-11.19,-11.18,-11.17,-11.15,-11.14, &
     &-11.13,-11.12,-11.11,-11.10,-11.09,-11.08,-11.07/)
!     lyman alpha quasi h2 cutoff
!     delta waveno = -22000+200*(n-1)  n=1,91  -4000

    cutoffh2=(/&
     &-13.64,-13.52,-13.39,-13.27,-13.14,-13.01,-12.87,-12.74, &
     &-12.63,-12.56,-12.51,-12.48,-12.47,-12.49,-12.52,-12.55, &
     &-12.57,-12.61,-12.65,-12.69,-12.72,-12.76,-12.79,-12.82, &
     &-12.84,-12.85,-12.87,-12.90,-12.93,-12.94,-12.93,-12.95, &
     &-12.95,-12.96,-12.97,-12.96,-12.96,-12.95,-12.95,-12.96, &
     &-12.98,-12.99,-12.95,-12.96,-13.00,-13.00,-12.98,-12.97, &
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, &
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-13.00, &
     &-13.00,-13.00,-13.00,-13.00,-13.00,-13.00,-12.89,-12.88, &
     &-12.87,-12.86,-12.85,-12.84,-12.83,-12.81,-12.80,-12.79, &
     &-12.78,-12.76,-12.74,-12.72,-12.70,-12.68,-12.65,-12.62, &
     &-12.59,-12.56,-12.53/)
!
     asumlyman=(/&
     & 0.000e+00, 6.265e+08, 1.897e+08, 8.126e+07, 4.203e+07, 2.450e+07, &
     & 1.236e+07, 8.249e+06, 5.782e+06, 4.208e+06, 3.158e+06, 2.430e+06, &
     & 1.910e+06, 1.567e+06, 1.274e+06, 1.050e+06, 8.752e+05, 7.373e+05, &
     & 6.269e+05, 5.375e+05, 4.643e+05, 4.038e+05, 3.534e+05, 3.111e+05, &
     & 2.752e+05, 2.447e+05, 2.185e+05, 1.959e+05, 1.763e+05, 1.593e+05, &
     & 1.443e+05, 1.312e+05, 1.197e+05, 1.094e+05, 1.003e+05, 9.216e+04, &
     & 8.489e+04, 7.836e+04, 7.249e+04, 6.719e+04, 6.239e+04, 5.804e+04, &
     & 5.408e+04, 5.048e+04, 4.719e+04, 4.418e+04, 4.142e+04, 3.888e+04, &
     & 3.655e+04, 3.440e+04, 3.242e+04, 3.058e+04, 2.888e+04, 2.731e+04, &
     & 2.585e+04, 2.449e+04, 2.322e+04, 2.204e+04, 2.094e+04, 1.991e+04, &
     & 1.894e+04, 1.804e+04, 1.720e+04, 1.640e+04, 1.566e+04, 1.496e+04, &
     & 1.430e+04, 1.368e+04, 1.309e+04, 1.254e+04, 1.201e+04, 1.152e+04, &
     & 1.105e+04, 1.061e+04, 1.019e+04, 9.796e+03, 9.419e+03, 9.061e+03, &
     & 8.721e+03, 8.398e+03, 8.091e+03, 7.799e+03, 7.520e+03, 7.255e+03, &
     & 7.002e+03, 6.760e+03, 6.530e+03, 6.310e+03, 6.100e+03, 5.898e+03, &
     & 5.706e+03, 5.522e+03, 5.346e+03, 5.177e+03, 5.015e+03, 4.860e+03, &
     & 4.711e+03, 4.569e+03, 4.432e+03, 4.300e+03/)


     asum=(/&
     & 0.000e+00, 4.696e+08, 9.980e+07, 3.017e+07, 1.155e+07, 5.189e+06, & 
     & 2.616e+06, 1.437e+06, 8.444e+05, 5.234e+05, 3.389e+05, 2.275e+05, & 
     & 1.575e+05, 1.120e+05, 8.142e+04, 6.040e+04, 4.560e+04, 3.496e+04, & 
     & 2.719e+04, 2.141e+04, 1.711e+04, 1.377e+04, 1.119e+04, 9.166e+03, & 
     & 7.572e+03, 6.341e+03, 5.338e+03, 4.523e+03, 3.854e+03, 3.302e+03, & 
     & 2.844e+03, 2.460e+03, 2.138e+03, 1.866e+03, 1.635e+03, 1.438e+03, & 
     & 1.269e+03, 1.124e+03, 9.983e+02, 8.894e+02, 7.947e+02, 7.120e+02, & 
     & 6.396e+02, 5.759e+02, 5.198e+02, 4.703e+02, 4.263e+02, 3.873e+02, & 
     & 3.526e+02, 3.215e+02, 2.938e+02, 2.689e+02, 2.465e+02, 2.264e+02, & 
     & 2.082e+02, 1.918e+02, 1.769e+02, 1.634e+02, 1.512e+02, 1.400e+02, & 
     & 1.298e+02, 1.206e+02, 1.121e+02, 1.043e+02, 9.720e+01, 9.066e+01, & 
     & 8.465e+01, 7.912e+01, 7.403e+01, 6.933e+01, 6.498e+01, 6.097e+01, & 
     & 5.725e+01, 5.381e+01, 5.061e+01, 4.765e+01, 4.489e+01, 4.232e+01, & 
     & 3.994e+01, 3.771e+01, 3.563e+01, 3.369e+01, 3.188e+01, 3.019e+01, & 
     & 2.860e+01, 2.712e+01, 2.572e+01, 2.442e+01, 2.319e+01, 2.204e+01, & 
     & 2.096e+01, 1.994e+01, 1.898e+01, 1.808e+01, 1.722e+01, 1.642e+01, & 
     & 1.566e+01, 1.495e+01, 1.427e+01, 1.363e+01/)



!------------------------------------------------------------------------!
!------------------------------------------------------------------------!
!------------------------------------------------------------------------!



    if(itemp.eq.itemp1)go to 20
!     set up depth vectors
      itemp1=itemp
    do  k=1,nrhox
      xne16=xne(k)**.1666667
      pp(k)=xne16*.08989/sqrt(t(k))
      fo(k)=xne16**4*1.25e-9
      y1b(k)=2./(1.+.012/t(k)*sqrt(xne(k)/t(k)))
      t4=t(k)/10000.
      t43=t4**.3
      y1s(k)=t43/xne16
!      t3nhe(k)=t43*xnfphe(k,1)
      t3nhe(k)=t43*xnfhe(k,1)
      t3nh2(k)=t43*xnfh2(k)
      c1d(k)=fo(k)*78940./t(k)
      c2d(k)=fo(k)**2/5.96e-23/xne(k)
      gcon1(k)=.2+.09*sqrt(t4)/(1.+xne(k)/1.e13)
      gcon2(k)=.2/(1.+xne(k)/1.e15)
    end do    
!     set up for this line
   20 if(n.eq.n1.and.m.eq.m1)go to 30
      n1=n
      m1=m
      mmn=m-n
      xn=dble(n) 
      xn2=xn*xn
      xm=dble(m) 
      xm2=xm*xm
      xmn2=xm2*xn2
      xm2mn2=xm2-xn2
      gnm= dble(xm2mn2)/dble(xmn2) 
      if(mmn.le.3.and.n.le.4)xknm=xknmtb(n,mmn)
      if(mmn.gt.3.or.n.gt.4)xknm=5.5e-5/gnm*xmn2/(1.+.13/float(mmn))
      y1num=320.
      if(m.eq.2)y1num=550.
      if(m.eq.3)y1num=380.
      y1wht=1.e13
      if(mmn.le.3)y1wht=1.e14
      if(mmn.le.2.and.n.le.2)y1wht=y1wtm(n,mmn)
      freqnm=rydh*gnm
      dbeta=2.99792458e18/freqnm**2/xknm
      wavenm=2.99792458e18/freqnm
      c1con=xknm/wavenm*gnm*dble(xm2mn2) 
      c2con=(xknm/wavenm)**2
!      radamp=1.389e9/xm**4.53/(1.+5./xm2/xm)
!      if(n.ne.1)radamp=radamp+1.389e9/xn**4.53/(1.+5./xn2/xn)
      radamp=asum(n)+asum(m)
      if(n.eq.1)radamp=asumlyman(m)
      radamp=radamp/12.5664
      radamp=radamp/freqnm
!      resont=hfnm(1,m)/xm/(1.-1./xm2)
!      if(n.ne.1)resont=resont+hfnm(1,n)/xn/(1.-1./xn2)
      resont=hfnm(1,m)/(1.-1./dble(xm2))
      if(n.ne.1)resont=resont+hfnm(1,n)/(1.-1./xn2)
!      fudge to baschek*2
      resont=resont*3.579e-24*0.5773/gnm
      vdw=4.45e-26/gnm*(dble(xm2)*(7.*dble(xm2)+5.))**.4
!     guess that h2 is twice as strong as he in txnxn
      hwvdw=vdw*t3nhe(j)+2.*vdw*t3nh2(j)
      hwrad=radamp
      stark=1.6678e-18*freqnm*xknm
!     fine structure components
!
!      if(n.gt.4)then
      if(n.gt.4.or.m.gt.10)then
      ifins=1
      finest(1)=0.
      finswt(1)=1.
      go to 30
      endif
!
      if(mmn.eq.1)go to 22
!     use m.eq.inf structure
      ifins=lncomp(n)
      do 21 i=1,ifins
      finest(i)=stcomp(i,n)*1.e7
      finswt(i)=stcpwt(i,n)/xn2
   21 continue
      go to 30
!     for alpha lines
   22 ifins=lnghal(n)
      ipos=istal(n)
      do 23 i=1,ifins
      k=ipos-1+i
      finest(i)=stalph(k)*1.e7
      finswt(i)=stwtal(k)/xn2/3.
   23 continue
!     now do this depth
   30 wl=wl0*10.+delw*10.
      freq=2.99792458e18/wl
      freq0=2.99792458e17/wl0
!      del=abs(freq-freqnm)
      del=abs(freq-freq0)
!     wl in nm
      wl=wl/10.
!     these half-widths are really dnu/nu
      hwstk=stark*fo(j)
      hwvdw=vdw*t3nhe(j)+2.*vdw*t3nh2(j)
      hwrad=radamp
!     xnfph(j,1)*2 is the number in the ground state
      hwres=resont*xnfph(j,1)*2.
      hwlor=hwres+hwvdw+hwrad
      hwdop=dopph(j)
!     specify largest half width in case of core calc
!     nwid=1, doppler  =2, lorentz  =3, stark
      nwid=1
      if(hwdop.ge.hwstk.and.hwdop.ge.hwlor)go to 31
      nwid=2
      if(hwlor.ge.hwstk)go to 31
      nwid=3
   31 hfwid=freqnm*max(hwdop,hwlor,hwstk)
!     sets flag if in a line core
!     hprofl=0.
      hprof4p=0.
      ifcore=0
      if(abs(del).le.hfwid)ifcore=1
!      dop=freqnm*hwdop
      dop=freq0*hwdop
      if(ifcore.eq.1)go to (32,40,50),nwid
!     do doppler
!     put fine structure in doppler core
   32 do 33 i=1,ifins
!      d=abs(freq-freqnm-finest(i))/dop
      d=abs(freq-freq0-finest(i))/dop
!     same normalization as voigt function
!     change in dfsynthe for speed
      if(d.le.7.)hprof4p=hprof4p+h0tab(int(200.*d+1.5))*finswt(i)
   33 continue
      if(ifcore.eq.1)return
!     do lorentz
   40 if(n.ne.1)go to 48
      if(m.ne.2)go to 48
!     lyman alpha
!     near center
!     modify old resonance broadening to match at 4000 cm-1
      hwres=hwres*4.
      hwlor=hwres+hwvdw+hwrad
      hhw=freqnm*hwlor
!     if (lambda<1277a)
      if(freq.gt.(82259.105-4000.)*2.99792458e10)then
      hprofres=hwres*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
      go to 44
      endif
!     only far red wing
!     data from n.f. allard, march 1997
!     insert lyman alpha cutoff a la n.f. allard and d. koester, a&a, 258,
!     464-468. 1992.
      cutoff=0.
      if(freq.lt.50000.*2.99792458e10)go to 43
!     tabulated at 200 cm-1 spacing
      spacing=200.*2.99792458e10
      freq22000=(82259.105-22000.)*2.99792458e10
      if(freq.lt.freq22000)then
      cutoff=(cutoffh2(2)-cutoffh2(1))/spacing*(freq-freq22000)+ & 
     &      cutoffh2(1)
      else
      icut=(freq-freq22000)/spacing
      icut=min(icut,89)
      cutfreq=icut*spacing+freq22000
      cutoff=(cutoffh2(icut+2)-cutoffh2(icut+1))/spacing* & 
     &  (freq-cutfreq)+cutoffh2(icut+1)
      endif
      cutoff=(10.**(cutoff-14.))*xnfph(j,1)*2./2.99792458e10
   43 hprofres=cutoff*1.77245*dop
!
   44 hprofrad=0.
!     rayleighscattering except near doppler core 
!     crossover from absorption to rayleigh scattering in hrayop
!     lambda>/=1217a
      if(freq.gt.2.4190611e15.and.freq.lt..77*3.28805e15) & 
     & hprofrad=hwrad*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
!     correction to lorentz profile   aller p.164   not used

      hprofvdw=hwvdw*freqnm/3.14159/(del**2+hhw**2)*1.77245*dop
!     van der waal cutoff for he and for h2
!     guess both 60000 cm-1=1.8e15 hz
      if(freq.lt.1.8e15)hprofvdw=0.
      hproflor=hprofres+hprofrad+hprofvdw
      hprof4p=hprof4p+hproflor
      if(ifcore.eq.1)return
      go to 50
!
!     not lyman alpha
   48 hhw=freqnm*hwlor
      top=hhw
      if(n.ne.1.or.m.gt.5)go to 49

!     lyman beta
      if(m.eq.3.and.freq.gt..885*3.288051e15.and. & 
     & freq.lt..890*3.288051e15)top=hhw-freqnm*hwrad
!     lyman gamma
      if(m.eq.4.and.freq.gt..936*3.288051e15.and. & 
     & freq.lt..938*3.288051e15)top=hhw-freqnm*hwrad
!     lyman delta
      if(m.eq.5.and.freq.gt..959*3.288051e15.and. & 
     & freq.lt..961*3.288051e15)top=hhw-freqnm*hwrad
   49 hproflor=top/3.14159/(del**2+hhw**2)*1.77245*dop
      hprof4p=hprof4p+hproflor
      if(ifcore.eq.1)return
!     do stark
   50 wty1=1./(1.+xne(j)/y1wht)
      y1scal=y1num*y1s(j)*wty1+y1b(j)*(1.-wty1)
      c1=c1d(j)*c1con*y1scal
      c2=c2d(j)*c2con
      g1=6.77*sqrt(c1)
      gnot=g1*max(0.0d0,.2114d0+log(sqrt(c2)/c1))*(1.-gcon1(j)-gcon2(j))
      beta=abs(del)/fo(j)*dbeta
      y1=c1*beta
      y2=c2*beta**2
      gam=gnot
!     if(y2.le..001)go to 51
      if(y2.le.1.e-4.and.y1.le.1.e-5)go to 51
!     change for dfsynthe for speed
      gam=g1*(.5*extab(nint(100.*min(80.0d0,y1)))+ & 
     &  vcse1f(y1)-.5*vcse1f(y2))* & 
     &  (1.-gcon1(j)/(1.+(90.*y1)**3)-gcon2(j)/(1.+2000.*y1))

      if(gam.le.1.e-20)gam=0.
   51 prqs=sofbet(beta,pp(j),n,m)
      if(m.gt.2)go to 53
!     assume quasistatic  profile is half protons, half electrons
      prqs=prqs*.5
      cutoff=0.

!     lyman alpha quasi h2+ cutoff
!     data from n.f. allard, march 1997
      if(freq.lt.(82259.105-20000.)*2.99792458e10)go to 53
      if(freq.gt.(82259.105-4000.)*2.99792458e10)go to 52
!     tabulated at 100 cm-1 spacing
      freq15000=(82259.105-15000.)*2.99792458e10
      spacing=100.*2.99792458e10
      if(freq.lt.freq15000)then
      cutoff=(cutoffh2plus(2)-cutoffh2plus(1))/spacing* & 
     &(freq-freq15000)+cutoffh2plus(1)
      else
      icut=(freq-freq15000)/spacing
      icut=min(icut,109)
      cutfreq=icut*spacing+freq15000
      cutoff=(cutoffh2plus(icut+2)-cutoffh2plus(icut+1))/ & 
     & spacing*(freq-cutfreq)+cutoffh2plus(icut+1)
      endif

!     xnfph(j,2)=xnfh(j,2)
      cutoff=(10.**(cutoff-14.))/2.99792458e10*xnfph(j,2)
      hprof4p=hprof4p+cutoff*1.77245*dop
      go to 53
   52 beta4000=4000.*2.99792458e10/fo(j)*dbeta
      prqsp4000=sofbet(beta4000,pp(j),n,m)*.5/fo(j)*dbeta
      cutoff4000=(10.**(-11.07-14.))/2.99792458e10*xnfph(j,2)
      hprof4p=hprof4p+cutoff4000/prqsp4000*prqs/fo(j)*dbeta*1.77245*dop
   53 f=0.
      if(gam.gt.0.)f=gam/3.14159/(gam**2+beta**2)
      p1=(.9*y1)**2
      fns=(p1+.03*sqrt(y1))/(p1+1.)
!     same normalization as voigt function
      hprof4p=hprof4p+(prqs*(1.+fns)+f)/fo(j)*dbeta*1.77245*dop
      return
      end function 


double precision  function hfnm(n,m)
!     calculates hydrogen oscillator strengths
  use types

  implicit none
      
  integer(kind=4) nstr, mstr, n, m
  real(kind=dp) xn, xm, xmn
  real(kind=dp) ginf, xmn12, wt, fnm, gca, fk, fkn, wtc 

  nstr=0
  mstr=0
  hfnm=0.0d0 

  if( m .gt. n) then 

    if(n .ne. nstr) then 
      xn=dble(n) 
      ginf=.2027/xn**.71
      gca=.124/xn
      fkn=xn*1.9603
      wtc=.45-2.4/xn**3*(xn-1.)
      nstr=n
!
      xm=dble(m) 
      xmn=dble(m-n) 
      fk=fkn*(xm/(xmn*(xm+xn)))**3
      xmn12=xmn**1.2
      wt=(xmn12-1.)/(xmn12+wtc)
      fnm=fk*(1.-wt*ginf-(.222+gca/xm)*(1.-wt))
      mstr=m
      hfnm=fnm

     else   

    
      if(m .ne. mstr) then 
         xm=dble(m) 
         xmn=dble(m-n) 
         fk=fkn*(xm/(xmn*(xm+xn)))**3
         xmn12=xmn**1.2
         wt=(xmn12-1.)/(xmn12+wtc)
         fnm=fk*(1.-wt*ginf-(.222+gca/xm)*(1.-wt))
         mstr=m
         hfnm=fnm
      else 

         hfnm=fnm
 
      end if      

     end if 
  end if
  return

end function 

!-----------------------------------------------------------
double precision function vcse1f(x)
!     rough, but arranged to be fast.  x.ge.0
    use types
    use dfcomm
  
    implicit none
!---- internal

    
    real(kind=dp) x, xn 
    vcse1f=0.

    if( x .gt. 0.) then 
      if(x .gt. .01) then 
        if(x .le. 1.) then 

          vcse1f=-log(x)-.57721566+x*(.99999193+x*(-.24991055+x*(.05519968+ & 
     &     x*(-.00976004+x*.00107857))))
        else 
          if(x.le.30.) then 
            vcse1f=(x*(x+2.334733)+.25062)/(x*(x+3.330657)+1.681534)/x*exp(-x)
          end if 
        end if 

      else 
        vcse1f=-log(x)-.577215+x
      end if
     end if 

end function 




  double precision  function sofbet(b,p,n,m)
!     generates s(beta,p) for hydrogen lines.  the alpha and beta lines
!     of the first three series are explicitly included and the h18
!     profile is used for the rest.
!
!     these profiles have been renormalized to full oscillator strength
!
!     storage for corrections (p,beta,ind),(p,ind),(p,ind)

   use types
   use dfcomm

   implicit none

   integer(kind=4) msave, nsave, m, n,  im, ip, index, j, jm, jp 
   real(kind=dp) p, psave
   real(kind=dp) b,  indexsave, wtpp, wtpm, cc, dd, b2, sb, wtbp 
   real(kind=dp) corr, wtbm, pr1, pr2, wt


      real(kind=dp) propbm(5,15,7),c(5,7),d(5,7)
      real(kind=dp) pp(5),beta(15),probeta(15)
      real(kind=dp) prob1(75),prob2(75),prob3(75),prob4(75),prob5(75)
      real(kind=dp) prob6(75),prob7(75)

      real(kind=dp) c1(5),c2(5),c3(5),c4(5),c5(5),c6(5),c7(5)
      real(kind=dp) d1(5),d2(5),d3(5),d4(5),d5(5),d6(5),d7(5)

      equivalence (propbm(1,1,1),prob1(1)),(propbm(1,1,2),prob2(1))
      equivalence (propbm(1,1,3),prob3(1)),(propbm(1,1,4),prob4(1))
      equivalence (propbm(1,1,5),prob5(1)),(propbm(1,1,6),prob6(1))
      equivalence (propbm(1,1,7),prob7(1))
      equivalence (c(1,1),c1(1)),(c(1,2),c2(1)),(c(1,3),c3(1)),(c(1,4),c4(1))
      equivalence (c(1,5),c5(1)),(c(1,6),c6(1)),(c(1,7),c7(1))
      equivalence (d(1,1),d1(1)),(d(1,2),d2(1)),(d(1,3),d3(1)),(d(1,4),d4(1))
      equivalence (d(1,5),d5(1)),(d(1,6),d6(1)),(d(1,7),d7(1))

!----------------------------------------------------------------------
!     lyman alpha
      prob1=(/ & 
     &-.980,-.967,-.948,-.918,-.873,-.968,-.949,-.921,-.879,-.821, & 
     &-.950,-.922,-.883,-.830,-.764,-.922,-.881,-.830,-.770,-.706, & 
     &-.877,-.823,-.763,-.706,-.660,-.806,-.741,-.682,-.640,-.625, & 
     &-.691,-.628,-.588,-.577,-.599,-.511,-.482,-.484,-.514,-.568, & 
     &-.265,-.318,-.382,-.455,-.531,-.013,-.167,-.292,-.394,-.478, & 
     & .166,-.056,-.216,-.332,-.415, .251, .035,-.122,-.237,-.320, & 
     & .221, .059,-.068,-.168,-.247, .160, .055,-.037,-.118,-.189, & 
     & .110, .043,-.022,-.085,-.147/) 

      c1=(/-18.396, 84.674,-96.273,  3.927, 55.191/) 
      d1= (/ 11.801,  9.079,  -.651,-11.071,-26.545/) 


!--- restore in propbm, c and d

!     do l = 1, 15
!      do i = 1, 5
!       propbm(i,l,1) = prob1((l-1)*5 +i)
!      end do 
!     end do 

!     lyman beta
      prob2=(/ & 
     &-.242, .060, .379, .671, .894, .022, .314, .569, .746, .818, & 
     & .273, .473, .605, .651, .607, .432, .484, .489, .442, .343, & 
     & .434, .366, .294, .204, .091, .304, .184, .079,-.025,-.135, &
     & .167, .035,-.082,-.189,-.290, .085,-.061,-.183,-.287,-.374, & 
     & .032,-.127,-.249,-.344,-.418,-.024,-.167,-.275,-.357,-.420, & 
     &-.061,-.170,-.257,-.327,-.384,-.047,-.124,-.192,-.252,-.306, & 
     &-.043,-.092,-.142,-.190,-.238,-.038,-.070,-.107,-.146,-.187, & 
     &-.030,-.049,-.075,-.106,-.140/) 
     

       c2 =(/ 95.740, 18.489, 14.902, 24.466, 42.456/) 
      d2 = (/ -6.665, -7.136,-10.605,-15.882,-23.632/) 


!     balmer alpha
      prob3=(/ & 
     &-.484,-.336,-.206,-.111,-.058,-.364,-.264,-.192,-.154,-.144, & 
     &-.299,-.268,-.250,-.244,-.246,-.319,-.333,-.337,-.336,-.337, & 
     &-.397,-.414,-.415,-.413,-.420,-.456,-.455,-.451,-.456,-.478, & 
     &-.446,-.441,-.446,-.469,-.512,-.358,-.381,-.415,-.463,-.522, & 
     &-.214,-.288,-.360,-.432,-.503,-.063,-.196,-.304,-.394,-.468, & 
     & .063,-.108,-.237,-.334,-.409, .151,-.019,-.148,-.245,-.319, & 
     & .149, .016,-.091,-.177,-.246, .115, .023,-.056,-.126,-.189, & 
     & .078, .021,-.036,-.091,-.145/) 



       c3=(/-25.088,145.882,-50.165,  7.902, 51.003/)
       d3=(/  7.872,  5.592, -2.716,-12.180,-25.661/) 


!     balmer beta
     prob4=(/ & 
     &-.082, .163, .417, .649, .829, .096, .316, .515, .660, .729, & 
     & .242, .393, .505, .556, .534, .320, .373, .394, .369, .290, & 
     & .308, .274, .226, .152, .048, .232, .141, .052,-.046,-.154, & 
     & .148, .020,-.094,-.200,-.299, .083,-.070,-.195,-.299,-.385, &
     & .031,-.130,-.253,-.348,-.422,-.023,-.167,-.276,-.359,-.423, & 
     &-.053,-.165,-.254,-.326,-.384,-.038,-.119,-.190,-.251,-.306, & 
     &-.034,-.088,-.140,-.190,-.239,-.032,-.066,-.103,-.144,-.186, & 
     &-.027,-.048,-.075,-.106,-.142/) 


       c4=(/ 93.783, 10.066,  9.224, 20.685, 40.136/) 
       d4=(/ -5.918, -6.501,-10.130,-15.588,-23.570/) 


!     paschen alpha
      prob5=(/ & 
     &-.819,-.759,-.689,-.612,-.529,-.770,-.707,-.638,-.567,-.498, & 
     &-.721,-.659,-.595,-.537,-.488,-.671,-.617,-.566,-.524,-.497, & 
     &-.622,-.582,-.547,-.523,-.516,-.570,-.545,-.526,-.521,-.537, & 
     &-.503,-.495,-.496,-.514,-.551,-.397,-.418,-.448,-.492,-.547, & 
     &-.246,-.315,-.384,-.453,-.522,-.080,-.210,-.316,-.406,-.481, & 
     & .068,-.107,-.239,-.340,-.418, .177,-.006,-.143,-.246,-.324, & 
     & .184, .035,-.082,-.174,-.249, .146, .042,-.046,-.123,-.190, & 
     & .103, .036,-.027,-.088,-.146/) 


      c5=(/-19.819, 94.981,-79.606,  3.159, 52.106/) 
      d5=(/ 10.938,  8.028, -1.267,-11.375,-26.047/) 
!     paschen beta

      prob6=(/ & 
     &-.073, .169, .415, .636, .809, .102, .311, .499, .639, .710, &
     & .232, .372, .479, .531, .514, .294, .349, .374, .354, .279, & 
     & .278, .253, .212, .142, .040, .215, .130, .044,-.051,-.158, & 
     & .141, .015,-.097,-.202,-.300, .080,-.072,-.196,-.299,-.385, & 
     & .029,-.130,-.252,-.347,-.421,-.022,-.166,-.275,-.359,-.423, &
     &-.050,-.164,-.253,-.325,-.384,-.035,-.118,-.189,-.252,-.306, & 
     &-.032,-.087,-.139,-.190,-.240,-.029,-.064,-.102,-.143,-.185, & 
     &-.025,-.046,-.074,-.106,-.142/) 



     c6=(/111.107, 11.910,  9.857, 21.371, 41.006/) 
     d6=(/ -5.899, -6.381,-10.044,-15.574,-23.644/) 


!     balmer 18
      prob7=(/ & 
     & .005, .128, .260, .389, .504, .004, .109, .220, .318, .389, & 
     &-.007, .079, .162, .222, .244,-.018, .041, .089, .106, .080, & 
     &-.026,-.003, .003,-.023,-.086,-.025,-.048,-.087,-.148,-.234, & 
     &-.008,-.085,-.165,-.251,-.343, .018,-.111,-.223,-.321,-.407, & 
     & .032,-.130,-.255,-.354,-.431, .014,-.148,-.269,-.359,-.427, & 
     &-.005,-.140,-.243,-.323,-.386, .005,-.095,-.178,-.248,-.307, & 
     &-.002,-.068,-.129,-.187,-.241,-.007,-.049,-.094,-.139,-.186, & 
     &-.010,-.036,-.067,-.103,-.143/) 


     c7=(/511.318,  1.532,  4.044, 19.266, 41.812/) 
     d7=(/ -6.070, -4.528, -8.759,-14.984,-23.956/) 
     pp=(/0.,.2,.4,.6,.8/) 
     beta=(/1.,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.310,7.943 & 
     &,10.,12.59,15.85,19.95,25.12/) 

!




      msave= 0 
      nsave= 0
      psave= 0.0d0 
      indexsave = 0.0d0 



      if(b.lt.500.)go to 1
      sofbet=(1.5/sqrt(b)+27./b**2)/b**2
      return
    1 if(m.eq.msave.and.n.eq.nsave.and.p.eq.psave)go to 5
      msave=m
      nsave=n
      psave=p
      index=7
      if(n.le.3.and.m-n.le.2)index=n+m-2
      if(index.eq. int(indexsave))go to 5
      im=min(int(5.*p)+1,4)
      ip=im+1
      wtpp=5.*(p-pp(im))
      wtpm=1.-wtpp
      cc=c(ip,index)*wtpp+c(im,index)*wtpm
      dd=d(ip,index)*wtpp+d(im,index)*wtpm
      do  j=1,15
       probeta(j)=propbm(ip,j,index)*wtpp+propbm(im,j,index)*wtpm
      end do

      indexsave=dble(index) 
    5 b2=b*b
      sb=sqrt(b)
      if(b.lt.25.12)go to 7
      corr=1.+dd/(cc+b*sb)
      sofbet=(1.5/sb+27./b2)/b2*corr
      return
    7 do j=2,15
      if(b.le.beta(j))go to 20
      end do
   20 jm=j-1
      jp=j
      wtbp=(b-beta(jm))/(beta(jp)-beta(jm))
      wtbm=1.-wtbp
      corr=1.+probeta(jp)*wtbp+probeta(jm)*wtbm
!      get inner approximate profile
      pr1=0.
      pr2=0.
      wt=max(min(.5d0*(10.0d0-b),1.0d0),0.0d0)
      if(b.le.10.)pr1=8./(83.+(2.+.95*b2)*b)
      if(b.ge.8.)pr2=(1.5/sb+27./b2)/b2
      sofbet=(pr1*wt+pr2*(1.-wt))*corr
      return
      end function 








END MODULE HPROF
