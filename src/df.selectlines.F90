  subroutine selectlines(n12,n19, ntemp, workdir )
  
! This subroutine opens linelists: 
! -  low-lines
! -  high-lines (T > 20 000K )
! -  TIO -lines (T < 5 000K ) 
! -  H2O - lines (T< 5 000K ) 
! -  Diatomic lines (T< 10 000K )
! -  NLTE lines
! 
!    if xnfdopmx = (1/rho  N_j (element)/ U_j(element) 1/delta V) / (cutoff * k_c(min))
!    is > 1
!    and  if two other conditions apply (see DFSYNTHE: how to use it by Castelli, http://sait.oat.ts.astro.it/MSAIS/8/PDF/34.pdf, the line is considered.  


    use types
    use dfcomm
    
    implicit none

    include 'common.teffbl'
    include 'common.odfnlte'
    include 'common.tempbl'

    integer, intent(in) :: ntemp
    character(len=256),intent(in) ::workdir
    integer, intent(out) :: n12, n19 

    integer(kind=4) itape12(6,10000)
    real(kind=4) tape12(6,10000),  tablog(32768) 

    real(kind=dp) wl, wlvac   

    integer(kind=2) ielion,ielo,igflog,igr,igs,igw
    integer(kind=4) iwltab(nufreq3)
    integer(kind=4) itape11(4,2000)
    integer(kind=4) gnnnggg(11)
    integer(kind=2) i2tape11(8,2000)

! ------------------------------
    real(kind=dp) hkt100, hckt100 
    real(kind=dp) wavetab(nufreq3), tabbolt(32768), xnfdopmax(377,594) 

    integer nbuff, i, j, nt,  n22, n32, n42, n44, n52, n62, l, line 
    integer(kind=4)  nu, numnu,iwl,   l11,  nelion, nlines,  ixwl 
    integer(kind=4) iwlstart, iwlstop, nzero 

!---- variables that where implicit in the old version -----------------------
    real(kind=4) wlvac4,freq4 , ratiolg4, wlbeg04 , frq4pi 
    real(kind=4) elo
    real(kind=4) gf, congf, boltstim
!
!----------------------------------------------------------------------------------
    equivalence (gnnnggg(1), gf) 
    equivalence (tape12(1,1),itape12(1,1))
    equivalence (itape11(1,1),i2tape11(1,1)) 

!------------------------------------------------------------------
!  set a iwlstart and iwlstop depending on the lowest and highest wavelength 
!  on the bin-grid
!--------------------------------------------------------------------------
    nzero=log(wlbeg)/log(1.d0+1.d0/500000.d0) 
    iwlstart = int(log(wavebig(1)-1.5d0)/ratiolg)-nzero +1 
    iwlstop  = int(log(wavebig(nsizebig+1)+1.5d0)/ratiolg)-nzero+1
 
    if (iwlstart .lt. 1) iwlstart = 1
    if (iwlstop .gt. length) iwlstop = length 
    print*,'iwlstart and stop = ', iwlstart, iwlstop

!-------------------------------------------------------------
    do nu=1, 377
      do j= 1,6
        do l= 1,99
         i=(l-1)*6+j 
         xnfdopmax(nu,i) = dopmax(nu,j,l)
        end do
       end do
    end do  

!-------------------------------------------------------------
      
  do  i=1,32768
      tablog(i)=10.**(float(i-16384)*.001)
  end do    
!--------------------------------------------------------------
        
  ratiolg4=ratiolg
  wlbeg04=wlbeg0
  if(dexp(ixwlbeg*ratiolg).lt.wlbeg)ixwlbeg=ixwlbeg+1

  hkt100=hkt(ntemp)*100.
  hckt100=hckt(ntemp)*100.
  
 
  do  i=1,nufreq3
      wavetab(i)=abs(wledge(i))
      ixwl=dlog(wavetab(i))/ratiolg
      ixwl=ixwl+1
      iwltab(i)=ixwl-ixwlbeg+1
  end do

  do  i=1,32768
     tabbolt(i)=exp(-(i-1)*hckt(ntemp))
  end do

  n12=0
  n22=0
  n32=0
  n42=0
  n44=0
  n52=0
  n62=0
  nu=1
  nbuff=0
!---------------------------------------------------
!     lowlines
!-----------------------------------------------------------
#ifdef GFORT
  open(unit=11, status='old', file=trim(workdir)//'/fort.11', form='unformatted',action='read')
#else
  open(unit=11, status='old', file=trim(workdir)//'/fort.11', form='unformatted',shared,readonly, recordtype='fixed',blocksize=32000,recl=8000)
#endif

  rewind 12
  l=0
  l11=2000


!  do  line=1,31108567
  do  line=1,50000000

    l11=l11+1
   ! read line and all information 

   if(l11.eq.2001)then
     read(11, end = 980 )itape11
     l11=1
   endif
   ! get position of the line on the high-res op grid
    iwl=itape11(1,l11)

    if (iwl .eq. 0 ) go to 980

!--- check if line is within requiered intervall

   if ((iwl .ge. iwlstart) .and.(iwl .le. iwlstop)) then

    do while (iwl .ge. iwltab(nu))
      nu=nu+1
    end do

    nelion=i2tape11(3,l11)

    if(xnfdopmax(nu,nelion) .ge. 1.d0 ) then
       
      congf=dble(.026538/1.77245*tablog(i2tape11(5,l11)))
      if (i2tape11(5,l11) .lt. 1) congf = 0.0

      if(congf*xnfdopmax(nu,nelion) .ge. 1.d0) then

        elo=tablog(i2tape11(4,l11))
        congf=dble(congf*exp(-elo*hckt(ntemp)))

        if(congf*xnfdopmax(nu,nelion).ge. 1.0d0 ) then 
        
          if(iwl.ne.nbuff) wlvac4=exp(iwl*ratiolg4)*wlbeg04

          nbuff=iwl
          l=l+1
          itape12(1,l)=nbuff
          freq4=2.99792458e17/wlvac4

          tape12(2,l)=(congf/freq4*(1.-exp(-freq4*hkt(ntemp))))
          itape12(3,l)=nelion
          frq4pi=freq4*12.5664

          tape12(4,l)=tablog(i2tape11(6,l11))/(frq4pi)
          tape12(5,l)=tablog(i2tape11(7,l11))/(frq4pi) 
          tape12(6,l)=tablog(i2tape11(8,l11))/(frq4pi) 

          n12=n12+1

          if(l.eq.10000)then
            write(12)tape12
            l=0
          endif
       end if
      end if
    end if
! check end if line is in intervall
   end if

  end do



! jumps out of previous loop if the line list is finished, either no more to read or line has iwl = 0 

980 continue      


  write(6,8881)n12
 8881 format(i10,' lines from lowlines')
  close(unit=11)



#ifdef ORIGL
!---------------------------------------------------------------
!     hilines
!--------------------------------------------------------------
      
  if (tsave(ntemp) .ge. 20000.0d0 ) then 
#ifdef GFORT
   open(unit=21,status='old',form='unformatted', file=trim(workdir)//'/fort.21', action='read')
#else
   open(unit=21,status='old',form='unformatted', file=trim(workdir)//'/fort.21', shared,readonly, recordtype='fixed',blocksize=32000,recl=8000)
#endif 
   nu=1
   nbuff=0
   l11=2000

   do line=1,10057574 ! 234077    

    l11=l11+1

    if(l11.eq.2001)then
      read(21)itape11
      l11=1
    endif

    iwl=itape11(1,l11)
    if ((iwl .ge. iwlstart) .and. (iwl .le. iwlstop)) then

       do while (iwl .ge. iwltab(nu))
          nu = nu +1
       end do
       
      nelion=i2tape11(3,l11)

      if(xnfdopmax(nu,nelion) .ge. 1.0d0 )  then 
         congf=dble(.026538/1.77245*tablog(i2tape11(5,l11)))

        if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then 
          elo=tablog(i2tape11(4,l11))
          congf=congf*dble(exp(-elo*hckt(ntemp)))

          if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
            if(iwl.ne.nbuff)wlvac4=exp(iwl*ratiolg4)*wlbeg04
             nbuff=iwl
             l=l+1
             itape12(1,l)=nbuff
             freq4=2.99792458d17/wlvac4
             tape12(2,l)=(congf/freq4*(1.-exp(-freq4*hkt(ntemp))))
             itape12(3,l)=nelion
             frq4pi=freq4*12.5664
             tape12(4,l)=tablog(i2tape11(6,l11))/(frq4pi) 
             tape12(5,l)=tablog(i2tape11(7,l11))/(frq4pi) 
             tape12(6,l)=tablog(i2tape11(8,l11))/(frq4pi)
             n22=n22+1

             if(l.eq.10000)then
               write(12)tape12
               l=0
             endif
          endif
        endif
      endif
   endif
   end do
 
   write(6,8882)n22
 8882 format(i10,' lines from hilines')
   close(unit=21)

  endif
#endif 

   
!------------------------------------------------------------------
!     diatomics
!------------------------------------------------------------------

  if (tsave(ntemp) .lt. 10000.0d0 ) then

#ifdef GFORT
   open(unit=31,status='old',form='unformatted',file=trim(workdir)//'/fort.31',action='read') 
#else
   open(unit=31,status='old',form='unformatted',file=trim(workdir)//'/fort.31',shared,readonly, recordtype='fixed',blocksize=32000,recl=8000)
#endif 
   nu=1
   nbuff=0
   l11=2000

   do line=1,7771848
      l11=l11+1
      if(l11.eq.2001)then
        read(31)itape11
        l11=1
      endif
      iwl=itape11(1,l11)
!  check if line in in the wavelength interval !

   if ((iwl .ge. iwlstart) .and. (iwl .le. iwlstop)) then 
 
      do while (iwl .ge. iwltab(nu)) 
         nu=nu+1
      end do
      
      nelion=i2tape11(3,l11)
     

      if(xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
         congf=dble(.026538/1.77245*tablog(i2tape11(5,l11)))

         if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0) then
 
           elo=tablog(i2tape11(4,l11))
           congf=congf*dble(exp(-elo*hckt(ntemp)))

           if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
 
              if(iwl.ne.nbuff)wlvac4=exp(iwl*ratiolg4)*wlbeg04

              nbuff=iwl
              l=l+1
              itape12(1,l)=nbuff
              freq4=2.99792458d17/wlvac4
              tape12(2,l)=(congf/freq4*(1.-exp(-freq4*hkt(ntemp))))
              itape12(3,l)=nelion
              frq4pi=freq4*12.5664
              tape12(4,l)=tablog(i2tape11(6,l11))/(frq4pi)
              tape12(5,l)=tablog(i2tape11(7,l11))/(frq4pi) 
              tape12(6,l)=tablog(i2tape11(8,l11))/(frq4pi)
              n32=n32+1

              if(l.eq.10000)then
                 write(12)tape12
                 l=0
              endif
           endif
        endif
     endif
!  check interval
   end if 
  end do
    write(6,8883)n32
 8883 format(i10,' lines from diatomics')
    close(unit=31)
 endif
   
!------------------------------------------------------------------
!     tio
!-----------------------------------------------------------------
 if (tsave(ntemp) .lt. 5000.0d0 ) then
#ifdef GFORT
    open(unit=41,status='old',form='unformatted', file=trim(workdir)//'/fort.41', action='read' )
#else
    open(unit=41,status='old',form='unformatted', file=trim(workdir)//'/fort.41', shared,readonly, recordtype='fixed',blocksize=32000,recl=8000)
#endif
    nu=1
    nbuff=0
    l11=2000

    do  line=1,36979284
      l11=l11+1
      if(l11.eq.2001)then
        read(41)itape11
        l11=1
      endif
   
      iwl=itape11(1,l11)

   if ((iwl .ge. iwlstart) .and. (iwl .le. iwlstop)) then  

      do while (iwl .ge. iwltab(nu))
         nu=nu+1
      end do
     
      nelion=i2tape11(3,l11)
      if(xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
        congf=dble(.026538/1.77245*tablog(i2tape11(5,l11)))
        if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
          elo=tablog(i2tape11(4,l11))
          congf=congf*dble(exp(-elo*hckt(ntemp)))
          if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
            if(iwl.ne.nbuff)wlvac4=exp(iwl*ratiolg4)*wlbeg04
              nbuff=iwl
              l=l+1
              itape12(1,l)=nbuff
              freq4=2.99792458d17/wlvac4
              tape12(2,l)=(congf/freq4*(1.-exp(-freq4*hkt(ntemp))))
              itape12(3,l)=nelion
              frq4pi=freq4*12.5664
              tape12(4,l)=tablog(i2tape11(6,l11))/(frq4pi) 
              tape12(5,l)=tablog(i2tape11(7,l11))/(frq4pi) 
              tape12(6,l)=tablog(i2tape11(8,l11))/(frq4pi) 
              n42=n42+1
              if(l.eq.10000)then
                 write(12)tape12
                 l=0
              endif
           endif
        endif
     endif     
!  check interval 
   end if 
   end do
   write(6,8884)n42
 8884 format(i10,' lines from tiolines')
   close(unit=41)
 endif
   
!------------------------------------------------------------
!
!     h2o
!-------------------------------------------------------------
      
 if(tsave(ntemp).lt.5000.) then
#ifdef GFORT
    open(unit=43,status='old',form='unformatted', file=trim(workdir)//'/fort.43', action='read')
#else
    open(unit=43,status='old',form='unformatted', file = trim(workdir)//'/fort.43',shared,readonly, recordtype='fixed',blocksize=32000,recl=8000)
#endif
    nu=1
    nbuff=0
    l11=2000
    do line=1,48999843
      l11=l11+1

      if(l11.eq.2001)then
        read(43)itape11
        l11=1
      endif

      iwl=itape11(1,l11)
! check if line is in interval
    if ((iwl .ge. iwlstart) .and. (iwl .le. iwlstop)) then  

      do while (iwl .ge. iwltab(nu)) 
         nu=nu+1
      end do
      
      nelion=i2tape11(3,l11)
      if(xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
        congf=dble(.026538/1.77245*tablog(i2tape11(5,l11)))
        if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then

          congf=congf*tabbolt(i2tape11(4,l11)+1)
          if(congf*xnfdopmax(nu,nelion) .ge. 1.0d0 ) then
             if(iwl.ne.nbuff)wlvac4=exp(iwl*ratiolg4)*wlbeg04
             nbuff=iwl
             l=l+1
             itape12(1,l)=nbuff
             freq4=2.99792458d17/wlvac4
             tape12(2,l)=(congf/freq4*(1.-exp(-freq4*hkt(ntemp))))
             itape12(3,l)=nelion
             frq4pi=freq4*12.5664
             tape12(4,l)=tablog(i2tape11(6,l11))/(frq4pi) 
             tape12(5,l)=tablog(i2tape11(7,l11))/(frq4pi)
             tape12(6,l)=tablog(i2tape11(8,l11))/(frq4pi) 
             n44=n44+1
             if(l.eq.10000)then
               write(12)tape12
               l=0
            endif
         endif
      endif
     endif
!  end interval check
    end if 

    enddo

    write(6,1884)n44
 1884 format(i10,' lines from h2olines')
    close(unit=43)
  end if
   
!------------------------------------------------------------------
!     nltelinesa
!------------------------------------------------------------------
#ifdef GFORT
   open(unit=51,status='old',form='unformatted', file=trim(workdir)//'/fort.51', action='read')
   open(unit=19,status='new',form='unformatted', file=trim(workdir)//'/fort.19')
#else
   open(unit=51,status='old',form='unformatted', file=trim(workdir)//'/fort.51', shared,readonly, recordtype='fixed',blocksize=32000,recl=16)
   open(unit=19,status='new',form='unformatted', file=trim(workdir)//'/fort.19', recordtype='fixed',blocksize=28000,recl=14)
#endif

   nu=1

   do line=1,999999

      read(51, end = 981)wl,elo, gnnnggg
       nelion = gnnnggg(4)
       nbuff = gnnnggg(11) 


      if ((nbuff .ge. iwlstart) .and. (nbuff .le. iwlstop)) then  

        do while (nbuff .ge. iwltab(nu))
          nu=nu+1
        end do
   
        if(xnfdopmax(nu,nelion) .ge. 1.0d0) then 
 
          boltstim=exp(-(elo)*hckt(ntemp))*(1.-exp(-1.d7/wl*hckt(ntemp)))
          if(boltstim*xnfdopmax(nu,nelion) .ge. 1.d0) then 
             gf = gf*boltstim
             write(19)wl,gnnnggg
 
             n52=n52+1
          endif
        endif

      end if 

   end do
981 continue   
! end nltelines
!------------------------------------------------------
  
   write(6,8885)n52
  8885 format(i10,' lines from nltelinesa')
   n19=n52
   close(unit=51)

 
   n12=n12+n22+n32+n42+n44+n62
   nlines=n12+n19
  write(6,8888)nlines
 8888 format(i10,' lines total')
  write(6,*) 'n12, n19, nTemp :', n12, n19, ntemp 
      
   l=l+1
   itape12(3,l)=0
   write(12)tape12
    
  end subroutine

