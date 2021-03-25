subroutine dfcalc_filter(iteff,j,iv,t,p,vt)
!!ADD COMMENTS HERE!!
  use types
  use dfcomm
  implicit none

  include 'common.ifopbl'


  integer,  intent(in) :: iteff, j, iv 
  real(kind=dp),  intent(in):: t,p, vt
!-------------------------------------------------------
! local variables
!-------------------------------------------------------
  integer*8 ::  ntot,  nbuff, idf, nbegin, nstop
  integer :: i, l, jj, l1, nu,  istep, isub
  integer*4 :: ibuffer(lenbuff), ibuffsort(223938)

  integer*4 :: nnsub, idfnum
  integer*8 :: isum


  real(kind=dp):: tlog, plog !careful the same are globale in common.tempbl that is not used here!
!--------------------------------------------------------
  real(kind=dp):: weights(nsubbin+1),sum1,err0,sumo,wt(nsubbin)
  real(kind=dp):: xdlam,xlambda0,xlambda,xden,xnum,wavg,xsum,sbsum,fsum,xlam,fn,sump,sumpb,fnb,dsum,wsum
  integer*4 :: nsum,j1,jj1,ii,ii1,idiff,jo,ifilter1,dcnt,xcnt
  integer*4 :: isum1,jsum,isum2,jsum2,isum3 
!-----------------------------------------------------------------------
     tlog=log10(t)
     plog=log10(p)
     ifilter1=1 

!!!i_AM=0: GM, i_AM=1: AM; i_AM=2: HM
!!changed it to ameanODF
   !i_AM=0
   !!if ameanODF is true, i_AM=1 
   !if(ameanODF .eq. 1) i_AM=1

   do nbuff=1, lenbuff ! length as far as I can see you do not need length here. 
     ibuffer(nbuff)=nint(log10(buffer(nbuff))*1000.)
     ! ibuffer(nbuff)=nint(log10(buffer(nbuff))*1000000.)
! debug !     write(97,*) 'nbuff, ibuffer', nbuff, buffer(nbuff)
   end do 

!  print*,'NEW idf loop'
!!!!!!!!!!!!!!!!start bin by bin computation!!!!!!!!!!!!!!
  !!high resolution lambda and lambdatilda grids
  lambda=0.
  lambdatilda=0.
  jo=0
  
  nnsub = nsubbin

  do idf=1, nsizelit+nsizebig 
    insteps = 0
    sumpb=0.

    weights(1)=0.
    do i=1,nnsub!+1 -- so you declared weights(nnsub+1), but binbound(nnsub), thus seg fault.. 
    weights(i+1)=binbound(idf,i) !bec binbound are only upper boundaries!!
!    print*,'weights=',weights(i)
    enddo

      
    nbegin=int(nbeg(idf))
    nstop=int(nend(idf))
    ntot=int(nstop-nbegin+1) 
   !print_if_needed
   !print*,'idf=',idf,'ntot=',ntot

    do i=nbegin,nstop
       ibuffsort(i-nbegin+1)=ibuffer(i)
!test
!      ibuffsort(i-nbegin+1)=nint(log10(buffer(i))*dlamtilda(idf,i-nbegin+1)*1000.)
      indexx(i-nbegin+1)=i
!   write(*,"(1(1x,i4),1(1x,i4),2(1x,e15.6))")idf,int(i),(buffer(i))
    end do  

    call bisort_filter(ibuffsort,ntot,ibuffsort(100000))

!   print*,'called bisort'
    do i=1,ntot
    jj1=indexx(i)
!   write(*,"(1(1x,i4),1(1x,i4),2(1x,e15.6))")idf,i,(buffer(jj1))
    end do  

      if(ifilter1==0) then
          sbwith(1,idf) = binbound(idf,1)
       do istep=2, nsubbin
          sbwith(istep,idf) = binbound(idf,istep)-binbound(idf,istep-1)
       enddo
      endif !

!      print*,'ntot=',ntot

     if(ifilter1==0) then

      do istep=1, nnsub
       nsteps(istep,idf)=0
       nsteps(istep+1,idf)=0
!       print*,'sbwith-orig=',sbwith(istep,idf)
!       sbwith(istep,idf)=0.
      enddo
      dsum=0.
      dcnt=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           print*,'idf=',idf,'nsizebig=',nsizebig,'nsizelit=',nsizelit
           if(idf.le.nsizebig) then 
               xlam=(wavebig(idf+1)-wavebig(idf)) !width of the interval
               xlambda0=wavebig(idf)
            if(wavebig(idf+1).gt.0) then
!! first end point 
             nsteps(1,idf)=1
!! subsequent steps
            nsum=0
         do j1=1,ntot !high resolution grid per bin
             jj1=indexx(j1)
             jj=indexx(j1)-nbegin+1 !! high resoultion original index
            
!            print*,'j1=',j1
            if(j1==1)xlambda=xlambda0
            if(j1>1) then
              xlambda=xlambda+dlambda(idf,jj) !non-equidistant high res lambda grid
            endif
            if(j1==ntot) xlambda=wavebig(idf+1)

            do isub=2,2

             if( xlambda.ge.(xlam*weights(isub-1)+xlambda0).and. &
              xlambda.le.(xlam*weights(isub)+xlambda0) )then

              nsum=nsum+1
              nsteps(isub,idf)=nsteps(isub,idf)+1
             endif
            enddo
            do isub=3,nnsub+1

              if( xlambda.gt.(xlam*weights(isub-1)+xlambda0).and. &
              xlambda.le.(xlam*weights(isub)+xlambda0) )then

               nsum=nsum+1
               nsteps(isub,idf)=nsteps(isub,idf)+1

              endif
            enddo

         enddo ! j=1,ntot loop 
!         print*,'big-bin: nsum=',nsum,'ntot',ntot
           ii=1
           !print_if_needed
           !print*,' nsteps-nsteps',ii,nsteps(ii,idf)
           ii=2
           !print_if_needed
           !print*,' nsteps-nsteps',ii,nsteps(ii,idf)
           do ii=3,nnsub+1!do not add the first boundary value 1 to end of 1st bin//start from 2nd bin//
!           print*,'weights',ii,weights(ii)
           if(weights(ii).ne.0.) then
           nsteps(ii,idf)=nsteps(ii-1,idf)+nsteps(ii,idf)
           endif
           if(weights(ii).eq.0.) then
           nsteps(ii,idf)=1
           endif
           
           !print_if_needed
           !print*,'nsteps-nsteps',ii,nsteps(ii,idf)
           enddo
         endif !wavbig gt 0
        endif ! idf le nsizebig 

        if(idf.gt.nsizebig) then
               xlam=(wavelit(idf-nsizebig+1)-wavelit(idf-nsizebig)) !width of the interval
               xlambda0=wavelit(idf-nsizebig)
            if(wavelit(idf-nsizebig+1).gt.0) then
!! first end point 
             nsteps(1,idf)=1
!! subsequent steps
            nsum=0
          do j1=1,ntot !high resolution grid per bin
           jj=indexx(j1)-nbegin+1 !! high resoultion original index

!            print*,'j1=',j1
            if(j1==1)xlambda=xlambda0
            if(j1>1) xlambda=xlambda+dlambda(idf,jj) !non-equidistant high res lambda grid
            if(j1==ntot) xlambda=wavelit(idf-nsizebig+1)

            do isub=2,2
             if( xlambda.ge.(xlam*weights(isub-1)+xlambda0).and. &
              xlambda.le.(xlam*weights(isub)+xlambda0) )then
                nsum=nsum+1
                nsteps(isub,idf)=nsteps(isub,idf)+1
              endif
            enddo
            do isub=3,nnsub+1
              if( xlambda.gt.(xlam*weights(isub-1)+xlambda0).and. &
               xlambda.le.(xlam*weights(isub)+xlambda0) )then
                nsum=nsum+1
                nsteps(isub,idf)=nsteps(isub,idf)+1
              endif
           enddo
         enddo ! j=1,ntot loop 

!         print*,'lit-bin: nsum=',nsum,'ntot',ntot
           ii=1
!           print*,'lit-bin: nsteps-nsteps',ii,nsteps(ii,idf)
           ii=2
!           print*,'lit-bin: nsteps-nsteps',ii,nsteps(ii,idf)
           do ii=3,nnsub+1 !do not add the first boundary value 1 to end of 1st bin//start from 2nd bin//
           if(weights(ii).ne.0.) then
           nsteps(ii,idf)=nsteps(ii-1,idf) &
                             +nsteps(ii,idf)
           endif
           if(weights(ii).eq.0.) then
           nsteps(ii,idf)=1
           endif
!           print*,'lit-bin: nsteps-nsteps',ii,nsteps(ii,idf)
           enddo
         endif !wavlit gt 0
        endif ! idf gt nsizebig 
     endif !!!ifilter1==0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(ifilter1==1) then

      do istep=1, nnsub
       nnsteps(istep,idf)=0
       nnsteps(istep+1,idf)=0
!       print*,'sbwith-orig=',sbwith(istep,idf)
       sbwith(istep,idf)=0.
      enddo
      dsum=0.
      dcnt=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           print*,'idf=',idf,'nsizebig=',nsizebig,'nsizelit=',nsizelit
           if(idf.le.nsizebig) then 
               xdlam=(wavbigtilda(idf+1)-wavbigtilda(idf)) !width of the interval
               xlam=(wavebig(idf+1)-wavebig(idf)) !width of the interval
               xlambda0=wavbigtilda(idf)
            if(wavbigtilda(idf+1).gt.0) then
!            print*,'idf=',idf,'xdlam=',xdlam,'ntot=',ntot
!! first end point 
             nnsteps(1,idf)=1
!! subsequent steps
            nsum=0
         do j1=1,ntot !high resolution grid per bin
             jj1=indexx(j1)
             jj=indexx(j1)-nbegin+1 !! high resoultion original index
            
!            print*,'j1=',j1
            if(j1==1)xlambda=xlambda0
            if(j1>1) then
              xlambda=xlambda+dlamtilda(idf,jj) !non-equidistant high res lambdatilda grid
            endif
            if(j1==ntot) xlambda=wavbigtilda(idf+1)

            do isub=2,2

             if( xlambda.ge.(xdlam*weights(isub-1)+xlambda0).and. &
              xlambda.le.(xdlam*weights(isub)+xlambda0) )then

!             if(xlambda.gt.0.) write(*,"(6(1x,e15.6))") xlambda,dlamtilda(idf,jj),phi_hr(idf,jj),(buffer(jj1))
!             write(*,"(3(1x,i4),2(1x,e15.6))")idf,isub,j1,xdlam*weights(isub-1)+xlambda0,xdlam*weights(isub)+xlambda0

              nsum=nsum+1
              nnsteps(isub,idf)=nnsteps(isub,idf)+1
!             print*,'isub',isub,'nsum=',nsum,'j1',j1
             endif
            enddo
            do isub=3,nnsub+1

              if( xlambda.gt.(xdlam*weights(isub-1)+xlambda0).and. &
              xlambda.le.(xdlam*weights(isub)+xlambda0) )then

!             if(xlambda.gt.0.)write(*,"(6(1x,e15.6))")xlambda,dlamtilda(idf,jj),phi_hr(idf,jj),(buffer(jj1))
!             write(*,"(3(1x,i4),2(1x,e15.6))")idf,isub,j1,xdlam*weights(isub-1)+xlambda0,xdlam*weights(isub)+xlambda0

               nsum=nsum+1
               nnsteps(isub,idf)=nnsteps(isub,idf)+1
!             print*,'isub',isub,'nsum=',nsum,'j1',j1

!               if(isub .eq. 3) then
!               dsum=dsum+dlamtilda(idf,jj)
!               dcnt=dcnt+1
!               print*,'sb3',j1,dsum,xdlam*(weights(isub)-weights(isub-1))
!               endif

              endif
            enddo

         enddo ! j=1,ntot loop 
!         print*,'big-bin: nsum=',nsum,'ntot',ntot
           ii=1
           !print_if_needed
          !print*,' nnsteps-nsteps',ii,nnsteps(ii,idf),nsteps(ii,idf)
           ii=2
           !print_if_needed
          !print*,' nnsteps-nsteps',ii,nnsteps(ii,idf),nsteps(ii,idf)
           do ii=3,nnsub+1!do not add the first boundary value 1 to end of 1st bin//start from 2nd bin//
!           print*,'weights',ii,weights(ii)
           if(weights(ii).ne.0.) then
           nnsteps(ii,idf)=nnsteps(ii-1,idf)+nnsteps(ii,idf)
           endif
           if(weights(ii).eq.0.) then
           nnsteps(ii,idf)=1
           endif
           
           !print_if_needed
           !print*,'nnsteps-nsteps',ii,nnsteps(ii,idf),nsteps(ii,idf)
           enddo
         endif !wavbigtilda gt 0
        endif ! idf le nsizebig 

        if(idf.gt.nsizebig) then
               xdlam=(wavlittilda(idf-nsizebig+1)-wavlittilda(idf-nsizebig)) !width of the interval
               xlam=(wavelit(idf-nsizebig+1)-wavelit(idf-nsizebig)) !width of the interval
               xlambda0=wavlittilda(idf-nsizebig)
            if(wavlittilda(idf-nsizebig+1).gt.0) then
!! first end point 
             nnsteps(1,idf)=1
!! subsequent steps
            nsum=0
          do j1=1,ntot !high resolution grid per bin
           jj=indexx(j1)-nbegin+1 !! high resoultion original index

!            print*,'j1=',j1
            if(j1==1)xlambda=xlambda0
            if(j1>1) xlambda=xlambda+dlamtilda(idf,jj) !non-equidistant high res lambdatilda grid
            if(j1==ntot) xlambda=wavlittilda(idf-nsizebig+1)

            do isub=2,2
             if( xlambda.ge.(xdlam*weights(isub-1)+xlambda0).and. &
              xlambda.le.(xdlam*weights(isub)+xlambda0) )then
                nsum=nsum+1
                nnsteps(isub,idf)=nnsteps(isub,idf)+1
!             print*,'isub=',isub,'nsum=',nsum,'j1',j1,'w=',weights(isub-1),weights(isub),'xdlam=',xdlam
              endif
            enddo
            do isub=3,nnsub+1
              if( xlambda.gt.(xdlam*weights(isub-1)+xlambda0).and. &
               xlambda.le.(xdlam*weights(isub)+xlambda0) )then
                nsum=nsum+1
                nnsteps(isub,idf)=nnsteps(isub,idf)+1
!             print*,'isub=',isub,'nsum=',nsum,'j1',j1,'w=',weights(isub-1),weights(isub),'xdlam=',xdlam
              endif
           enddo
         enddo ! j=1,ntot loop 

!         print*,'lit-bin: nsum=',nsum,'ntot',ntot
           ii=1
!           print*,'lit-bin: nnsteps-nsteps',ii,nnsteps(ii,idf),nsteps(ii,idf)
           ii=2
!           print*,'lit-bin: nnsteps-nsteps',ii,nnsteps(ii,idf),nsteps(ii,idf)
           do ii=3,nnsub+1 !do not add the first boundary value 1 to end of 1st bin//start from 2nd bin//
           if(weights(ii).ne.0.) then
           nnsteps(ii,idf)=nnsteps(ii-1,idf) &
                             +nnsteps(ii,idf)
           endif
           if(weights(ii).eq.0.) then
           nnsteps(ii,idf)=1
           endif
!           print*,'lit-bin: nnsteps-nsteps',ii,nnsteps(ii,idf),nsteps(ii,idf)
           enddo
         endif !wavlittilda gt 0
        endif ! idf gt nsizebig 
     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       sbsum=0.
       if(ifilter1==0) then
        if(idf.le.nsizebig) then
        xlam=(wavebig(idf+1)-wavebig(idf)) !width of the interval
        endif !idf le nsizebig
        if(idf.gt.nsizebig) then
        xlam=(wavelit(idf-nsizebig+1)-wavelit(idf-nsizebig)) !width of the interval
        endif !idf gt nsizebig
       endif
       

!       write(3,*) idf,ntot
      do istep=1, nnsub
        isum=0
        isum1=0
        isum2=0
        isum3=0

        sum1=0.
        sumo=0.
        err0=0.
        fsum=0.
        sump=0.

      if(ifilter1==0) then


        do i=nsteps(istep,idf),nsteps(istep+1,idf)-1 !!starts with 1 for each bin and ends with ntot

!          isum=isum+ibuffsort(i) !original

           j1=indexx(i)!! high resoultion original index
           jj=indexx(i)-nbegin+1 !! high resoultion original index
           if(ameanODF .eq. 0) then
           !!!GM
           fsum=fsum+log10(buffer(j1))
           else if(ameanODF .eq. 1) then
           !!!AM
           fsum=fsum+(buffer(j1))
           else if(ameanODF .eq. 2) then
            !HM
            if(buffer(j1) .ne. 0.) then 
            fsum=fsum+(1.0/buffer(j1))
            endif
           endif

           sump=sump+phi_hr(idf,jj)
           sumpb=sumpb+phi_hr(idf,jj)
!        print*,'idf=',idf,'i=',i,'jj=',jj,'phi_hr=',phi_hr(idf,jj)

!!!!!!!!!!print phi_hr to a file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if(idf.le.nsizebig) then
         xdlam=(wavbigtilda(idf+1)-wavbigtilda(idf)) !width of the interval
         if(i==1) lambda=wavebig(idf)
         if(i>1) lambda=lambda+dlambda(idf,i) 
         if(i==ntot) lambda=wavebig(idf+1)
         endif
         if(idf.gt.nsizebig) then
         xdlam=(wavlittilda(idf+1)-wavlittilda(idf)) !width of the interval
         if(i==1) lambda=wavelit(idf-nsizebig)
         if(i>1) lambda=lambda+dlambda(idf,i) 
         if(i==ntot) lambda=wavelit(idf+1-nsizebig)
         endif
!       write(3,"(6(1x,e15.6))") &
!                      & lambda,& !!new grid
!                       phi_hr(idf,jj) 
!       if(idf==2)write(*,"(1(1x,i3),6(1x,e15.6))") &
!                      & jj,lambda,& !!new grid
!                       phi_hr(idf,jj) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                       

        enddo!i loop
!          insteps(istep)=(isum+int(nsteps(istep+1,idf)-nsteps(istep,idf),8)/2)/ &
!                   & int(nsteps(istep+1,idf)-nsteps(istep,idf),8)

            fn=dble(nsteps(istep+1,idf)-nsteps(istep,idf))
            if(ameanODF .eq. 0) then
             !!!GM
            fsum=fsum/fn
            !isum=nint(fsum*1000.)
            isum=int(fsum*1000.,8)
            else if(ameanODF .eq. 1) then
             !!!AM
            fsum=fsum/fn
            !isum=nint(log10(fsum)*1000.)
            isum=int(log10(fsum)*1000.,8)
            elseif(ameanODF .eq. 2) then 
            fsum=fn/fsum
            isum=int(log10(fsum)*1000.,8)
            endif

            !insteps(istep)=isum
            insteps(istep)=int(isum,2)

            sump=sump/fn !filter averaged over sub-bin
!           sbwith(istep,idf) = sbwith(istep,idf)*xlam
!           sbwith(istep,idf) = sbwith(istep,idf)*sump*xlam
!           sbwith(istep,idf) = xdlam*(weights(istep+1)-weights(istep))
           
!   print*,'idf=',idf,'istep=',istep,xdlam,sbwith(istep,idf),weights(istep),weights(istep+1),(weights(istep+1)-weights(istep))

!!!!!!!write to files!!!!!!!!!
!       do i=nsteps(istep,idf),nsteps(istep+1,idf)-1 !!starts with 1 for each bin and ends with ntot
!
!          j1=indexx(i)!! high resoultion original index
!          jo=nbegin-i+1     !!original index
!          jj=indexx(i)-nbegin+1 !!re-numbered from 1 to ntot
!
!         if(idf.le.nsizebig) then
!         if(i==1) lambda=wavebig(idf)
!         if(i>1) lambda=lambda+dlambda(idf,i) 
!         if(i==ntot) lambda=wavebig(idf+1)
!         endif
!         if(idf.gt.nsizebig) then
!         if(i==1) lambda=wavelit(idf-nsizebig)
!         if(i>1) lambda=lambda+dlambda(idf,i) 
!         if(i==ntot) lambda=wavelit(idf+1-nsizebig)
!         endif

!       if(idf==657) then
!         if(idf.le.nsizebig) then
!         write(1,"(6(1x,e15.6))") &
!                      & lambda,& !!new grid
!                      & buffer(jo),&  !!un-sorted opacity
!                      & buffer(j1),&  !!sorted opacity
!                      & 10.**(dble(insteps(istep))/1000.)
!         else
!         write(3,"(6(1x,e15.6))") &
!                      & lambda,& !!new grid
!                      & buffer(jo),&  !!un-sorted opacity
!                      & buffer(j1),&  !!sorted opacity
!                      & 10.**(dble(insteps(istep))/1000.)
!        endif
!     endif

!          if(idf==657) then
!            & print*,'lambda',i,jo,j1,lambda,dlambda(idf,i),&
!            & buffer(jo),10.**(dble(ibuffer(jo))/1000.),& !!un-sorted
!            & buffer(j1),10.**(dble(ibuffsort(i))/1000.) !!sorted
!             & (sum1)/(1.0*ntot)
!           endif
!        enddo!i loop
           !print_if_needed
      !    write(*,"(A8,I7,3(1x,A8,I7,1x,A8,e15.6))") "istep=",istep,&
      !                & "isum=",isum,"odf=",exp(dble(isum)/1000.)



      else if(ifilter1==1) then
!test to mimic odf
!      if(idf.eq.2.or.idf.eq.5 .or.idf.eq.6.or.idf.eq.7.or.idf.eq.8) then
!      if(idf.eq.2.or.idf.eq.6.or.idf.eq.7.or.idf.eq.8) then
!      nnsteps(istep,idf)=nsteps(istep,idf)
!      nnsteps(istep+1,idf)=nsteps(istep+1,idf)
!      print*,'IDF=',idf,'changed'
!      endif
!!!!! 
         idiff=nnsteps(istep+1,idf)-nnsteps(istep,idf)

        if(idiff.gt.0) then
!         print*,' '
!         print*,'nnsteps(istep,idf)=',nnsteps(istep,idf)
!         print*,'nnsteps(istep+1,idf)=',nnsteps(istep+1,idf)
!         print*,'idiff!=0',' ','idf=',idf,'istep',istep

        xden=0.
        xsum=0.
        xcnt=0
         do i=nnsteps(istep,idf),nnsteps(istep+1,idf)-1
!         do i=nnsteps(istep,idf)+1,nnsteps(istep+1,idf)
         !!!!!! connect high resolution grid to nnsteps !!!
          jj=indexx(i)-nbegin+1 !! high resoultion original index
          
         xden=xden+dlamtilda(idf,jj)
         
!test to mimic odf
!         xden=xden+dlambda(idf,i)

         xsum=xsum+dlambda(idf,i)
        xcnt=xcnt+1
        enddo!i loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!             sbwith(istep,idf) = xden
            sbwith(istep,idf) = xdlam*(weights(istep+1)-weights(istep))
!If I divide by xdlam, we see last weight equlas one, but we need to multiply while adding the I over bins..
!          sbwith(istep,idf) = xden/xlam
!
!            xden=xdlam*sbwith(istep,idf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!          sbwith(istep,idf) = xden/xsum !!this is wrong!! bec xsum is just delta_lambda summed over only that sub-bin numbers!!
          sbsum=sbsum+sbwith(istep,idf)

!         print*,'sbwith',istep,idf,sbwith(istep,idf),xdlam*(weights(istep+1)-weights(istep))

          wsum=0.
         do i=nnsteps(istep,idf),nnsteps(istep+1,idf)-1
          j1=indexx(i)!! high resoultion original index
          jo=nbegin-i+1     !!original index
          jj=indexx(i)-nbegin+1 !!re-numbered from 1 to ntot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WGM
           !wsum=wsum+(dlamtilda(idf,jj)/xden)
           if(ameanODF .eq. 0) then
            if(idf.le.nsizebig)xnum=(log10(buffer(j1)))*(dlamtilda(idf,jj))
            if(idf.gt.nsizebig)xnum=(log10(buffer(j1)))*(dlamtilda(idf,jj))
!GM
!           if(idf.le.nsizebig)xnum=(log10(buffer(j1)))
!           if(idf.gt.nsizebig)xnum=(log10(buffer(j1)))
!            fn=dble(nsteps(istep+1,idf)-nsteps(istep,idf))
            else if(ameanODF .eq. 1) then 

!AM
           if(idf.le.nsizebig)xnum=((buffer(j1)))*(dlamtilda(idf,jj))
           if(idf.gt.nsizebig)xnum=((buffer(j1)))*(dlamtilda(idf,jj))
          endif

!          if(idf.le.nsizebig)xnum=(log10(buffer(j1)))*(dlamtilda(idf,jj))*sbwith(istep,idf)
!          if(idf.gt.nsizebig)xnum=(log10(buffer(j1)))*(dlamtilda(idf,jj))*sbwith(istep,idf)
!            fn=dble(nsteps(istep+1,idf)-nsteps(istep,idf))
!          if(idf.le.nsizebig)xnum=0.5*(log10(buffer(j1)))*((1./fn)+dlamtilda(idf,jj)/xden)
!          if(idf.gt.nsizebig)xnum=0.5*(log10(buffer(j1)))*((1./fn)+dlamtilda(idf,jj)/xden)
!          if(idf.le.nsizebig)xnum=(log10(buffer(j1)))*((1./fn))
!          if(idf.gt.nsizebig)xnum=(log10(buffer(j1)))*((1./fn))
!          xden=1.


!test to mimic odf
!          if(idf.le.nsizebig)xnum=(log10(buffer(j1)))*dlambda(idf,i)
!          if(idf.gt.nsizebig)xnum=(log10(buffer(j1)))*dlambda(idf,i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!test to mimic odf
!!!for unit filter, the following gives exact results!!
!          if(idf.le.nsizebig)xnum=(log10(buffer(j1)))
!          if(idf.gt.nsizebig)xnum=(log10(buffer(j1)))
!          xden=1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         if(abs(xden).gt.0.)then
!!GM
!          xden=fn
!!WGM
          wavg=xnum/xden
         else
          wavg=0.
         endif
         if(abs(wavg).gt.0.) then
         !jsum=nint(wavg*1000.) 
         !jsum2=nint(wavg*10000.) 
         jsum=int(wavg*1000.,8) 
         jsum2=int(wavg*10000.,8) 

         else
         jsum=0 
         jsum2=0 
         endif

         isum=isum+jsum
         !isum1=isum1+ibuffsort(i)
         isum1=isum1+int(ibuffsort(i),8)
         isum2=isum2+jsum2
         
!!!ADDING THE FLOAT VALUES FOR WEIGHTD GM AND THEN TAKING THE INTEGER PART IS VERY ACCURATE
!!!THIS HAS BEEN TESTED FOR UNIT FILTER AND GIVES VERY LITTLE REL ERR IN FLUX.!!!!!
!!!THUS THE FOLLOWING IS A VERY VERY IMP STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          sum1=sum1+wavg
          sumo=sum1
          sum1=sumo+wavg
!!!Kahan summation
!          err0=err0+(wavg-(sum1-sumo))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
!         write(*,"(A4,I4,3(1x,A4,e15.6))")"i=",i,"buf=",buffer(j1),"lgb=",log10(buffer(j1)),"dlam=",dlamtilda(idf,jj) 
!
!         write(*,"(A4,I4,3(1x,A4,e15.6))")'i=',i,'xnum=',xnum,'xden=',xden,'wavg=',wavg 
!         write(*,"(A4,I4,4(1x,A6,I6))")'i=',i,'jsum=',jsum,'isum=',isum,"jsum2=",jsum2,"isum2=",isum2

!           if(idf .eq. 5 .and. istep.eq.6) then 
!            print*,'j1=',j1,'buffer=',buffer(j1),'dlamt=',dlamtilda(idf,jj)
!            print*,'i=',i,'xnum=',xnum,'xden=',xden,'sum1=',sum1,((dble(isum1)/(fn*1000.)))
!            print*,10.**sum1,10.**((dble(isum1)/(fn*1000.)))
!           endif

        enddo !!i loop
!         print*,'wsum=',wsum
!!!THUS THE FOLLOWING IS A VERY VERY IMP STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         sum1=sum1+err0
!         isum3=nint(sum1*10000.)
!         isum3=nint(log10(sum1)*1000.)

          !!WGM
         if(ameanODF .eq. 0) then
         !isum3=nint(sum1*1000.)
         isum3=int(sum1*1000.,8)
          !!AM
         else if(ameanODF .eq. 1) then
         !isum3=nint(log10(sum1)*1000.)
         isum3=int(log10(sum1)*1000.,8)
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         write(*,"(A4,I4,3(1x,A4,e15.6))")'i=',i,'odf=',10.**(sum1)

!!The following step is necessary only when using geometric mean; NOT when using weighted geometric mean.
!!This is because, in GM we divide by N, and in weighted GM we divide by float sum of weights.
!!!for unit filter, the following gives exact results!!!!!!!!!!!!!!!!!!!
!!!for unit filter, there is slight difference is nsteps and nnsteps, 
!!!!but it is negligible for flux rel.err.!!!!
!         isum=(isum+int(nsteps(istep+1,idf)-nsteps(istep,idf))/2,8)/ &
!                   & int(nsteps(istep+1,idf)-nsteps(istep,idf),8)
!
!         isum2=(isum2+int(nsteps(istep+1,idf)-nsteps(istep,idf))/2,8)/ &
!                   & int(nsteps(istep+1,idf)-nsteps(istep,idf),8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!filtered odfs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         isum=(isum+int(nnsteps(istep+1,idf)-nnsteps(istep,idf))/2,8)/ &
!                   & int(nnsteps(istep+1,idf)-nnsteps(istep,idf),8)

!         isum2=(isum2+int(nnsteps(istep+1,idf)-nnsteps(istep,idf))/2,8)/ &
!                   & int(nnsteps(istep+1,idf)-nnsteps(istep,idf),8)

!test to mimic odf
!         isum3=(isum3+int(nnsteps(istep+1,idf)-nnsteps(istep,idf),8)/2)/ &
!                   & int(nnsteps(istep+1,idf)-nnsteps(istep,idf),8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!Normal odf to comapre!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !isum1=(isum1+(nsteps(istep+1,idf)-nsteps(istep,idf))/2)/ &
         !          & (nsteps(istep+1,idf)-nsteps(istep,idf))
         isum1=(isum1+int(nsteps(istep+1,idf)-nsteps(istep,idf),8)/2)/ &
                   & int(nsteps(istep+1,idf)-nsteps(istep,idf),8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        insteps(istep)=int(isum,2)
!        insteps(istep)=int(isum2,2)
!
!test to mimic odf!int sum instead of float sum as in orig odf
!         insteps(istep)=int(isum1,2)
!!!THUS THE FOLLOWING IS A VERY VERY IMP STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !insteps(istep)=isum3
         insteps(istep)=int(isum3,2)
!integer sum
!         insteps(istep)=isum
!         insteps(istep)=int(isum,2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!write to files!!!
!         do i=nnsteps(istep,idf),nnsteps(istep+1,idf)-1
!          j1=indexx(i)!! high resoultion original index
!          jo=nbegin-i+1     !!original index
!          jj=indexx(i)-nbegin+1 !!re-numbered from 1 to ntot
!         if(idf.le.nsizebig) then
!         if(i==1) lambdatilda=wavbigtilda(idf)
!         if(i>1) lambdatilda=lambdatilda+dlamtilda(idf,jj) 
!         if(i==ntot) lambdatilda=wavbigtilda(idf+1)

!         if(i==1) lambdatildaus=wavbigtilda(idf)
!         if(i>1) lambdatildaus=lambdatildaus+dlamtilda(idf,i) 
!         if(i==ntot) lambdatildaus=wavbigtilda(idf+1)
!         endif
!         if(idf.gt.nsizebig) then
!         if(i==1) lambdatilda=wavlittilda(idf-nsizebig)
!         if(i>1) lambdatilda=lambdatilda+dlamtilda(idf,jj) 
!         if(i==ntot) lambdatilda=wavlittilda(idf+1-nsizebig)
!
!         if(i==1) lambdatildaus=wavlittilda(idf-nsizebig)
!         if(i>1) lambdatildaus=lambdatildaus+dlamtilda(idf,i) 
!         if(i==ntot) lambdatildaus=wavlittilda(idf+1-nsizebig)
!         endif

!      if(idf==657) then
!         if(idf.le.nsizebig) then
!         write(2,"(6(1x,e15.6))") &
!                      & lambdatildaus,& !!new grid, unsorted
!                      & lambdatilda,& !!new grid
!                      & buffer(jo),&  !!un-sorted opacity
!                      & buffer(j1),&  !!sorted opacity
!                      & 10.**(sum1) 
!         else
!         write(*,"(6(1x,e15.6))") &
!                      & lambdatildaus,& !!new grid, unsorted
!                      & lambdatilda,& !!new grid
!                      & buffer(jo),&  !!un-sorted opacity
!                      & buffer(j1),&  !!sorted opacity
!                      & 10.**(sum1) 
!         
!         write(4,"(6(1x,e15.6))") &
!                      & lambdatildaus,& !!new grid, unsorted
!                      & lambdatilda,& !!new grid
!                      & buffer(jo),&  !!un-sorted opacity
!                      & buffer(j1),&  !!sorted opacity
!                      & 10.**(sum1) 
!         endif
!      endif
!        enddo !i loop

           !print_if_needed
          !write(*,"(A8,I7,3(1x,A8,I7,1x,A8,e15.6))") "istep=",istep,&
          !            & "isum1=",isum1,"odf=",exp(dble(isum1)/1000.),&
          !            & "isum=",isum,"fodf1=",exp(dble(isum)/1000.),&
          !            & "isum2=",isum2,"fodf2=",exp(dble(isum2)/10000.),&
          !            & "isum3=",isum3,"fodf3=",exp(dble(isum3)/1000.)

        endif !idiff gt 0
       endif !ifilter1==0 or 1
      end do !istep loop
!       print*,'before: sbsum=',sbsum
!!!!!!!re-normalize sub-bin-weights to make it one!!
!!
!      do istep=1, nnsub
!       sbwith(istep,idf)=sbwith(istep-1,idf)+sbwith(istep,idf) !This is wrong!!
!       sbwith(istep,idf)=sbwith(istep,idf)/sbsum
!       print*,'sub-bin',sbwith(istep,idf),(weights(istep+1)-weights(istep))
!      enddo
!       sbsum=0.
!      do istep=1, nnsub
!       sbsum=sbsum+sbwith(istep,idf)
!      enddo
!       print*,'after: sbsum=',sbsum
      

!      idfout(1)=int(vt*100)*1000000+idf*10000+iteff*100+j
!       
!      if(iv.eq.1)write(15)idfout
!      if(iv.eq.2)write(16)idfout
!      if(iv.eq.3)write(17)idfout
!      if(iv.eq.4)write(18)idfout
!      if(iv.eq.5)write(20)idfout

!!!!!TEST!!!
     do i = 1, 20
!     do i = 1,1 
      if(iv.eq.i) then
       do istep = 1, nnsub
        iodfstep(istep,idf,j,iteff,i)=insteps(istep)
       end do
       end if
     end do

   fnb=dble(ntot)
   sumpb=sumpb/fnb
!   print*,'idf=',idf,'sumpb=',sumpb
   !!xdlam is the sub-bin averaged filter weight to each bin!!
   !!xdlam is equilvalent to bin-averaged filter-function*delta_lambda that we multiply to 
   !!bin averaged intensity!!!
!   if(idf.le.nsizebig .and. iteff.eq.1 .and. j.eq.1) write(112,*) xdlam
  end do !idf loop
!   STOP
  return
end subroutine
