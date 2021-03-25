subroutine dfintervals_filter
!
!     The intervals big and little  remain  hardcoded for now 
!     The amount of sub-bins can/should be changed here
!------------------------------------------------------------------------------
      use types
      use dfcomm
      implicit none
      
      include 'common.ifopbl'

      real(kind=dp) :: nzero
      integer :: i, j, litbig, check 
      integer :: ip, jj 

!-----------------------------------------------------------------------------

      call def_binsize()
!   if not Kuruzc grid then wavelit, wavebig and nsizes are defined in def_binsize
 
! high resolution gird is defined as:    l_n = l_0 *(1 + 1/res)**(n-1)
! lg(l_n) / lg(1 + 1/res) = lg(l_0)/lg(1+1/res) + n -1
!-------------------------------------------------------------------------

! first the bin's upper boundary is found using the equation above, where l_0 = wlbeg
!---------------------------------------------------------------------------
      do  i=1,nsizebig 
         nbigend(i)=log(wavebig(i+1))/log(1.d0+1.d0/500000.d0)
      end do
      
      nzero=log(wlbeg)/log(1.d0+1.d0/500000.d0)
      nzero=nzero+1 
      wbegin=(1.d0+1.d0/500000.d0)**nzero

      do  i=1,nsizebig
       nbigend(i)=nbigend(i)-nzero+1
      end do
!-------------- then the end of each bin is used to determine the beginning 
! of the next bin ----------------------------------------------------------      
      do  i=2,nsizebig
       nbigbeg(i)=nbigend(i-1)+1
      end do
!--- the first bin's beginning is determined
      
      nbigbeg(1)= int(log(wavebig(1))/log(1.0d0 + 1.0d0/500000.0d0)-nzero) + 1 
      if (nbigbeg(1) .lt. 1 ) nbigbeg(1) = 1

      wend=(1.d0+1.d0/500000.d0)**(nbigend(nsizebig)+nzero-1)

!--- amount of grid points on the high-resolution grid in each bin is calculated:

      do  i=1,nsizebig
        nn(i)=nbigend(i)-nbigbeg(i)
      end do

!  same procedure as above for Kurucz little  gird
      
      do  i=1,nsizelit
        nlitend(i)=log(wavelit(i+1))/log(1.d0+1.d0/500000.d0)
      end do

      nzero=log(wlbeg)/log(1.d0+1.d0/500000.d0)
      nzero=nzero+1
      wbegin=(1.d0+1.d0/500000.d0)**nzero

      do i=1,nsizelit
        nlitend(i)=nlitend(i)-nzero+1 
      end do

      if (ifkbin) then 
        do i=2,nsizelit
          nlitbeg(i)=nlitend(i-1)+1 
        end do
      end if 
      
      nlitbeg(1)= int(log(wavelit(1))/log(1.0d0 + 1.0d0/500000.0d0)-nzero)+1
      if (nlitbeg(1) .lt. 1) nlitbeg(1) = 1

      wend=(1.d0+1.d0/500000.d0)**(nlitend(nsizelit)+nzero-1)

      do  i=1,nsizelit
         nn(i+nsizebig)=nlitend(i)-nlitbeg(i)
      end do



!-------------------------------------------------------------------
!      Here the sub-bin sizes are defined !!!
!------------------------------------------------------------------

     if (ifkbin ) then 


      do  i=1, nsizebig+ nsizelit !1540 
!     1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/20.1/30,1/60
!     1/10,2/10,3/10,4/10,5/10,6/10.7/10.8/10,9/10,19/20,59/60,60/60
!
         nsteps(1,i)=1
         nsteps(2,i)=(nn(i)*6)/60+1
         nsteps(3,i)=(nn(i)*12)/60+1
         nsteps(4,i)=(nn(i)*18)/60+1
         nsteps(5,i)=(nn(i)*24)/60+1
         nsteps(6,i)=(nn(i)*30)/60+1
         nsteps(7,i)=(nn(i)*36)/60+1
         nsteps(8,i)=(nn(i)*42)/60+1
         nsteps(9,i)=(nn(i)*48)/60+1
         nsteps(10,i)=(nn(i)*54)/60+1
         nsteps(11,i)=(nn(i)*57)/60+1
         nsteps(12,i)=(nn(i)*59)/60+1
         nsteps(13,i)=nn(i)+1
      end do
     
     else 

       nsteps = 0 
! --- for user defined sub-bin sizes the input file subbin_info.dat is read and the nsteps 
! --- are calculated.
        

      open (unit = 78, status = 'old', file='./INPUT/subbin_info.dat', form = 'formatted')
!      open (unit = 78, status = 'old', file='./INPUT/sub_bin_sizes.dat', form = 'formatted')
      read(78,*) nsubbin, check 
!      write(6,*) 'in intervals, read subbin_info' 
!      write(6,*) 'nsubbin, check', nsubbin, check      
       
      do i = 1, nsizebig+ nsizelit 
         read(78,*) (binbound(i,j), j=1, nsubbin) 
         ! here the upper boundaries are read, for example for a bin with 3 subbins:
         ! 0.3, 0.6, 1.0 ---------------------------------------------------------- 
 
          nsteps(1,i)=1
          sbwith(1,i) = binbound(i,1)
        do j = 2, nsubbin 
          sbwith(j,i) = binbound(i,j)-binbound(i,j-1)       
          nsteps(j,i) = (nn(i)*int(1000*binbound(i,j-1)))/1000 +1

!         print*,'binbound=',binbound(i,j-1)
!         print*,'nsteps=',nsteps(j,i),nn(i) 
        end do 
         nsteps(nsubbin+1,i) = (nn(i)*int(1000*binbound(i,nsubbin)))/1000 +1

!         print*,'binbound=',binbound(i,nsubbin)
!         print*,'nsteps=',nsteps(nsubbin+1,i),nn(i) 
      end do 

      close(unit =78) 


     end if 
!--- finally the global arrays nbeg, nend which are equivlent to nbig and nlit end & begin are 
!--- assigned. These arrays are used in df.calc together with nsteps to calculated the opacity averages

      litbig = nsizebig+nsizelit

      do i = 1, nsizebig
       nbeg(i) = nbigbeg(i)
       nend(i) = nbigend(i)
!     print*,'i=',i,'nbeg=',nbeg(i),'nend=',nend(i)
      end do 


      do i = nsizebig+1, litbig
       nbeg(i) = nlitbeg(i-nsizebig)
       nend(i) = nlitend(i-nsizebig)
!     print*,'i=',i,'nbeg=',nbeg(i),'nend=',nend(i)
      end do 
!clsa
!!clsa
!!CREATE dlamtilda HERE!!
         if(ifilter==1) then
         do i=1,litbig
          do j=1,lenbuff
          dlamtilda(i,j)=0.
          enddo
         enddo

        do i=1,nsizebig+1
        wavbigtilda(i)=0.
!        print*,'i=',i,'wavbig=',wavebig(i)
        enddo

        jj=0
        do i=1,nsizebig
        nntot=nn(i)+1
        !write(3,*)i,nntot
        dlam=(wavebig(i+1)-wavebig(i))/(1.0*(nntot-1))

        lambdatilda=wavbigtilda(i)
        totdlam=0. 
          ic=0
!        print*,'i=',i,'nnot=',nntot
!        print*,'wavbig:',wavebig(i),wavebig(i+1)
!        print*,'phi:',xlambda_phi(1),xlambda_phi(Nlambda_phi)
        do j=1,nntot
        if(j.eq.1)lambda=wavebig(i)
        if(j>1)lambda=lambda+dlam
        if(j.eq.nntot)lambda=wavebig(i+1)
        dlambda(i,j)=dlam !! high resolution lambda grid

!        jj=jj+1
!        dlambda(i,jj)=dlam !! high resolution lambda grid

        phi_lp=0.
        phi_avg=0.
!        print*,'Nlambda_phi',Nlambda_phi

        do ip=1,Nlambda_phi-1 !bec of ip+1
!            print*,'xlambda_phi=',xlambda_phi(ip)
          if(lambda>=xlambda_phi(ip) .and.  &
           lambda<=xlambda_phi(ip+1) ) then
!interpolated
         phi_lp=((filter(ip+1)-filter(ip))/ &
              (xlambda_phi(ip+1)-xlambda_phi(ip)) )* &
              (lambda-xlambda_phi(ip))+filter(ip)
!arithmatic average
         phi_avg=0.5*(filter(ip)+filter(ip+1))

!         print*,'xlambda_phi=',xlambda_phi(ip),xlambda_phi(ip+1)
!         print*,'ip',ip,'filter',filter(ip+1),filter(ip)
!         print*,i,j,ip,'phi_lp',phi_lp
          ic=ic+1
          endif
        enddo !~ip
          dlamtilda(i,j)=phi_lp*dlam 
!         dlamtilda(i,j)=phi_avg*dlam 
         phi_hr(i,j)=phi_lp
!       write(3,"(6(1x,e15.6))") &
!                      & lambda,& !!new grid
!                       phi_hr(i,j)
!          if(i==2) print*,'i=',i,'j=',j,'phi=',phi_hr(i,j)
!!!!!!!!!!!!!!!end points must be defined based on filter profile!!!!
!!!!!!!!!!the following are wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         if(j==1) wavbigtilda(i)=phi_lp*wavebig(i) 
!         if(j==nntot) wavbigtilda(i)=phi_lp*wavebig(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        lambdatilda=lambdatilda+dlamtilda(i,j)
        totdlam=totdlam+dlamtilda(i,j)
!        if(dlamtilda(i,j)>0.) then
!         print*,'i=',i,'j=',j,phi_lp,dlamtilda(i,j),dlam
!        endif

!for filter=1 test!
!        if(abs(dlamtilda(i,j)-dlam)>0.) then
         !print*,'dlam=',i,j,lambda,phi_lp,dlamtilda(i,j),dlam,totdlam
!        endif


        enddo !j loop
       if(ic>=nsubbin) then
       wavbigtilda(i+1)=wavbigtilda(i)+totdlam
      !print*,'ic=',ic,'wavbigt',i+1,wavbigtilda(i),wavbigtilda(i+1)
       endif
!       else
!       print*,'wavbigt',i+1,wavbigtilda(i),wavbigtilda(i+1),'ic=',ic
!       endif
!       if(totdlam.ne.0) then
!       print*,'i=',i,'ic',ic,'nntot',nntot
!       endif
!       else
!       print*,'i=',i,'ic',ic,'nntot',nntot,'totdlam=',totdlam
!       endif

       enddo !i loop
        
       do i=1,nsizelit+1
       wavlittilda(i)=0.
       enddo

!       print*,' '

       do i=1,nsizelit
       nntot=nn(i+nsizebig)+1
!       print*,'nntot',nntot
!       print*,' '
       dlam=(wavelit(i+1)-wavelit(i))/(1.0*(nntot-1))

       lambdatilda=wavlittilda(i)
       totdlam=0.
        ic=0
       do j=1,nntot
        if(j.eq.1)lambda=wavelit(i) !initialize
        if(j>1)lambda=lambda+dlam
        if(j.eq.nntot)lambda=wavelit(i+1)
        dlambda(i+nsizebig,j)=dlam !! high resolution lambda grid

!!        jj=jj+1
!!        dlambda(i+nsizebig,jj)=dlam !! high resolution lambda grid

!         if(i.eq.334)print*,j,'lambda=',lambda
         phi_lp=0.
         phi_avg=0.
         do ip=1,Nlambda_phi-1 !bec of ip+1
          if(lambda>=xlambda_phi(ip) .and. &
           lambda<=xlambda_phi(ip+1) ) then
!interpolated
         phi_lp=((filter(ip+1)-filter(ip))/ &
                 (xlambda_phi(ip+1)-xlambda_phi(ip)) )* &
                 (lambda-xlambda_phi(ip))+filter(ip)
!arithmatic average
         phi_avg=0.5*(filter(ip)+filter(ip+1))
!         print*,'xlambda_phi=',xlambda_phi(ip),xlambda_phi(ip+1)
!         print*,'ip',ip,'filter',filter(ip+1),filter(ip)
!         print*,'phi_lp',phi_lp
           ic=ic+1
          endif
         enddo !~ip
         dlamtilda(i+nsizebig,j)=phi_lp*dlam
!         dlamtilda(i+nsizebig,j)=phi_avg*dlam
!         phi_hr(i,j)=phi_lp

!!!!!!!!!!!!!!!end points must be defined based on filter profile!!!
!!!!!!!!!!the following are wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         if(j==1) wavlittilda(i)=phi_lp*wavelit(i)
!         if(j==nntot) wavlittilda(i)=phi_lp*wavelit(i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        lambdatilda=lambdatilda+dlamtilda(i+nsizebig,j)
        totdlam=totdlam+dlamtilda(i+nsizebig,j)

!        if(i+nsizebig>=657 .and. i+nsizebig<=662) then
!        print*,'j=',j,'dlam',dlam,'dlamtilda',dlamtilda(i+nsizebig,j),'phi_lp',phi_lp,'wavlit',wavlittilda(i)
!        endif

        enddo !j loop

       if(ic>=nsubbin) then
       wavlittilda(i+1)=wavlittilda(i)+totdlam
!       print*,'dlam',i+nsizebig+1,(wavelit(i+1)-wavelit(i)),totdlam
!       print*,'ic=',ic,'wavlitt',i+nsizebig+1,wavelit(i),wavlittilda(i)
       endif
       

!       if(totdlam.ne.0) then
!       print*,'ic',ic,'nntot',nntot
!       endif
       enddo !i loop
     
!       do i=1,nsizebig+1
!       if(wavbigtilda(i)>0 .or. wavbigtilda(i+1)>0) then
!      ! print*,'wavbig',wavebig(i),wavbigtilda(i)
!       endif
!       enddo

!       do i=1,nsizelit+1
!       if(wavlittilda(i)>0 .or. wavlittilda(i+1)>0) then
!       print*,'wavlit',wavelit(i),wavlittilda(i)
!       endif
!       enddo
!clsa
!clsa

      endif 
            
      return
end subroutine dfintervals_filter


!----------------------------------------------------------------------

!----------------------------------------------------------------------
