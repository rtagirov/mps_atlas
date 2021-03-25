subroutine dfintervals
!
!     The intervals big and little  remain  hardcoded for the standard Kurucz bin/sub-bin configuration  
!     The sub-bins  borders for the user defined configuration are read in here
!--------------------------------------------------------------------------------------------------------
      use types
      use dfcomm
      implicit none
      
      include 'common.ifopbl'

      real(kind=dp) :: nzero
      integer :: i, j, litbig, check 

!-----------------------------------------------------------------------------

      call def_binsize()
!   if not Kuruzc grid then wavelit, wavebig and nsizes are defined in def_binsize
! 
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
!- then the end of each bin is used to determine the beginning of the next bin ----------!

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

!--- same procedure again as above for Kurucz little grid
      
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
 
!       ifbin = true --> standard Kurucz bin/sub-bin is used:
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
! --- for user defined sub-bins the input file subbin_info.dat is read and the nsteps are calculated.

      open (unit = 78, status = 'old', file='./INPUT/subbin_info.dat', form = 'formatted')
      read(78,*) nsubbin, check 
       
      do i = 1, nsizebig+ nsizelit 
         read(78,*) (binbound(i,j), j=1, nsubbin) 
         ! here the upper boundaries are read, for example for a bin with 3 subbins:
         ! 0.3, 0.6, 1.0 ---------------------------------------------------------- 

 
          nsteps(1,i)=1
          sbwith(1,i) = binbound(i,1)
        do j = 2, nsubbin 
          sbwith(j,i) = binbound(i,j)-binbound(i,j-1)       
          nsteps(j,i) = (nn(i)*int(1000*binbound(i,j-1)))/1000 +1
        end do 
         nsteps(nsubbin+1,i) = (nn(i)*int(1000*binbound(i,nsubbin)))/1000 +1
 
      end do 

      close(unit =78) 


     end if 
!--- finally the global arrays nbeg, nend which are equivlent to nbig and nlit end & begin are 
!--- assigned. These arrays are used in df.calc together with nsteps to calculated the opacity averages

      litbig = nsizebig+nsizelit

      do i = 1, nsizebig
       nbeg(i) = nbigbeg(i)
       nend(i) = nbigend(i)
      end do 

      do i = nsizebig+1, litbig
       nbeg(i) = nlitbeg(i-nsizebig)
       nend(i) = nlitend(i-nsizebig)
      end do 
      
            
      return
end subroutine dfintervals 






subroutine dfinterhighres( istart, ifin)  
!
!     For the high-resolution opacity tables :
!------------------------------------------------------------------------------
      use types
      use dfcomm
      implicit none

      include 'common.ifopbl'

      real(kind=dp) :: l1, l2, delta2  
      integer :: i, j, litbig, check
      integer(kind=4), intent(in)::  istart, ifin
!-----------------------------------------------------------------------------
! high resolution gird is defined as:    l_n = l_0 *(1 + 1/res)**(n-1)
! lg(l_n) / lg(1 + 1/res) = lg(l_0)/lg(1+1/res) + n -1
!-------------------------------------------------------------------------
!     we need to define wavebig boundaries! 
        l1 = wlbeg*(1.0d0 + 1.0d0/resolu)**(istart-2) 

      do i = istart, ifin+1
        j = i - istart +1

        l2 = wlbeg*(1.0d0 + 1.0d0/resolu)**(i-1)
        delta2 = (l2-l1)*0.5d0
        wavebig(j) = l1+delta2
        sbwith(1,j) = 1.0d0
        l1 = l2 
      end do 




      return
end subroutine dfinterhighres 


!----------------------------------------------------------------------

!----------------------------------------------------------------------
