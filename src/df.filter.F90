      SUBROUTINE filter_profile()
      use types
      use dfcomm
      implicit none
      double precision :: dxlambda_phi,x0,x1,bb
      integer :: i 
      integer :: if_data;

       CALL filterdata()

!     xlambda_phi(1)=299.6-4.4
!     xlambda_phi(Nlambda_phi)=299.6+4.4

      xlambda_phi(1)=xlambda_min 
      xlambda_phi(Nlambda_phi)=xlambda_max 
      dxlambda_phi=(xlambda_phi(Nlambda_phi)-xlambda_phi(1)) &
                   /(Nlambda_phi-1)
       
      do i=2,Nlambda_phi
      xlambda_phi(i)=xlambda_phi(i-1)+dxlambda_phi
      enddo      
      xlambda_phi(Nlambda_phi)=xlambda_max 

!!!lazy method for linear filter only!!!!!!!!!!
!      x0=xlambda_phi(1)
!      x1=xlambda_phi(Nlambda_phi)
!      bb=(xlambda_phi(Nlambda_phi)-xlambda_phi(1))
!      do i=1,Nlambda_phi
!      filter(i)=filter(i)*(x1-xlambda_phi(i))/bb
!      print*,'i=',i,'x0=',x0,'filter=',filter(i)
!      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      END SUBROUTINE filter_profile
!
!
       subroutine filterdata()
       use types
       use dfcomm
       implicit none
       double precision :: norm, ff
       integer :: i

! parameters of the atomic data file atom_2996_1202.dat (3rd line)
!WLREF           = 2996.000 ; reference wavelength in Angström
!WL1             = -44.00   ; lower boundary of the wavelength range [WLREF-WL1 , WLREF+WL2]
!WL2             = +44.00   ; upper boundary of the wavelength range [WLREF-WL1 , WLREF+WL2]
!NWL             = 4401     ; number of wavelength points


!
      open(unit = 102, file = 'filterdata.txt', form = 'formatted', &
     &      status = 'old')
      read(102,*) Nlambda_phi 
      print*,Nlambda_phi
      read(102,*) xlambda_min
      print*,xlambda_min
      read(102,*) xlambda_max
      print*,xlambda_max
      read(102,*) norm
      print*,norm

      allocate(xlambda_phi(Nlambda_phi))
      allocate(filter(Nlambda_phi))

      do i=1,Nlambda_phi 
      read(102,*) ff
      filter(i)=ff
      enddo

!1111  format('e14.7') 

      do i=1,Nlambda_phi
      filter(i)=filter(i)*norm
      print*,'filter:',i,filter(i)
      enddo
      close(102)

      return 
      END SUBROUTINE filterdata
