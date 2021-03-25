subroutine dfcalc(iteff,j,iv,t,p,vt)

! This subroutine, takes the high-resolution opacities (globally stored in buffer
! calculates the log_10 () and multiplies by 1000
! then it sorts each bin and calculates the geometric mean in each of the sub-bins
! the geometric mean is multiplied by sqrt(10)
!------------------------------------------------------

  use types
  use dfcomm
  implicit none
  include 'common.ifopbl'


  integer,  intent(in) :: iteff, j, iv 
  real(kind=dp),  intent(in):: t,p, vt
!-------------------------------------------------------
! local variables
!-------------------------------------------------------
  integer(kind=8) ::  ntot,  nbuff, idf, nbegin, nstop
  integer :: i,l,jj, l1, nu,  istep  
  integer(kind=4) :: oldbufi,  ibuffer(lenbuff), ibuffsort(223938)
  integer(kind=4) ::  nnsub, idfnum 
  integer(kind=8) :: isum, isum2
  real(kind=dp):: tlog, plog !careful the same are globale in common.tempbl that is not used here!
!--- for arithmetic mean
  real(kind=8) ::dbsum

!-----------------------------------------------------

     tlog=log10(t)
     plog=log10(p)
     length = 3507859    ! length is total amount of frequency points on high resolution 

! takes high-res opacity,  calculates the log_10 () and multiplies by 1000
   do  nbuff=1,length
     ibuffer(nbuff)=nint(log10(buffer(nbuff))*1000.)   
   end do 

! Loop over the bin grids (for now big and little in one go)
  nnsub = nsubbin 

  do idf=1, nsizebig+nsizelit     
    insteps = 0 

! --- begin, end and amount of points in the bin with index idf
      
    nbegin=int(nbeg(idf))
    nstop=int(nend(idf))
    ntot=int(nstop-nbegin+1)   
!--- take only the part of the high-resolution which corresponds to the bin idf

    do i=nbegin,nstop
      ibuffsort(i-nbegin+1)=ibuffer(i)   
    end do  

!  -   sort ibuffsort, ibuffsort(100000) is used as buffer-space
    call bisort(ibuffsort,ntot,ibuffsort(100000))  
!    	 					 

    if ( ameanODF .eq. 1) then 
! ----- get arithmetic mean in each sub-bin
      do istep=1, nnsub
        dbsum = 0.0d0

! -- nsteps give the points of the sub-bin boundaries

        do i=nsteps(istep,idf),nsteps(istep+1,idf)-1
            dbsum = dbsum +10.0**(0.001d0*dble(ibuffsort(i)))
        end do

         dbsum = dbsum / (dble(nsteps(istep+1,idf)-nsteps(istep,idf)))
         insteps(istep) = nint(log10(dbsum)*1000.)

      end do


    else 
!   --- get geometric mean in each sub-bin
 
      do istep=1, nnsub 
        isum=0
! -- nsteps give the points of the sub-bin boundaries

        do i=nsteps(istep,idf),nsteps(istep+1,idf)-1       
           isum=isum+int(ibuffsort(i), 8)                       
        end do
    
        isum2 = (isum+int((nsteps(istep+1,idf)-nsteps(istep,idf)),8)/2)/ & 
                   & int((nsteps(istep+1,idf)-nsteps(istep,idf)),8)
        insteps(istep) = int(isum2,2) 
      end do 

! --- isteps is the averaged opacity in each sub-bin

     end if 
      
     do i = 1, 20 
      if(iv.eq.i) then 
       do istep = 1, nnsub
        iodfstep(istep, idf, j,iteff, i)= insteps(istep)
       end do
       end if 

     end do    
  end do 

  return
end subroutine


subroutine dfcalchighres(iteff,j,iv,t,p,vt, nbegin)
!--- this is for high-resolution opacity tables, not for ODF tables:
!
! This subroutine takes the high-resolution opacities 
! and stores it into the right memory of the iodfstep
!------------------------------------------------------

  use types
  use dfcomm
  implicit none

  integer,  intent(in) :: iteff, j, iv
  integer,  intent(in) :: nbegin
  real(kind=dp),  intent(in):: t,p, vt


!-------------------------------------------------------
! local variables
!-------------------------------------------------------
  integer(kind=8) ::    nbuff, idf
  integer :: i 
  integer(kind=4) :: ibuffer(lenbuff)


!-------------------------------------
! takes the high-resolution opacities (globally stored in buffer
! calculates the log_10 () and multiplies by 1000

   do idf = 1, nsizebig+1   
     nbuff = nbegin+(idf-1)
     ibuffer(idf)=nint(log10(buffer(nbuff))*1000.)

     do i = 1, 20
      if(iv.eq.i) then
        iodfstep(1, idf, j,iteff, i)= ibuffer(idf)
      end if
     end do
  end do

  return
end subroutine



