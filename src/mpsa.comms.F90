#ifdef MPI
module comms
   use types
   use atlcomm

   implicit none
   include 'mpif.h'
   integer :: ierr, comm, rank, sizee, error
   integer :: status




   Contains
!-----------------------------------------------
subroutine initialise_comm
  implicit none 

!       Initialize MPI.
!
      call MPI_Init ( ierr )
      comm = MPI_COMM_WORLD
!
!  Determine this process's ID.
!
      call MPI_Comm_rank ( comm, rank, ierr )
!
!  Find out how many processes are active.
!
      call MPI_Comm_size ( comm, sizee,  ierr )





end subroutine initialise_comm


!------------------------------------------------
   subroutine begin_comm
   implicit none 

   include 'common.rhoxbl'
   include 'common.freqbl'
   include 'common.teffbl'
   include 'common.fluxbl'
   include 'common.turbpr'
   include 'common.steplg' 

!   integer, intent(in) :: comm, rank
   integer:: i, j 
!   integer :: ierr,  error
!pradk0, pzero, teff, flux, grav, glog, vturb
 
! --- initilise, and allocate arrays for communications
! --- allocate only dim(nrhox, 8) 

! --- allocate only dim(8) for the scalars
   nrhox = krhox 
   call alloc_comm(nrhox) 

! --- assigne the scalars: 
   if (rank .eq. 0) then 
      sendscalar(1) = teff
      sendscalar(2) = grav
      sendscalar(3) = glog
      sendscalar(4) = flux
      sendscalar(5) = pzero
      sendscalar(6) = pradk0
      sendscalar(7) = vturb(1) 
   end if 

! --- broad cast:
 
   call MPI_bcast(sendscalar, 8, MPI_double_precision, 0, comm, ierr)

   if (rank .ne. 0) then 
      teff =    sendscalar(1)
      grav =    sendscalar(2) 
      glog =    sendscalar(3)
      flux =    sendscalar(4)
      pzero =   sendscalar(5)
      pradk0 =  sendscalar(6) 
      vturb(1) =   sendscalar(7)

      do i=2,nrhox
        vturb(i) = vturb(1)
      end do 
   end if 

end subroutine begin_comm



!-------------------------------------------------
subroutine startiteration !(rank,comm)
! t, rhox, p, tau, accrad, pradk, prad,  ptotal,( pturb ?) 
     implicit none

     include 'common.rhoxbl'
     include 'common.freqbl'
     include 'common.stateb'
     include 'common.tempbl' 
     include 'common.abross'



!     integer, intent(in) :: comm, rank
     integer:: i, j, ncount 
!     integer :: ierr,  error

     ncount = nrhox *9 
     if (rank .eq. 0 ) then 


       sdarrays(1:nrhox,1) = rhox(1:nrhox)
       sdarrays(1:nrhox,2) = t(1:nrhox)
       sdarrays(1:nrhox,3) = p(1:nrhox)
       sdarrays(1:nrhox,4) = accrad(1:nrhox)

       sdarrays(1:nrhox,5) = ptotal(1:nrhox)
       sdarrays(1:nrhox,6) = prad(1:nrhox)
       sdarrays(1:nrhox,7) = pradk(1:nrhox)
       sdarrays(1:nrhox,8) = tauros(1:nrhox)
       sdarrays(1:nrhox,9) = xne(1:nrhox)
     endif 


!    use  bcast to get all the atmosphere information to all cores :
    call MPI_bcast(sdarrays,ncount,MPI_double_precision,0,comm,ierr) 

     if (rank .ne. 0) then 

      rhox(1:nrhox)   =  sdarrays(1:nrhox,1) 
         t(1:nrhox)   =  sdarrays(1:nrhox,2) 
       p(1:nrhox)     =  sdarrays(1:nrhox,3)
       accrad(1:nrhox)=  sdarrays(1:nrhox,4) 

      ptotal(1:nrhox) =  sdarrays(1:nrhox,5) 
      prad(1:nrhox)   =  sdarrays(1:nrhox,6) 
      pradk(1:nrhox)  =  sdarrays(1:nrhox,7) 
      tauros(1:nrhox)    =  sdarrays(1:nrhox,8) 
         xne(1:nrhox) =  sdarrays(1:nrhox,9)



     endif 


end subroutine startiteration


subroutine finiteration!(rank,comm)


     implicit none
     include 'common.rhoxbl'
     include 'common.freqbl'
     include 'common.abross'
     include 'common.fluxbl'


!     integer, intent(in) :: comm, rank
     integer:: i, j, ncount
!     integer :: ierr,  error
     real(kind=8) :: sends, recs
     ncount = nrhox *9
       sdarrays(1:nrhox,1) = abross(1:nrhox)
       sdarrays(1:nrhox,2) = raden(1:nrhox)
       sdarrays(1:nrhox,3) = flxrad(1:nrhox)
       sdarrays(1:nrhox,4) = accrad(1:nrhox)
!
       sdarrays(1:nrhox,5) = hflux(1:nrhox)
       sdarrays(1:nrhox,6) = rdabh (1:nrhox)
       sdarrays(1:nrhox,7) = rjmins(1:nrhox)
       sdarrays(1:nrhox,8) = rdiagj (1:nrhox)

!    use MPI_allreduce :
    call MPI_allreduce(sdarrays, recarrays,  ncount, & 
  &            MPI_double_precision,  MPI_sum, comm, ierr)


       abross(1:nrhox) = recarrays(1:nrhox,1)
       raden(1:nrhox)  = recarrays(1:nrhox,2) 
       flxrad(1:nrhox) = recarrays(1:nrhox,3)
       accrad(1:nrhox) = recarrays(1:nrhox,4) 


       hflux(1:nrhox) = recarrays(1:nrhox,5)
       rdabh(1:nrhox) = recarrays(1:nrhox,6)
       rjmins(1:nrhox) =recarrays(1:nrhox,7)
       rdiagj(1:nrhox) = recarrays(1:nrhox,8)

!--- scalar pradk0

     sends = pradk0

     call MPI_allreduce(sends, recs,1, MPI_double_precision, MPI_sum, & 
  &                      comm, ierr) 
     pradk0 = recs


end subroutine finiteration 
end module 
#endif 
