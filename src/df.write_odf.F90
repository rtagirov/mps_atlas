subroutine write_odf(iv, comm, myrank, ntemp, ODFfile)
!     gets the ODFs computed for each T
!     separates ODFs in BIG and LITTLE ODFs
!     new: writes out the whole big and little to .nc files, gives the frequency grid!
      use types
      use dfcomm

      implicit none

      include 'common.rhoxbl'
      include 'common.turbpr'
      include 'common.tempbl'
      include 'common.stateb'
      include 'common.ifopbl'
      include 'common.junkbl'

      integer, intent(in) :: iv, comm, myrank
      integer, intent(in) :: ntemp
      character(30), intent(in) :: ODFfile

#ifdef MPI
      include 'mpif.h'
      integer mtag, mstatus(MPI_STATUS_SIZE)  
 
#endif


      integer    :: i,  ip, inu, it, itape, istep
      integer    :: filen
      integer    :: ibin, isubbin
      integer    :: ncid, ier, idfnum
      integer    :: ncountl, ncountb


      real(kind=8), allocatable   :: logt(:), logp(:), freqgrid(:)
      character(1) velf
      character(6) tempstr


!     first get the right frequency bins!

      if (ifkbin ) then

       ibin= int(nsizebig+nsizelit,kind(ibin))
       isubbin = 12
       allocate(freqgrid(ibin+2))
! -- pass on the frequency grid
       freqgrid(1:nsizebig+1) = wavebig
       freqgrid(nsizebig+2:ibin+2) = wavelit

      else

       ibin = int(nsizebig, kind(ibin))
       isubbin = int(nsubbin,kind(isubbin))

       allocate(freqgrid(ibin+1))

       freqgrid(1:nsizebig+1) = wavebig

      end if

!--- need the t-p gird in log10

       allocate(logp(nrhox))
       allocate(logt(ntemp))

       do i=1,nrhox

        logp(i) = log10(p(i))
       end do

       do i=1, ntemp
        logt(i) = log10(tsave(i))
       end do


! ---- get the right subset of iodfsteps!

#ifdef MPI


     ncountb = int(nsizebig,kind(ncountb))*ntemp*nrhox*isubbin


       do inu = 1, nsizebig
        do it = 1, ntemp
         do ip = 1, nrhox
           do  i = 1, isubbin
             iodfsendb(i,inu, ip, it) = iodfstep(i,inu, ip, it, iv)
           end do
         end do
        end do
       end do

      ncountl = int(nsizelit,kind(ncountl))*ntemp*nrhox*isubbin
       do inu = 1, nsizelit
        do it = 1, ntemp
         do ip = 1, nrhox
           do  i = 1, isubbin
             iodfsendl(i,inu, ip, it) = iodfstep(i,nsizebig+inu, ip, it, iv)
           end do
         end do
        end do
       end do



!  --- send information accross all cores

      call MPI_allreduce(iodfsendb, iodfrecvb, ncountb, MPI_integer2, MPI_sum, comm, ier)

      call MPI_allreduce(iodfsendl, iodfrecvl, ncountl, MPI_integer2, MPI_sum, comm, ier)

! --- if filters are used the sbwiths are changed later during calculations, they need to be send to rank 0
      if (ifilter .eq. 1 ) then 
        ncountl = int((nsizelit+nsizebig),kind(ncountl))*isubbin
        mtag = 2345
        if (myrank .eq. 1 ) then
          sb_weight_send =sbwith 
          call mpi_send(sb_weight_send, ncountl, MPI_DOUBLE_PRECISION, 0, mtag, comm, ier)        
        else if (myrank .eq. 0) then 
          call mpi_recv(sb_weight_recv, ncountl, MPI_DOUBLE_PRECISION, mpi_any_source, mpi_any_tag,  comm, mstatus, ier) 
          sbwith = sb_weight_recv
        end if 
      endif 

#else
       do inu = 1, nsizebig
        do it = 1, ntemp
         do ip = 1, nrhox
          do i = 1, isubbin
            iodfrecvb(i,inu, ip, it) = iodfstep(i,inu, ip, it, iv)
          end do
         end do
        end do
       end do

       do inu = 1, nsizelit
        do it = 1, ntemp
         do ip = 1, nrhox
          do i = 1, isubbin
            iodfrecvl(i ,inu, ip, it) = iodfstep(i, nsizebig+inu, ip, it, iv)
          end do
         end do
        end do
       end do

#endif



!--- need the turb-vel as character for the filename

       velf = char(int(iv) +48)

!--- only one core writes all the data out---!
! --- note in the non MPI version the only core has been given myrank = 0 !

if (myrank .eq. 0) then

      if (ifkbin) then

! --- old fashioned odfs in binary are still produced

        open(unit=2, file = 'p00big'//trim(velf)//'.bdf', form='unformatted',status='new')
        open(unit=3, file = trim(ODFfile) , form='unformatted',status='new')

        do inu=1,nsizebig
          do it=1,ntemp
           write(2)((iodfrecvb(istep, inu, ip,it),istep=1,12),ip=1,25)
         end do
       end do

       do inu=1,nsizelit
          do it=1,ntemp
           write(3)((iodfrecvl(istep, inu, ip,it),istep=1,12),ip=1,25)
         end do
       end do

      close(unit=2)
      close(unit=3)

!old standart binary finished

      else

      if(isubbin .gt. 2) then 
       open(unit=2, file = 'p00big'//trim(velf)//'.bdf', form='formatted',status='new')
         do inu=1,nsizebig
           write(2,*) freqgrid(inu), freqgrid(inu+1)
           do ip = 1, nrhox
            do it=1,ntemp
              write(2,4)(iodfrecvb(istep, inu, ip,it),istep=1,isubbin)
4       FORMAT(1X,12(I6,1X))
            end do
           end do
          end do
       close(unit=2)

      end if 


!       open(unit=2, file = 'p00big'//trim(velf)//'.bdf', form='formatted',status='new')
!         do inu=1,nsizebig
!          do it=1,ntemp
!           write(2,*)((iodfrecvb(istep, inu, ip,it),istep=1,isubbin),ip=1,nrhox )
!          end do
!          end do
!       close(unit=2)



      call CreateODF(0,ncid,trim(ODFfile),title,ibin,isubbin,nrhox,ntemp ,ivt(iv),ier)

      call WriteODF(ncid,0,comm,iodfrecvb,freqgrid,logt,logp,sbwith, ibin, nrhox, ntemp, isubbin, ier)

      call CloseNetCDF( 0, ncid ,ier )



      end if
end if

!    clean up the arrays

      deallocate(freqgrid)
      deallocate(logt)
      deallocate(logp)

      return
end subroutine

