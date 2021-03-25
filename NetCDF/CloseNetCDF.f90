!--------------------------------------------------------------
! Closes off NetCDF file if on the master processor (myrank=0)
!--------------------------------------------------------------
subroutine CloseNetCDF( myrank, ncid ,ier )
      implicit none

      include 'netcdf.inc'

      integer         ncid
      integer         myrank
      integer         ier
      integer         status


      ier = 0 
!      if ( myrank .ne. 0 ) return 

! Synchronises buffers etc...
      status = nf_sync(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)

! Closes off named NetCDF file
      status = nf_close(ncid)
      if ( status .ne. NF_NOERR ) call handle_error(status)


      return
end subroutine CloseNetCDF
