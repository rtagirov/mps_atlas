!--------------------------------------------------------------
! This routine uses the input "status" (integer) to determine 
! whether or not an error has occurred. NF_NOERR is an integer 
! which implies success. NF_STRERROR returns a string indicating 
! that there has been an error
!--------------------------------------------------------------

subroutine handle_error(status)
      implicit none

      integer  status 
      include 'netcdf.inc'

      if ( status .ne. NF_NOERR ) then
        print*, NF_STRERROR(STATUS)
      end if
end  subroutine handle_error

