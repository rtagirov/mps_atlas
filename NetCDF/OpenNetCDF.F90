!--------------------------------------------------------------
! OpenNetCDF:
!         Open a netcdf file, read the attributes and
!         array dimensions.
!--------------------------------------------------------------

      subroutine OpenContinum(ncid, name, nrhox, ntemp, nufreq, ier)

!
      implicit none
!
      include 'netcdf.inc'

      integer          ier

      character(len=*)   name
      integer::        ncid, nrhox, ntemp, nufreq
      integer::        status, varid
      integer::        nrhoxid, ntempid, nufreqid

      ier = -1
      ncid = -1
! NF_OPEN reads named file, NF_WRITE implies file is writeable
! ncid is returned NetCDF file ID.

      status =  NF_OPEN( trim(name), NF_WRITE, ncid )
      if ( status .ne. NF_NOERR ) then
         print*,'ERROR: Unable to open file ', name
         call  handle_error(status)
         return
      end if

 

!.......Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid
 
      status = nf_inq_dimid( ncid, 'nrhox', nrhoxid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'ntemp', ntempid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'nufreq', nufreqid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

!.......Read the actual dimensions 
! nf_inq_dimlen returns length of dimension
!  
      status = nf_inq_dimlen( ncid, nrhoxid, nrhox )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, ntempid, ntemp )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, nufreqid, nufreq )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      ier = 0

      return 

      end subroutine  OpenContinum

#ifdef MPI 
!--------------------------------------------------------------
! OpenNetCDF:
!         Open a netcdf file, read the attributes and
!         array dimensions.
!--------------------------------------------------------------

subroutine OpenContpar (ncid, myrank,name, nrhox, ntemp, nufreq, ier)

!
      implicit none
!
      include 'netcdf.inc'

      integer          ier

      character(len=*)    name
      integer::        ncid, nrhox, ntemp, nufreq
      integer::        status, varid
      integer::        nrhoxid, ntempid, nufreqid
      integer:: myrank

      ier = -1
      ncid = -1
! NF_OPEN reads named file, NF_WRITE implies file is writeable
! ncid is returned NetCDF file ID.

      status =  NF_OPEN( trim(name), NF_WRITE, ncid )
      if ( status .ne. NF_NOERR ) then
         print*,'ERROR: Unable to open file ', name
         call  handle_error(status)
         return
      end if

 
 if (myrank .eq. 0 ) then 
!.......Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid
 
      status = nf_inq_dimid( ncid, 'nrhox', nrhoxid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'ntemp', ntempid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'nufreq', nufreqid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

!.......Read the actual dimensions 
! nf_inq_dimlen returns length of dimension
!  
      status = nf_inq_dimlen( ncid, nrhoxid, nrhox )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, ntempid, ntemp )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, nufreqid, nufreq )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

  end if 

      ier = 0

      return 

      end subroutine  OpenContpar


#endif
 

!      subroutine OpenODF (ncid, name, title,  nbin, nsubbin,  np, nt, vturb, ier)
     subroutine OpenODF (ncid, name, nbin, nsubbin, np, nt, vturb, ier)

!
      implicit none
!
      include 'netcdf.inc'

      character(len=*)    name
!      character(len=74)   title
      integer::        ier
      integer::        ncid, np, nt,  nbin, nsubbin 
      integer::        status, varid
      integer::        npid, ntid,  nbinid, nsbinid, nnfid 
      double precision :: vturb
      ier = -1
      ncid = -1
! NF_OPEN reads named file, NF_WRITE implies file is writeable
! ncid is returned NetCDF file ID.

      status =  NF_OPEN( trim(name), NF_WRITE, ncid )
      if ( status .ne. NF_NOERR ) then
         print*,'ERROR: Unable to open file ', name
         call  handle_error(status)
         return
      end if

!      status = nf_get_att_text( ncid,NF_GLOBAL,'title', title )
!      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_get_att_double( ncid, NF_GLOBAL, 'vturb', vturb)
      if ( status .ne. NF_NOERR ) call  handle_error(status)


!.......Read the dimension id's
! nf_inq_dimid returns dimension ID (integer) of "string" in ncid

      status = nf_inq_dimid( ncid, 'np', npid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'nt', ntid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'nbins', nbinid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimid( ncid, 'nsubbins', nsbinid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_inq_dimid( ncid, 'numfp', nnfid )
      if (status .ne. NF_NOERR ) call handle_error(status)


!.......Read the actual dimensions 
! nf_inq_dimlen returns length of dimension
!  
      status = nf_inq_dimlen( ncid, npid, np )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, ntid, nt )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, nbinid, nbin )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
      status = nf_inq_dimlen( ncid, nsbinid, nsubbin )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      ier = 0

      return

      end subroutine  OpenODF


