!--------------------------------------------------------------
! CreateNetCDF:
!         Create a netcdf file, dimension the arrays and 
!         put it into data mode
!--------------------------------------------------------------

      subroutine CreateContinum ( myrank, ncid, name,  nrhox, ntemp, &
                                nufreq,  ier)

      implicit none 
     
      include 'netcdf.inc'

      integer:: ier, ncid 
      integer:: nrhox, ntemp, nufreq, nufreq3 
      integer:: atmodim, numdim , twodim , mdim, m99dim 
      integer:: status
      integer:: myrank ! just for later parallisation purposes


      integer:: numions, onedim, numatoms, nugrid, freqdim
      integer:: dims1, dims2(2), dims3(3)
      integer::  varid, sixtydim
      character(len=*)     name

      nufreq3=nufreq/3

      ier = 0
      if ( myrank .ne. 0 ) return

      status = nf_create( trim(name), OR(NF_CLOBBER,NF_64BIT_OFFSET) , ncid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'nrhox', nrhox, atmodim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'ntemp', ntemp, numdim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'nufreq', nufreq, freqdim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'merge',nrhox*ntemp , mdim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
    
      status = nf_def_dim( ncid, 'merge1', 99*ntemp , m99dim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)
  
!------------------------------------------------------------------
      status = nf_def_dim( ncid, 'six', 6, numions )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'ninenine', 99, numatoms )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'nugrid', nufreq3, nugrid  )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'one', 1, onedim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'two', 2, twodim )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      status = nf_def_dim( ncid, 'sixty', 60, sixtydim )
      if ( status .ne. NF_NOERR ) call  handle_error(status) 
!---------------------------------------------------------------!

!           CREATE VARIABLES ------------------------------------!

       
        dims1 = sixtydim

        status = nf_def_var( ncid, 'mID', NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)


        status = nf_def_var( ncid, 'mmass', NF_DOUBLE, 1,dims1,varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)



       dims1 = nugrid 

       status = nf_def_var( ncid, 'frqed',NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)


       status = nf_def_var( ncid, 'wled',NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       status = nf_def_var( ncid, 'cmed',NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)


!---------   variable sizes ------------------------------------------

       dims1 = freqdim

      status = nf_def_var( ncid, 'frqset',NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims1 = atmodim

      status = nf_def_var( ncid, 'pressure',NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims1 = numdim 
      
      status = nf_def_var( ncid, 'temp',NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims1 = mdim

        status = nf_def_var( ncid, 't4', NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status = nf_def_var( ncid, 'rho', NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status = nf_def_var( ncid, 'Xne', NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status = nf_def_var( ncid, 'xnatom', NF_DOUBLE, 1, dims1,varid)
        if ( status .ne. NF_NOERR ) call  handle_error(status)
  
        status = nf_def_var( ncid, 'vturb', NF_DOUBLE, 1, dims1, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status = nf_def_var( ncid, 'xnfh24', NF_DOUBLE, 1, dims1,varid)
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims2(1) = mdim
       dims2(2) = twodim  
 
        status = nf_def_var( ncid, 'xnfh4', NF_DOUBLE, 2, dims2, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status = nf_def_var( ncid, 'xnfhe4', NF_DOUBLE, 2, dims2,varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims2(1) = freqdim
       dims2(2) = mdim


        status = nf_def_var( ncid, 'contall', NF_DOUBLE, 2,dims2,varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status = nf_def_var( ncid, 'contsca', NF_DOUBLE, 2,dims2,varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims3(1) = numions
       dims3(2) = numatoms 
       dims3(3) = mdim

        status = nf_def_var( ncid, 'xnfpel4', NF_DOUBLE, 3,dims3,varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

        status =nf_def_var( ncid, 'dopple', NF_DOUBLE, 3, dims3, varid)
               if ( status .ne. NF_NOERR ) call  handle_error(status)

       dims3(1) = nugrid
       dims3(2) = numions
       dims3(3) = m99dim 

        status = nf_def_var( ncid, 'dopmax', NF_DOUBLE, 3,dims3, varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)

! nf_enddef puts named file out of "define" mode
      status = nf_enddef(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)

! nf_sync synchronises disk writes with memory buffers (immediately 
! available after writing)
      status = nf_sync(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      end subroutine CreateContinum
!--------------------------------------------------------------------------------
!
!  Subroutines to generate ODF files!
!--------------------------------------------------------------------------------
subroutine CreateODF (myrank,ncid,name, title, nbin,nsubbin,np,nt,vturb,ier)
!subroutine CreatODF (myrank, ncid, name, nbin, nsubbin, np, nt, vturb, ier)

      implicit none

      include 'netcdf.inc'


      integer:: ier, ncid
      integer, intent(in):: np,   nbin, nsubbin 
      integer, intent(in) :: nt
      integer :: nsubp1
      integer:: dimp, dimt, dimsubbin, dimbin, dimnumf 
      integer:: status
      integer:: myrank 
      integer:: dims4(4), dims2(2)
      integer::  varid
      integer:: len
      character(len=*)     name
      character(len=74)    title
     
      double precision :: vturb

      nsubp1 = nbin+ 1

      ier = 0
      if ( myrank .ne. 0 ) return
!--- create

      status = nf_create( trim(name), OR(NF_CLOBBER,NF_64BIT_OFFSET) , ncid )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

!-- define dimensions

      status = nf_def_dim( ncid, 'np', np, dimp )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      status = nf_def_dim( ncid, 'nt', nt, dimt )
      if ( status .ne. NF_NOERR ) call  handle_error(status)

      status = nf_def_dim( ncid, 'nbins', nbin, dimbin )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      status = nf_def_dim( ncid, 'nsubbins', nsubbin , dimsubbin )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      status = nf_def_dim(ncid, 'numfp', nsubp1, dimnumf )
      if ( status .ne. NF_NOERR) call handle_error(status)

!---------------------------------------------------------------!
      len = 73
      status = nf_put_att_text( ncid, NF_GLOBAL, 'title',     &
                                   len, title  )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


      status = nf_put_att_double( ncid, NF_GLOBAL, 'vturb',     &
                                   NF_DOUBLE, 1, vturb )
      if ( status .ne. NF_NOERR ) call  handle_error(status)


       dims4(1) = dimsubbin 
       dims4(2) = dimbin 
       dims4(3) = dimp 
       dims4(4) = dimt


        status = nf_def_var( ncid, 'ODF', NF_SHORT, 4,dims4, varid )

        status = nf_def_var( ncid, 'FreqG', NF_DOUBLE, 1, dimnumf , varid )
        if ( status .ne. NF_NOERR ) call  handle_error(status)


        status = nf_def_var( ncid, 'P', NF_DOUBLE, 1, dimp, varid) 
        if (status .ne. NF_NOERR ) call handle_error(status)


        status = nf_def_var( ncid, 'T', NF_DOUBLE, 1, dimt, varid) 
        if (status .ne. NF_NOERR ) call handle_error(status) 

        dims2(1) = dimsubbin
        dims2(2) = dimbin

        status = nf_def_var( ncid, 'subbin', NF_DOUBLE, 2, dims2, varid) 
        if (status .ne. NF_NOERR ) call handle_error(status)

! nf_enddef puts named file out of "define" mode
      status = nf_enddef(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)

! nf_sync synchronises disk writes with memory buffers (immediately 
! available after writing)
      status = nf_sync(ncid)
      if ( status .ne. NF_NOERR ) call  handle_error(status)

end subroutine CreateODF


