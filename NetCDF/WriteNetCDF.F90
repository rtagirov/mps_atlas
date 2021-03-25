!--------------------------------------------------------------
! WriteNetCDF:
!         Write output to netcdf file
!.....................................................................
! 
!--------------------------------------------------------------
!---------------------------------------------------------------

subroutine WriteContinum(  ncid, myrank, nproc, comm,   &
                        wledge, freqedge, cmedge, idmol, momass, freqset, &
                        temp, pressure, t4, rho, xne, xnatom, vturb, xnfh24, xnfhe4, &
                        xnfh4, contall, contsca, xnfpel4, dopmax, dopple4, &
                        nrhox, ntemp, nufreq, ier) 

!
      implicit none
!
      include 'netcdf.inc'

#ifdef MPI
      include 'mpif.h'
#else
! ---careful right now it is not parallel
#endif

      integer          istart1,  istart2(2), istart3(3), istart4(4)
      integer          icount1,  icount2(2), icount3(3), icount4(4)
      integer           i,j,k,l,m 



      integer     ier
      integer     iter
      integer     myrank 
      integer     comm
      integer     nproc ! tells the proccessor or the number of T_i value!
      integer     ncid
      

      integer     nrhox, ntemp, nufreq
      integer     nufreq3  


      double precision  wledge(nufreq/3), freqedge(nufreq/3), cmedge(nufreq/3)
      double precision  idmol(60), momass(60)
      double precision  freqset(nufreq)
      double precision  temp(ntemp), pressure(nrhox)

      double precision  t4(nrhox), rho(nrhox)
      double precision  xne(nrhox), xnatom(nrhox)
      double precision  vturb(nrhox), xnfh24(nrhox)
      double precision  xnfhe4(nrhox,  2), xnfh4(nrhox, 2)
!
      double precision  contall(nufreq, nrhox), contsca(nufreq, nrhox)
      double precision  xnfpel4(6, 99, nrhox), dopmax(nufreq/3, 6, 99 ) 
      double precision  dopple4(6, 99, nrhox) 


      integer           status

      integer           varid

      integer           ncount
      

      character(4)       :: varname

      ier = 0
      nufreq3= nufreq/3

! Non-MPI version.




!---- if iproc it equal to one then write out quantities that do not change 
!     during the calculations such as wledge, freqset and so on


      if (nproc .eq. 1) then 

      print*, 'Writing fixed quantities '
  
      
      status = nf_inq_varid( ncid, 'mID', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, 60, idmol(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_inq_varid( ncid, 'mmass', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, 60, momass(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)


      status = nf_inq_varid( ncid, 'frqed', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, nufreq3, freqedge(1)  )
      if ( status .ne. NF_NOERR) call handle_error(status)


      status = nf_inq_varid( ncid, 'wled', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, nufreq3, wledge(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)


      status = nf_inq_varid( ncid, 'cmed', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, nufreq3, cmedge(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)


      status = nf_inq_varid( ncid, 'frqset', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, nufreq, freqset(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)


      status = nf_inq_varid( ncid, 'pressure', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, nrhox, pressure(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_inq_varid( ncid, 'temp', varid )
      if ( status .ne. NF_NOERR) call handle_error(status)

      status = nf_put_vara_double(ncid,varid,1, ntemp, temp(1) )
      if ( status .ne. NF_NOERR) call handle_error(status)


     endif



     istart1 = 1+ (nrhox*(nproc-1))
     icount1 = nrhox

          status = nf_inq_varid( ncid, 't4', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart1, icount1, t4 )
          if ( status .ne. NF_NOERR) call handle_error(status)
 

          status = nf_inq_varid( ncid, 'rho', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart1, icount1, rho )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_inq_varid( ncid, 'Xne', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart1, icount1, xne )
          if ( status .ne. NF_NOERR) call handle_error(status)


          status = nf_inq_varid( ncid, 'xnatom', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart1, icount1, xnatom )
          if ( status .ne. NF_NOERR) call handle_error(status)



          status = nf_inq_varid( ncid, 'vturb', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart1, icount1, vturb )
          if ( status .ne. NF_NOERR) call handle_error(status)


          status = nf_inq_varid( ncid, 'xnfh24', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart1, icount1, xnfh24 )
          if ( status .ne. NF_NOERR) call handle_error(status)

          do i=1, 2

       istart2(1) = 1+(nrhox*(nproc-1))
       istart2(2) = i 
       icount2(1) = nrhox 
       icount2(2) = 1 


          status = nf_inq_varid( ncid, 'xnfh4', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart2, icount2, xnfh4(1,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)


                    status = nf_inq_varid( ncid, 'xnfhe4', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart2, icount2, xnfhe4(1,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)
         end do


        do i = 1, nrhox
       istart2(1) = 1
       istart2(2) = i+(nrhox*(nproc-1))
       icount2(1) = nufreq
       icount2(2) = 1

          status = nf_inq_varid( ncid, 'contall', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)
          status = nf_put_vara_double(ncid,varid,istart2, icount2, contall(1,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)


          status = nf_inq_varid( ncid, 'contsca', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)
          status = nf_put_vara_double(ncid,varid,istart2, icount2, contsca(1,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

        end do





       do i = 1, nrhox
         do j = 1, 99
          
         istart3(1) = 1
         istart3(2) = j
         istart3(3) = i+ (nrhox*(nproc-1))
         icount3(1) = 6
         icount3(2) = 1
         icount3(3) = 1


          status = nf_inq_varid( ncid, 'xnfpel4', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)
          status = nf_put_vara_double(ncid,varid,istart3, icount3, xnfpel4(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_inq_varid( ncid, 'dopple', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)
          status = nf_put_vara_double(ncid,varid,istart3, icount3, dopple4(1,j,i) )
          if ( status .ne. NF_NOERR) call handle_error(status)


         end do 
      end do 



          do j = 1, 99
             do l= 1, 6 
         istart3(1) = 1
         istart3(2) = l
         istart3(3) = j +(99*(nproc-1)) 

         icount3(1) = nufreq3 
         icount3(2) = 1
         icount3(3) = 1


          status = nf_inq_varid( ncid, 'dopmax', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

          status = nf_put_vara_double(ncid,varid,istart3, icount3, dopmax(1,l,j) )
          if ( status .ne. NF_NOERR) call handle_error(status)

           end do
        end do

!............................................................................

      return 

end  subroutine  WriteContinum


subroutine WriteODF(ncid, myrank, comm,   &
                      odf, freqgrid, tabt, tabp, sbwt, nbin,  np, nt, nsubbin,  ier)

      implicit none
      include 'netcdf.inc'

#ifdef MPI
      include 'mpif.h'
#else
!---no MPi implementation
#endif
      integer          istart1,  istart2(2), istart3(3), istart4(4)
      integer          icount1,  icount2(2), icount3(3), icount4(4)
      integer          i,j,k,l,m


      integer     ier
      integer     iter
      integer     myrank
      integer     comm
      integer     nproc ! tells the proccessor or the number of T_i value!
      integer     ncid

      integer     np, nt
      integer     nbin, nsubbin 
      integer(kind=2)     odf(nsubbin,nbin,np,nt) 
      double precision freqgrid(nbin+1), tabt(nt), tabp(np), sbwt(nsubbin, nbin)

      integer           status
      integer           varid
      integer           ncount


      ! --- write the boundaries of the subbins
          status = nf_inq_varid( ncid, 'subbin', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)
 
          do j = 1, nbin 

          istart2(1)  = 1
          istart2(2)  = j

          icount2(1)  = nsubbin
          icount2(2)  = 1

          status = nf_put_vara_double(ncid,varid,istart2, icount2, sbwt(1,j))
          if ( status .ne. NF_NOERR) call handle_error(status)

          end do 
      ! --- write the boundraries of the bins

          status = nf_inq_varid( ncid, 'FreqG', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)


          istart1 = 1
          icount1 = nbin+1 

          status = nf_put_vara_double(ncid,varid,istart1, icount1, freqgrid(1))
          if ( status .ne. NF_NOERR) call handle_error(status)

      ! --- write the temperature tab

          status = nf_inq_varid( ncid, 'T', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)


          istart1 = 1
          icount1 = nt 

          status = nf_put_vara_double(ncid,varid,istart1, icount1, tabt(1))
          if ( status .ne. NF_NOERR) call handle_error(status)




      ! --- write the pressure tab

          status = nf_inq_varid( ncid, 'P', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)


          istart1 = 1
          icount1 = np

          status = nf_put_vara_double(ncid,varid,istart1, icount1, tabp(1))
          if ( status .ne. NF_NOERR) call handle_error(status)




      ! --- now write the actual ODF 
          status = nf_inq_varid( ncid, 'ODF', varid )
          if ( status .ne. NF_NOERR) call handle_error(status)

        if (nsubbin .eq. 1) then

        do j = 1, nsubbin 
             do l= 1, np  
               do  m = 1, nt
         istart4(1) = j
         istart4(2) = 1 
         istart4(3) = l
         istart4(4) = m

         icount4(1) = 1
         icount4(2) = nbin
         icount4(3) = 1
         icount4(4) = 1


          status = nf_put_vara_int2(ncid,varid,istart4, icount4, odf(j,1,l,m))
          if ( status .ne. NF_NOERR) call handle_error(status)
             end do 
           end do
        end do

        else 

        do j = 1, nbin
             do l= 1, np
               do  m = 1, nt
         istart4(1) = 1
         istart4(2) = j
         istart4(3) = l
         istart4(4) = m

         icount4(1) = nsubbin 
         icount4(2) = 1
         icount4(3) = 1
         icount4(4) = 1


          status = nf_put_vara_int2(ncid,varid,istart4, icount4, odf(1,j,l,m))
          if ( status .ne. NF_NOERR) call handle_error(status)
             end do
           end do
        end do


        endif 







      return
end subroutine WriteODF 
