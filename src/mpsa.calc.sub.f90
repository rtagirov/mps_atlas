 subroutine calc_se 

    use population
    


    implicit none
    include 'common.rhoxbl'
    include 'common.xnfblk'
    include 'common.xnfpbl'
    include 'common.stateb'
    include 'common.ifblkk'
    include 'common.ionsbl'

    integer j
    double precision  part(maxd,6)


               if ( recalxne .or. ifpres )  then
                 call pops(0.0d0, 1, xne)
               end if



               call pops( 1.01d0, 11, xnfph)
               call pops( 2.02d0, 11, xnfphe)
               call pops( 6.03d0, 11, xnfpc)
               call pops( 7.04d0, 11, xnfpn)
               call pops( 8.05d0, 11, xnfpo)
               call pops(10.05d0, 11, xnfpne)
               call pops(12.01d0, 11, xnfpmg)
               call pops(13.00d0, 11, xnfpal)
               call pops(14.01d0, 11, xnfpsi)
               call pops(20.01d0, 11, xnfpca)
               call pops(26.00d0, 11, xnfpfe)
!*
               if (.not. ifmol) then
                  call pops( 1.01d0, 12, xnfh)
                  call pops( 2.02d0, 12, xnfhe)
                  call pops( 6.05d0, 12, xnfc)
                  call pops( 7.05d0, 12, xnfn)
                  call pops( 8.05d0, 12, xnfo)
                  call pops(10.05d0, 12, xnfne)
                  call pops(12.05d0, 12, xnfmg)
                  call pops(14.05d0, 12, xnfsi)
                  call pops(16.05d0, 12, xnfs)
                  call pops(26.04d0, 12, xnffe)

!
               else
! --- number density / partition function

                  call pops(106.00d0, 11, xnfpch)
                  call pops(108.00d0, 11, xnfpoh)
!
!                  call pops(107.00d0, 11, xnfpnh)
!----
! -- number density:
                  call pops(107.00d0, 12, xnfpnh) 
!                  call pops(106.00d0, 12, xnfpch)
!                  call pops(108.00d0, 12, xnfpoh) 


!
!...  THE POPS WILL NOT RETURN NUMBER DENSITIES WHEN MOLECULES ARE ON
!..  SO WE COMPUTE NUMBER DENSITIES/PART FUNCTIONS  AND PART FUNCTIONS
!
                  call pops( 6.05d0, 11, xnfc)
                  call pops( 7.05d0, 11, xnfn)
                  call pops( 8.05d0, 11, xnfo)
                  call pops(10.05d0, 11, xnfne)
                  call pops(12.05d0, 11, xnfmg)
                  call pops(14.05d0, 11, xnfsi)
                  call pops(16.05d0, 11, xnfs)
                  call pops(26.04d0, 11, xnffe)
!
                  do j = 1,nrhox
!
                     call pfsaha(j, 1, 1, 3, part)
                     xnfh(j, 1) = xnfph(j, 1) * part(j, 1)
                     xnfh(j, 2) = xnfph(j, 2)
!
                     call pfsaha(j, 2, 2, 13, part)
                     xnfhe(j, 1) = xnfphe(j, 1) * part(j, 1)
                     xnfhe(j, 2) = xnfphe(j, 2) * part(j, 2)
                     xnfhe(j, 3) = xnfphe(j, 3)
!
                     call pfsaha(j, 6, 6, 13, part)
                     xnfc(j, 1) = xnfc(j, 1) * part(j, 1)
                     xnfc(j, 2) = xnfc(j, 2) * part(j, 2)
                     xnfc(j, 3) = xnfc(j, 3) * part(j, 3)
                     xnfc(j, 4) = xnfc(j, 4) * part(j, 4)
                     xnfc(j, 5) = xnfc(j, 5) * part(j, 5)
                     xnfc(j, 6) = xnfc(j, 6) * part(j, 6)
!
                     call pfsaha(j, 7, 6, 13, part)
                     xnfn(j, 1) = xnfn(j, 1) * part(j, 1)
                     xnfn(j, 2) = xnfn(j, 2) * part(j, 2)
                     xnfn(j, 3) = xnfn(j, 3) * part(j, 3)
                     xnfn(j, 4) = xnfn(j, 4) * part(j, 4)
                     xnfn(j, 5) = xnfn(j, 5) * part(j, 5)
                     xnfn(j, 6) = xnfn(j, 6) * part(j, 6)
!
                     call pfsaha(j, 8, 6, 13, part)
                     xnfo(j, 1) = xnfo(j, 1) * part(j, 1)
                     xnfo(j, 2) = xnfo(j, 2) * part(j, 2)
                     xnfo(j, 3) = xnfo(j, 3) * part(j, 3)
                     xnfo(j, 4) = xnfo(j, 4) * part(j, 4)
                     xnfo(j, 5) = xnfo(j, 5) * part(j, 5)
                     xnfo(j, 6) = xnfo(j, 6) * part(j, 6)

                     call pfsaha(j, 10, 6, 13, part)
                     xnfne(j, 1) = xnfne(j, 1) * part(j, 1)
                     xnfne(j, 2) = xnfne(j, 2) * part(j, 2)
                     xnfne(j, 3) = xnfne(j, 3) * part(j, 3)
                     xnfne(j, 4) = xnfne(j, 4) * part(j, 4)
                     xnfne(j, 5) = xnfne(j, 5) * part(j, 5)
                     xnfne(j, 6) = xnfne(j, 6) * part(j, 6)

                     call pfsaha(j, 12, 6, 13, part)
                     xnfmg(j, 1) = xnfmg(j, 1) * part(j, 1)
                     xnfmg(j, 2) = xnfmg(j, 2) * part(j, 2)
                     xnfmg(j, 3) = xnfmg(j, 3) * part(j, 3)
                     xnfmg(j, 4) = xnfmg(j, 4) * part(j, 4)
                     xnfmg(j, 5) = xnfmg(j, 5) * part(j, 5)
                     xnfmg(j, 6) = xnfmg(j, 6) * part(j, 6)

                     call pfsaha(j, 14, 6, 13, part)
                     xnfsi(j, 1) = xnfsi(j, 1) * part(j, 1)
                     xnfsi(j, 2) = xnfsi(j, 2) * part(j, 2)
                     xnfsi(j, 3) = xnfsi(j, 3) * part(j, 3)
                     xnfsi(j, 4) = xnfsi(j, 4) * part(j, 4)
                     xnfsi(j, 5) = xnfsi(j, 5) * part(j, 5)
                     xnfsi(j, 6) = xnfsi(j, 6) * part(j, 6)


!
                     call pfsaha(j, 16, 6, 13, part)
                     xnfs(j, 1) = xnfs(j, 1) * part(j, 1)
                     xnfs(j, 2) = xnfs(j, 2) * part(j, 2)
                     xnfs(j, 3) = xnfs(j, 3) * part(j, 3)
                     xnfs(j, 4) = xnfs(j, 4) * part(j, 4)
                     xnfs(j, 5) = xnfs(j, 5) * part(j, 5)
                     xnfs(j, 6) = xnfs(j, 6) * part(j, 6)

                     call pfsaha(j, 26, 5, 13, part)
                     xnffe(j, 1) = xnffe(j, 1) * part(j, 1)
                     xnffe(j, 2) = xnffe(j, 2) * part(j, 2)
                     xnffe(j, 3) = xnffe(j, 3) * part(j, 3)
                     xnffe(j, 4) = xnffe(j, 4) * part(j, 4)
                     xnffe(j, 5) = xnffe(j, 5) * part(j, 5)
                  end do
!
               end if








 end subroutine calc_se

 subroutine calc_press_hydro

   use atlcomm

  implicit none
  include 'common.rhoxbl'
  include 'common.turbpr'
  include 'common.teffbl'
  include 'common.stateb'

  integer i, j


!.... INTEGRATE THE EQUATION OF HYDROSTATIC EQUILIBRIUM
!
                  pzero = pcon + pradk0 + pturb(1)
!
                  do j = 1, nrhox
                     ptotal(j) = grav * rhox(j) + pzero
                     p(j) = grav * rhox(j) - prad(j)-pturb(j)-pcon
                   if (p(j) .le. 0.0d0) then
                     p(j)= grav*rhox(j)-(prad(j)-prad(1))-pturb(j)-pcon
                       if (p(j) .le. 0.0d0) then
                          write (16, 9002)     & 
                          ' pressure le 0.0d0 at depth', j, & 
                          ' pzero = ', pzero, & 
                          'j', 'p(j)', 'accrad(j)', 'prad(j)',  &    
                          (i, p(i), accrad(i), prad(i), i = 1, j)
!
9002                    format (a, i5, a, 1pe15.4, / ,     & 
                          4x, a, 8x, a, 9x, a, 7x, a /     &  
                         (i5, 3e15.4))
                        stop
                      end if
                     end if
!
                  end do
!-----------------------------------------------------------*

 end subroutine calc_press_hydro




