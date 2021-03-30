      subroutine josh

      use atlcomm
      implicit none
*
*.... 1994 FEB - MODIFIED TO SOLVE THE TRANSFER EQUATION ONLY IF
*                TAUNU(1) .LT. 0.2
*
*----------------------------- COMMONS ---------------------------------
*
      include 'common.abtotb'
      include 'common.freqbl'
      include 'common.ifblkk'
      include 'common.meudbl'
      include 'common.musblk'
      include 'common.optotb'
      include 'common.parblk'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.taushj'
      include 'common.tcorrb'
      include 'common.tempbl'
*
*--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  extau(maxd), slope
      integer           j, mu, nc

!---- test
      real(kind=8) ao(maxd), bo(maxd), co(maxd)
      real(kind=8) b2ct(maxd), b2ct1(maxd), ctwo(maxd)


*
*------------------------------ EXTERNALS ------------------------------
*
      external integ, meudon, parcoe, iparcoe 
*
*------------------------------- EXECUTION -----------------------------
*
      ndepth = nrhox
*
      do j = 1, nrhox
         abtot(j) = acont(j)  + aline(j) + sigmac(j) + sigmal(j)
         alpha(j) = (sigmac(j) + sigmal(j)) / abtot(j)
         eps(j)   = (1.0d0 - alpha(j)) * bnu(j)
      end do
*
      call integ (rhox, abtot, taunu, (nrhox), (abtot(1) * rhox(1)))
! ---- in case of bad interpolation/ integration, use a linear interpolation at high tau
! --- first find the depth point where stuff starts to go wrong, then make a linear interpolation from there

      if (minval(taunu) .lt. 0.0d0) then 
        print*, 'Linear interpolation needs to be used'
        do j = nrhox, 2, -1 
          if (taunu(j) .lt. taunu(j-1)) nc = j-1
        end do
        print*, 'Linear interpolation starting from nc = ', nc 
        call integmod(rhox,abtot,taunu,nrhox, (abtot(1)*rhox(1)), nc)

      endif 
*
      if(taunu(1) .ge. 0.2d0) then
*
         do j = 1, nrhox
            jnu(j)   = 0.0d0
            hnu(j)   = 0.0d0
            knu(j)   = 0.0d0
            snu(j)   = 0.0d0
            jmins(j) = 0.0d0
         end do
*
         do mu = 1, nmu
            surfi(mu) = 0.0d0
         end do
*
      else
*
*.... THE OPACITY AT J+1/2 = 0.5d0 * (abtot(j)+abtot(j+1))
*.... THE MASS WITHIN THE LAYER FROM J TO J+1 = rhox(j+1)-rhox(j)
*.... THEREFORE, dtaunu IS BETWEEN J AND J+1, CENTERED AT J+1/2
*
         do j = 1, nrhox-1
            dtaunu(j) = 0.5d0 * (abtot(j) + abtot(j+1)) * 
     &                          (rhox(j+1) - rhox(j))
         end do
*
         dtaunu(nrhox) = taunu(nrhox) - taunu(nrhox-1)
         pl = bnu(nrhox)
         ddp = (pl - bnu(nrhox-1)) / dtaunu(nrhox)
         call meudon
*
         do j = 1, nrhox
            jnu(j)   = ej(j)
            hnu(j)   = eh(j)
            knu(j)   = ek(j)
            snu(j)   = eps(j) + alpha(j) * ej(j)
            jmins(j) = jnu(j) - snu(j)
         end do
*
         hnu(nrhox) = ddp / 3.0d0
*
         if (ifsurf .eq. 2) then
*
*.... COMPUTE THE SURFACE SPECIFIC INTENSITY.  
*.... FOR NOW, THIS IS THE ANALYTIC INTEGRATION FROM TAU(J) TO TAU(J+1)
*.... OF THE SOURCE FUNCTION REPRESENTED BY (A +B*TAU). 
*.... LATER I WILL DO A BETTER INTEGRATION METHOD.
*
        call iparcoe(snu, taunu, ao,bo,co, nrhox )
!
        do j = 1, nrhox 
           ctwo(j) = co(j) * 2.0
           b2ct(j) = bo(j) + ctwo(j) * taunu(j)
        end do
! 
        do j = 1, nrhox -1 
           b2ct1(j) = bo(j) + ctwo(j) * taunu(j+1)
        end do
!
!
            do mu = 1, nmu
!           
               do j = 1, nrhox
                  extau(j)  = exp(-taunu(j) / angle(mu))
               end do
!
               surfi(mu) = 0.0d0
!---------------------------------------
!
               do j = 1, nrhox -1
               surfi(mu) = surfi(mu) + extau(j) * (snu(j) + (b2ct(j) +
     &                              ctwo(j) * angle(mu)) * angle(mu))-
     &                              extau(j+1) * (snu(j+1) +
     &                              (b2ct1(j) + ctwo(j) * angle(mu)) *
     &                              angle(mu))
               end do 

               surfi(mu) = surfi(mu) + extau(nrhox) * (snu(nrhox) +
     &                           (b2ct(nrhox) + ctwo(nrhox) *
     &                            angle(mu)) * angle(mu))

!--------------------------------
!               do j = 1, nrhox-1
!                  slope = (snu(j+1) - snu(j)) / (taunu(j+1) - taunu(j))
!                  surfi(mu) = surfi(mu) +
!     &               extau(j) * (snu(j) + slope * angle(mu)) -
!     &               extau(j+1) * 
!     &                  (snu(j) + slope * (angle(mu) + taunu(j+1) -
!     &                   taunu(j)))
!               end do
!---------------------------------------------------------------------!
            end do
*
         end if
*
      end if
*
      end
*
*********** E N D   O F   S U B R O U T I N E   J O S H ****************
*
      subroutine meudon
      use atlcomm

      implicit none
*
*----------------------------- COMMONS ---------------------------------
*
      include 'common.abtotb'
      include 'common.meudbl'
      include 'common.tcorrb'
      include 'common.taushj'
*
*------------------------------ CONSTANT -------------------------------
*
      integer     nmu
      parameter ( nmu = 3 )
!      parameter ( nmu = 4 )
!      parameter ( nmu = 8 )
*
*-------------------------- LOCAL VARIABLES ----------------------------
*
      double precision  aa(nmu), bb(nmu, nmu), cc(nmu), 
     &                  dd(nmu, nmu, maxd), dt0, dtinv(maxd), factor, 
     &                  jmu(maxd, nmu), jgrid, mu(nmu), mu2(nmu), 
     &                  psi(nmu, maxd), q(nmu), 
     &                  wt(nmu), wtmu(nmu), wtmu2(nmu)
      integer           i, id, j
      logical           first
      save              first, jmu, mu, mu2, wt, wtmu, wtmu2
*
*----------------------------- INITIALIZATION --------------------------
*
      data first / .true. /
*
!--- for 3
      data mu   / 0.8872983346,  0.5000000000,  0.1127016654 /
!
!--- for 4 
!      data mu   / 0.9305681558,  0.6699905218,
!     &            0.3300094782,  0.0694318442 /
!------for 8 
!      data mu   / 0.9801449282,  0.8983332387,  0.7627662050,
!     &            0.5917173212,  0.4082826788,  0.2372337950,
!     &            0.1016667613,  0.0198550718 /
!

!---- for 3: 
      data wt / 0.2777777778, 0.44444444444, 0.27777777778 /
! ----
!--- for 4 : 
!      data wt / 0.1739274226, 0.3260725774,
!     &          0.3260725774, 0.1739274226  /
!--- for 8 : 
!      data wt / 0.0506142681,  0.1111905172,  0.1568533229,
!     &          0.1813418917,  0.1813418917,  0.1568533229,
!     &          0.1111905172,  0.050614281 /

*------------------------------- EXECUTION -----------------------------
*
      if (first) then
         first = .false.
*
         do i = 1, nmu
            mu2(i)   = mu(i) * mu(i)
            wtmu(i)  = wt(i) * mu(i)
            wtmu2(i) = wt(i) * mu2(i)
         end do
*
      end if
*
      do j = 1, ndepth - 1
         dtinv(j) = 1.0d0 / dtaunu(j)
      end do
*
*.... UPPER BOUNDARY.  SECOND ORDER CONDITION
*
      do i = 1, nmu
         factor = 0.5d0 * dtaunu(1) / mu(i)
         aa(i) = 0.0d0
         cc(i) = mu(i) * dtinv(1)
         q(i) = eps(1) * factor
*
*.... THE OFF-DIAGONAL TERMS FROM THE SCATTERING
*
         do j = 1, nmu
            bb(i, j) = -factor * alpha(1) * wt(j)
         end do
*
*.... THE DIAGONAL TERMS
*
         bb(i, i) = bb(i, i) + 1.0d0 + cc(i) + factor
      end do
*
      call matinv (bb, (nmu) )
*
      do i = 1, nmu
         psi(i, 1) = 0.0d0
*
         do j = 1, nmu
            psi(i, 1) = psi(i, 1) + bb(i, j) * q(j)
            dd(i, j, 1) = bb(i, j) * cc(j)
         end do
*
      end do
*      
*.... NORMAL DEPTH POINTS
*
      do id = 2, ndepth - 1
         dt0  =  2.0d0 / (dtaunu(id-1) + dtaunu(id))
*
         do i = 1, nmu
            aa(i) = mu2(i) * dtinv(id-1) * dt0
            cc(i) = mu2(i) * dtinv(id) * dt0
            q(i) = eps(id) + aa(i) * psi(i, id-1)
*
*.... SCATTERING TERMS
*
            do j = 1, nmu
               bb(i, j) = -wt(j) * alpha(id)
            end do
*
*.... ON-DIAGONAL TERMS
*
            bb(i, i) = bb(i, i) + 1.0d0 + aa(i) + cc(i)
         end do
*
         do i = 1, nmu
*
            do j = 1, nmu
               bb(i, j) = bb(i, j) - aa(i) * dd(i, j, id-1)
            end do
*
         end do
*
         call matinv (bb, (nmu) )
*
         do i = 1, nmu
            psi(i, id) = 0.0d0
*
            do j = 1, nmu
               psi(i, id) = psi(i, id) + bb(i, j) * q(j)
               dd(i, j, id) = bb(i, j) * cc(j)
            end do
*
         end do
*
      end do
*
*.... LOWER BOUNDARY.  SECOND ORDER CONDITION
*
      do i = 1, nmu
         factor = 0.5d0 * dtaunu(ndepth-1) / mu(i)
         aa(i) = mu(i) * dtinv(ndepth-1)
         q(i) = pl + mu(i) * ddp + factor * eps(ndepth) + 
     &          aa(i) * psi(i, ndepth - 1)
*
*.... SCATTERING TERMS
*
         do j = 1, nmu
            bb(i, j) = -factor * alpha(ndepth) * wt(j)
         end do
*
*.... ON-DIAGONAL TERMS
*
         bb(i, i) = bb(i, i) + 1.0d0 + aa(i) + factor
      end do
*
      do i = 1, nmu
*
         do j = 1, nmu
            bb(i, j) = bb(i, j) - aa(i) * dd(i, j, ndepth - 1)
         end do
*
      end do
*
      call matinv (bb, (nmu) )
*
      do i = 1, nmu
         psi(i, ndepth) = 0.0d0
*
         do j = 1, nmu
            psi(i, ndepth) = psi(i, ndepth) + bb(i, j) * q(j)
         end do
*
         jmu(ndepth, i) = psi(i, ndepth)
      end do
*
      ej(ndepth) = 0.0d0
*
      do i = 1, nmu
         ej(ndepth) = ej(ndepth) + wt(i) * jmu(ndepth, i)
      end do
*
*.... BACK SUBSTITUTION
*
      do id = ndepth - 1, 1, -1
*
         do i = 1, nmu
            jmu(id, i) = 0.0d0
*
            do j = 1, nmu
               jmu(id, i) = jmu(id, i) + dd(i, j, id) * jmu(id+1, j)
            end do
*
            jmu(id, i) = jmu(id, i) + psi(i, id)
         end do
*
      end do
*
*.... CHANGE THE TOP ONE THE WAY THAT DIMITRI SHOWS IN FRH, P 373
*
      ej(1) = 0.0d0
      eh(1) = 0.0d0
      ek(1) = 0.0d0
*
      do i = 1, nmu
         jgrid = jmu(1, i) / (1.0d0 + 0.5d0 * dtaunu(1) / mu(i))
         ej(1) = ej(1) + wt(i)    * jmu(1, i)
         eh(1) = eh(1) + wtmu(i)  * jgrid
         ek(1) = ek(1) + wtmu2(i) * jmu(1, i)
      end do
*
      do id = 2, ndepth - 1
         ej(id) = 0.0d0
         eh(id) = 0.0d0
         ek(id) = 0.0d0
         dt0  =  0.5d0 * (dtaunu(id-1) + dtaunu(id))
*
         do i = 1, nmu
            eh(id) = eh(id) + wtmu2(i) * (jmu(id,i) - jmu(id-1,i)) /
     &                                   dtaunu(id-1)
            ej(id) = ej(id) + wt(i)    * jmu(id, i)
            ek(id) = ek(id) + wtmu2(i) * jmu(id, i)
         end do
*
      end do
*
      end
*
*********** E N D   O F   S U B R O U T I N E   M E U D O N ************
*
      subroutine matinv(a, n)
      implicit none
*
*--------------------------- DUMMY ARGUMENTS ---------------------------
*
      integer n
      double precision a(n, n)
*
*--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  div, sum
      integer           i, j, l, k, k0
*
*------------------------------- EXECUTION -----------------------------
*
*.... LU DECOMPOSITION ,  STARTING WITH L
*
      do i  =  2,  n
*
         do j  =  1,  i - 1
            div = a(j, j)
            sum  =  0.0d0
*
            if ( j - 1 .ge. 1 ) then
*
               do l  =  1,  j - 1
                  sum  =  sum + a(i, l) * a(l, j)
               end do
*
            end if
*
            a(i, j) = (a(i, j) - sum) / div
         end do
*
         do j  =  i, n
            sum  =  0.0d0
*
            do l  =  1,  i - 1
               sum  =  sum + a(i, l) * a(l, j)
            end do
*
            a(i, j)  =  a(i, j) - sum
         end do
*
      end do
*
*.... INVERSION OF L
*
      do i  =  n, 2, -1
*
         if ( i - 1 .ge. 1 )  then
*
            do j  =  i - 1, 1, -1
               sum  =  0.0d0
*
               if ( j + 1 .le. i - 1) then
*
                  do k  =  j + 1,  i - 1
                     sum  =  sum + a(i, k) * a(k, j)
                  end do
*
               end if
*
               a(i, j)  =  - a(i, j) - sum
            end do
*
         end if
*
      end do
*
*.... U INVERSION
*
      do i  =  n, 1, -1
         div = a(i, i)
*
         if ( i + 1 .le. n ) then
*
            do j  =  n, i + 1, -1
               sum  =  0.0d0
*
               do k  =  i + 1,  j
                  sum  =  sum + a(i, k) * a(k, j)
               end do
*
               a(i, j) = -sum / div
            end do
*
         end if
*
         a(i, i)  =  1.0d0 / a(i, i)
      end do
*
*.... MULTIPLICATION OF U INVERSE AND L INVERSE
*
      do i  =  1,  n
*
         do j  =  1,  n
            k0  =  max0(i, j)
*
            if ( k0 .eq. j ) then
               sum  =  a(i, k0)
               k0  =  k0 + 1
*
            else
               sum  =  0.0d0
            end if
*
            do k  =  k0,  n
               sum  =  sum + a(i, k)*a(k, j)
            end do
*
            a(i, j)  =  sum
         end do
*
      end do
*
      end
*
*********** E N D   O F   S U B R O U T I N E   M A T I N V ************
*
