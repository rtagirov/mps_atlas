      module population

       use types
       use atlcomm

       implicit none


      contains


      subroutine pops (code, mode, number)


      implicit none
!
!----------------------------- COMMONS ---------------------------------
!
!      include 'common.sizebl'
      include 'common.ifblkk'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xabund'
*
*--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision, intent(in) ::  code
      double precision, intent(out):: number(maxd, *)
      integer, intent(in)::  mode
!
!------------------------------- EXTERNALS ----------------------------
!
!      external molec, nelect, nmolec, pfsaha
!
!--------------------------- LOCAL VARIABLES --------------------------
!
      integer  j, itemp1, ion, iz, nion, nnnn
      save     itemp1
*
*--------------------------- INITIALIZATION ---------------------------
*
      data itemp1 / 0 /
*
*------------------------------- EXECUTION -----------------------------
*
      if (ifmol) then
        
        if ( recalxne ) then     
          if(itemp .ne. itemp1) call nmolec((mode))
        else 
          if( ifpres .and. (itemp .ne. itemp1)) then 
             call nmolec((mode))
          end if 
        end if

         itemp1 = itemp
         if (code .gt. 0.0d0) call molec((code), (mode), number)
*
      else

        if (recalxne)  then 
         if ((itemp .ne. itemp1)) call nelect
        else
         if (ifpres .and. (itemp .ne. itemp1)) call nelect 
        end if 


         itemp1 = itemp
*
         if (code .gt. 0.0d0) then
*
            if (code .lt. 100.0d0) then
               iz   = code
               nion = (code - dble(iz)) * 100.0d0 + 1.5
*
               do j = 1, nrhox
                  call pfsaha((j), (iz), (nion), (mode), number)
*
*.... PFSAHA RETURNS IONIZATION FRACTIONS OR  
*....                IONIZATION FRACTIONS / PARTITION FUNCTIONS 
*.... SO CONVERT TO NUMBER DENSITIES
*
                  nnnn = nion
                  if ( mode .lt. 10) nnnn = 1
*
                  do ion = 1, nnnn
                     number(j, ion) = number(j, ion) * xnatom(j) *
     &                                xabund(iz)
                  end do
*
               end do
*
            else
               write (16, *) '1molecules off'
               print *, '1molecules off'
               stop
            end if
*
         end if
*
      end if
*
      end subroutine 
*
*********** E N D   O F   S U B R O U T I N E   P O P S ****************
*
      subroutine molec (codout, mode, number)
     
      implicit none
*
*----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.elembl'
      include 'common.ifblkk'
      include 'common.ifequa'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xabund'
      include 'common.xnmolb'
*
*--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision, intent(in) ::  codout
      double precision, intent(out):: number(maxd, *)
      integer, intent(in) ::     mode
*
*------------------------------- EXTERNALS ----------------------------
*
!      external pfsaha
*
*--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  c
      integer           i, id, ion, j, jmol, nn
      logical           iread, match, sec400
      save              iread
*
*--------------------------- INITIALIZATION ---------------------------
*
      data iread / .false. /
*
*------------------------------- EXECUTION -----------------------------
*
! exclude this as we do not read xnmol (modi 01)
!-----------------------------------------------------------
!      if ((.not. iread) .and. (.not. ifpres)) then
*
!         read (15, 151) nummol
!151      format (i5)
*
!         do jmol = 1, nummol
!            read (15, 152) code(jmol)
!152         format (f20.2)
!            read (15, 153) (xnmol(j, jmol), j = 1, nrhox)
!153         format (1p8d10.3)
!            write (16, 154) ' xnmol for code', code(jmol), 
!     &                      (xnmol(j, jmol), j=1,nrhox)
!154         format (a, f20.2 / (1p8d10.3))
!         end do
*
!         read  (15, 158)
!         read  (15, 158) (xnatom(j), j = 1, nrhox)
!         read  (15, 158)
!
!         read  (15, 158) (rho(j), j = 1, nrhox)
!         read  (15, 158)
!         read  (15, 158) (xne(j), j = 1, nrhox)
!158      format (1p8d10.3)
!         write (16, 159) ' xnatom', (xnatom(j), j = 1, nrhox)
!
!         write (16, 159) ' rho', (rho(j), j = 1, nrhox)
!         write (16, 159) ' xne', (xne(j), j = 1, nrhox)
!159      format (a / (1p8d10.3))
!         iread = .true.
!      end if
*
! -- for molecules ...
      if (codout .ge. 100.0d0) then
         jmol = 1
*
         do while (jmol .le. nummol .and. code(jmol) .ne. codout)
            jmol = jmol + 1
         end do
*
         if (jmol .gt. nummol) then
            write (16, 202) ' WARNING from subroutine molec:',
     &                      ' no data for molecule', codout,
     &                      ' - setting its number densities to 0.0'
202         format (a / 5x, a, f20.2 / 5x, a)
*
            do j = 1, nrhox
               number(j, 1) = 0.0d0
            end do
*
         else

           if (mode .eq. 12) then 
             do j =1, nrhox
              number(j,1) = xnmolc(j,jmol)
             end do 
           else 
            do j = 1, nrhox
               number(j, 1) = xnmol(j, jmol)
            end do
           end if  

*
         end if
! -- for atoms ...
      else
         c  = codout
*
         if (mode .eq. 11) then
            nn = (c - dble( int( c ))) * 100.0d0 + 1.5
*
         else
            nn = 1
         end if
*
         i = 0
         sec400 = .false.
! open while (1) 
         do while (i .lt. nn .and. .not. sec400)
            i = i + 1
            ion = nn - i + 1
            jmol = 1
            match = .false.
*
            do while (.not. match .and. jmol .le. nummol)
*
               if (abs(code(jmol) - c) .lt. 0.001d0) then
                  match = .true.
*
               else
                  jmol = jmol + 1
               end if
*
            end do
*
            if (match) then
*
               do j = 1, nrhox
                  number(j, ion) = xnmol(j, jmol)
               end do
*
            else
               id    = codout
               jmol  = 1
               match = .false.
*
               do while (.not. match .and. jmol .le. nummol)
*
                  if(int(code(jmol)) .eq. id) then
                     match = .true.
*
                  else
                     jmol = jmol + 1
                  end if
*
               end do
*
               if (match) then
*
                  do j = 1, nrhox
                     number(j, ion) = 0.0d0
                  end do
*
               else
                  sec400 = .true.
*
               end if
*
            end if
            if(.not. sec400) c = c - 0.01d0
         end do
! ---close while (1)
*
! -----
         if(sec400) then
            ion = (codout - dble(id)) * 100.0d0 + 1.5d0
*
            if(mode .eq. 1) then
               nn = 1
*
            else
               nn  = ion
            end if
*
            do j = 1, nrhox
               call pfsaha ((j), (id), (ion), (mode), number)
*
               do i = 1, nn
                  number(j, i) = number(j, i) * xnatom(j) * xabund(id)
               end do
*
            end do
*
         end if
! ------ *
      end if
*
      end subroutine 
*
*********** E N D   O F   S U B R O U T I N E   M O L E C **************
*
      subroutine nelect

      implicit none
*
*----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.constb'
*
      include 'common.elembl'
      include 'common.edenbl'
      include 'common.iterbl'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xabund'
*
*------------------------------ CONSTANTS ----------------------------
*
      double precision  factre, factrm, factrp
      integer           c4, c5, c12, nz
      parameter (c4 = 4, c5 = 5, c12 = 12, nz = 10, 
     &           factre = 1.0d0 / 0.999d0,
     &           factrm = 0.999d0 / 1.001d0,
     &           factrp = 1.001d0)
*
*------------------------------- EXTERNALS ----------------------------
*
!      external pfsaha
*
*--------------------------- LOCAL VARIABLES --------------------------
*
      character         ee(nz)*2
      double precision  e(nz), eion(7), elec(maxd), error, 
     &                  pfmin(13), pfplus(13), x(nz), x1, xnenew, 
     &                  xntot, xnumb(maxd, 3)
      integer           i, ion, iz, j, l, nelem, nelemz(nz), nion,
     &                  nionz(nz)
      logical           mask(nz), more
      save              nelemz, nionz
*
*--------------------------- INITIALIZATION ---------------------------
*
!                   H, He, C, Na, Mg, Al, Si, K, Ca, Fe
      data nelemz / 1, 2, 6, 11, 12, 13, 14, 19, 20, 26 /
      data nionz  / 2, 3, 2,  2,  2,  2,  2,  2,  2,  2 /
*
*------------------------------- EXECUTION -----------------------------
*
      if(.not. ifedns) then
*
         do i = 1, nz
            nelem = nelemz(i)
            ee(i) = elem(nelem)
         end do
*
         write (16, 402) 'electron contributions',
     &                  (ee(i), i = 1, nz), (ee(i), i = 1, nz)
402      format ( // 50x, a // 2x, 20(4x, a2) / )

          xne(1) = 0.5d0 * p(1) / tk(1)
      end if
*
      do j = 1, nrhox
         if (j .gt. 1 .and. .not. ifedns) xne(j) = xne(j-1) * 
     &      p(j) / p(j-1)
         xntot     = p(j) / tk(j)
         xnatom(j) = xntot - xne(j)
*
         do i = 1, nz
            mask(i) = .true.
         end do
*
         l = 0
         more = .true.
*
         do while (more .and. l .lt. 200)
            l = l + 1
            xnenew = 0.0d0
*
            do i = 1, nz
*
               if (mask(i)) then
                  iz     = nelemz(i)
                  nion   = nionz(i)
                  call pfsaha ((j), (iz), (nion), c4, elec)
                  e(i)   = elec(j)
                  x(i)   = elec(j) * xnatom(j) * xabund(iz)
                  xnenew = xnenew + x(i)
               end if
*
            end do
*
            xnenew    = 0.5d0 * (xnenew + xne(j))
            error     = abs(xne(j) - xnenew) / xnenew

            xne(j)    = xnenew
            xnatom(j) = xntot - xne(j)
*
            if (error .lt. 5.0d-4) then
               more = .false.
*
            else
*
               if (j .gt. 1) then
                  x1 = 1.0d-5 * xne(j)
                  if (error .lt. 0.05d0) x1 = x1 * 10.0d0
*
                  do i = 1, nz
                     if (x(i) .lt. x1) mask(i) = .false.
                  end do
*
               end if
*
            end if
*
         end do
*
         if (l .ge. 200) then
            write (16, 250) ' xne does not converge for depth', j
	print *, ' xne does not converge for depth', j
	print *, ' error is ', error
250         format (a, i3)
            stop
*
         else
            rho(j) = xnatom(j) * wtmole * amu
*
            if (ifedns) then
               edens(j) = 1.5d0 * xntot * tk(j)
*
               do i = 1, nz
                  iz      = nelemz(i)
                  nion    = nionz(i)
*
                  t(j)    = t(j)    * factrp
                  tk(j)   = tk(j)   * factrp
                  tkev(j) = tkev(j) * factrp
                  call pfsaha((j), (iz), (nion), c5, pfplus)
*
*.... THE FOLLOWING LOOP IS NEEDED TO AVOID AN EQUIVALENCE IN THE
*.... ORIGINAL.  THE VALUES FOR EION ARE RETURNED FROM PFSAHA IN
*.... THE VARIABLE PFPLUS(7+).
*
                  do ion = 1, nion
                     eion(ion) = pfplus(ion + 6)
                  end do
*
                  t(j)    = t(j)    * factrm
                  tk(j)   = tk(j)   * factrm
                  tkev(j) = tkev(j) * factrm
                  call pfsaha((j), (iz), (nion), c5, pfmin)
*
                  t(j)    = t(j)    * factre
                  tk(j)   = tk(j)   * factre
                  tkev(j) = tkev(j) * factre
                  call pfsaha((j), (iz), (nion), c12, xnumb)
*
                  do ion = 1, nion
                     xnumb(j, ion) = xnumb(j, ion) * xnatom(j) * 
     &                               xabund(iz)
                     edens(j) = edens(j) + xnumb(j, ion) * tk(j) * 
     &                          (eion(ion) / tkev(j) + 
     &                          (pfplus(ion) - pfmin(ion)) / 
     &                          (pfplus(ion) + pfmin(ion)) * 1000.0d0)
                  end do
*
               end do
*
               edens(j) = edens(j) / rho(j)
*
            end if 
*
            do i = 1, nz
               x(i) = x(i) / xne(j)
            end do
*
            write (16, 404) j, (x(i), i = 1, nz), (e(i), i = 1, nz)
404         format (i3, 20f6.3)
*
         end if
*
      end do
*
      end subroutine 
*
*********** E N D   O F   S U B R O U T I N E   N E L E C T ************
*
      subroutine nmolec (mode)

      implicit none
*
*--------------------------- DUMMY ARGUMENTS ---------------------------
*
      integer, intent(in) :: mode
*
*----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.edenbl'
      include 'common.elembl'
      include 'common.ifequa'
      include 'common.iterbl'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xabund'
      include 'common.xnmolb'
*
*------------------------------ CONSTANTS ----------------------------
*
      double precision  factre, factrm, factrp
      integer           c3, c5, c12
      parameter (c3 = 3, c5 = 5, c12 = 12,
     &           factre = 1.0d0 / 0.999d0,
     &           factrm = 0.999d0 / 1.001d0,
     &           factrp = 1.001d0 )
*
*------------------------------- EXTERNALS ----------------------------
*
      external solvit
*
*--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  amass, d, deq(maxmeq, maxmeq), eion(7), 
     &                  eq(maxmeq), eqold(maxmeq), equilj(maxmol), 
     &                  frac(maxd, 6), pfm(13), pfmin(maxd), pfp(13), 
     &                  pfplus(maxd), ratio, scale, term, tminus, 
     &                  tplus, x, xab(maxmeq), xn(maxmeq), xn100, 
     &                  xneq, xnsave(maxd, maxmeq), xntot, 
     &                  xnz(maxd, maxmeq), dterm(maxmeq), 
     &                  deqt(maxmeq*maxmeq)
      integer           id, ion, j, jmol, jmol1, jmol10, k,  
     &                  l, locj1, locj2, locm, lock, m, ncomp, nequa1
      integer           kl, k1, kk, mk, nequak
      logical           iferr, klimit


! --- bug fix for solving equilibirum
      integer           ncount 


      save xnsave
*
*------------------------------- EXECUTION -----------------------------
*
      nequa1 = nequa + 1
*
      do k = 2, nequa
         id = idequa(k)
         if(id .lt. 100) xab(k) = max(xabund(id), 1.0d-20)
      end do
*
      if(id .eq. 100) xab(nequa) = 0.0d0
      xntot = p(1) / tk(1)
      xn(1) = xntot * 0.5d0
      if(t(1) .lt. 4000.0d0) xn(1) = xntot

      x = xn(1) * 0.1d0
*
      do k = 2, nequa
         xn(k) = x * xab(k)
      end do
*
      if(id .eq. 100) xn(nequa) = x
      xne(1) = x
!   do over atmospher 
!  ---------------------*
      do j = 1, nrhox

        do ncount = 1 , 3

         xntot = p(j) / tk(j)
*
         if((j .gt. 1) .and. (ncount .eq. 1) ) then

            ratio  = p(j) / p(j-1)
            xne(j) = xne(j-1) * ratio
*
            do k = 1, nequa
               xn(k) = xn(k) * ratio
            end do

         end if
*
         if(ifedns) then
*
            do k = 1, nequa
               xn(k) = xnsave(j, k)
            end do
*
         end if
*
         do jmol = 1, nummol
            ncomp = locj(jmol+1) - locj(jmol)
*
            if(equil(1, jmol) .ne. 0.0d0) then
               ion = (code(jmol) - dble(int(code(jmol)))) * 
     &               100.0d0 + 0.5d0
               equilj(jmol) = 0.0d0
    
               if(t(j) .le. 10000.0d0) then 
                  if(code(jmol).eq.101.)then
                   equilj(jmol)=exp(4.478/tkev(j)-4.64584d1+
     &              (1.63660d-3+(-4.93992d-7+(1.11822d-10+(-1.49567d-14+
     &              (1.06206d-18-3.08720d-23*t(j))*t(j))*t(j))*t(j))*
     &              t(j))*t(j)- 1.5*tlog(j))
                  else
                   equilj(jmol)=exp(equil(1,jmol)/tkev(j)-equil(2,jmol)+
     &(equil(3,jmol)+(-equil(4,jmol)+(equil(5,jmol)+(-equil(6,jmol)+
     & equil(7,jmol)*
     &t(j))*t(j))*t(j))*t(j))*t(j)-1.5*(dble(ncomp-ion-ion-1))*tlog(j))
                  endif
                end if

!               if(t(j) .le. 10000.0d0) equilj(jmol) = 
!     &            exp(equil(1, jmol) / tkev(j) - 
!     &            equil(2, jmol) + (equil(3, jmol) + 
!     &            (-equil(4, jmol) + (equil(5, jmol) - equil(6, jmol) *
!     &            t(j)) * t(j)) * t(j)) * t(j) - 1.5d0 * 
!     &            (dble(ncomp - ion - ion - 1)) * tlog(j))
*
            else 
*
               if(ncomp .le. 1) then
                  equilj(jmol) = 1.0
*
               else
                  id  = code(jmol)
                  ion = ncomp-1
                  call pfsaha((j), (id), (ncomp), c12, frac)
                  equilj(jmol) = frac(j, ncomp) / frac(j, 1) * 
     &                           xne(j) ** ion
               end if
*
            end if
*
         end do

         do k = 1, nequa
            eqold(k) = 0.0d0
         end do
*
*.... SET UP 1ST ORDER EQUATIONS FOR THE CHANGE IN NUMBER DENSITY OF
*.... EACH ELEMENT
*
         iferr = .true.
*
! --- do while it did not convert
!---------------------------------
         do while (iferr)
*
            do l = 1, nequa
               do k = 1, nequa
                  deq(k,l) = 0.0d0
               end do
            end do
!
            eq(1) = -xntot

!
            do k = 2, nequa
               deq(1, k) = 1.0d0
               deq(k, k) = 1.0d0
               deq(k, 1) = -xab(k)
               eq(1)     = eq(1) + xn(k)
               eq(k)     = xn(k) - xab(k) * xn(1)

            end do
!
            if (idequa(nequa) .ge. 100) then
               deq(nequa, nequa) = -1.0d0
               eq(nequa) = -xn(nequa)
            end if


! --- do loop over all molecules !
!--------------------------------*
            do jmol = 1, nummol
               ncomp = locj(jmol + 1) - locj(jmol)
*
               if (ncomp .ne. 1) then
                  term  = equilj(jmol)
                  locj1 = locj(jmol)
                  locj2 = locj(jmol+1) - 1
*
                  do lock = locj1, locj2
                     k = kcomps(lock)
*
                     if (k .eq. nequa1) then
                        term = term / xn(nequa)
*
                     else
                        term = term * xn(k)
                     end if
*
                  end do
*
                  eq(1) = eq(1) + term
*
                  do lock = locj1, locj2
                     k = kcomps(lock)
*
                     if (k .lt. nequa1) then
                        d = term / xn(k)
*
                     else
                        k = nequa
                        d = -term / xn(k)
                     end if
*
                     eq(k)     = eq(k) + term
                     deq(1, k) = deq(1, k) + d
*
                     do locm = locj1, locj2
                        m = kcomps(locm)
                        if (m .eq. nequa1) m = nequa
                        deq(m, k) = deq(m, k) + d
!-----
                     end do
*
                  end do
*
*.... CORRECTION TO CHARGE EQUATION FOR NEGATIVE IONS
*
                  k = kcomps(locj2)
*
                  if (idequa(k) .eq. 100) then
*
                     do lock = locj1, locj2
                        k = kcomps(lock)
                        d = term / xn(k)
                        if (k .eq. nequa) eq(k) = eq(k) - term - term
*
                        do locm = locj1, locj2
                           m = kcomps(locm)
                           if (m .eq. nequa) then 
                             deq(m,k) = deq(m,k)-d-d
!-------- 
                           end if
                        end do
*
                     end do
*
                  end if
*
               end if
*
            end do
!---- do end over all molecules 
!-------------------*

            call solvit (deq, (maxmeq), (nequa), eq)

            iferr = .false.
            scale = 100.0d0
*
            do k = 1, nequa
               if (abs (eq(k) / xn(k)) .gt. 0.001d0) iferr = .true.
               if (eqold(k) * eq(k) .lt. 0.0d0) eq(k) = eq(k) * 0.69d0
               xneq = xn(k) - eq(k)
               xn100 = xn(k) * 0.01d0
*
               if (abs(xneq) .lt. xn100) then
                  xn(k) = xn(k) / scale
                  if (eqold(k) * eq(k) .lt. 0.0d0) scale = sqrt(scale)
*
               else
                  xn100 = xn(k) * 100.0d0
                  xn(k) = abs(xneq)
               end if
*
               eqold(k) = eq(k)
            end do
*
         end do
!--- end do converion loop !

*
         do k = 1, nequa
            xnz(j, k) = xn(k)
         end do
*
         xnatom(j) = xn(1)
         rho(j)    = xnatom(j) * wtmole * 1.660d-24
         if (idequa(nequa) .eq. 100) xne(j) = xn(nequa)
*
         do jmol = 1, nummol
            ncomp = locj (jmol + 1) - locj(jmol)
            xnmol(j, jmol) = equilj(jmol)
            locj1 = locj(jmol)
            locj2 = locj (jmol + 1) - 1
*
            do lock = locj1, locj2
               k = kcomps(lock)
*
               if (k .eq. nequa1) then
                  xnmol(j, jmol) = xnmol(j, jmol) / xn(nequa)
*
               else
                  xnmol(j, jmol) = xnmol(j, jmol) * xn(k)

               end if
*
            end do
            xnmolc(j, jmol) = xnmol(j,jmol) 
         end do
       end do ! ncount loop, to make sure the solution converges correctly 
      end do
!--- loop over the atmosphere is finished
!----------------------------------------*

      if(ifedns) then
*
         do j = 1, nrhox
            edens(j) = 1.5d0 * p(j)
         end do
*
         do jmol = 1, nummol
            ncomp = locj(jmol+1) - locj(jmol)
*
            if(equil(1,jmol) .eq. 0.0d0) then
               id = code(jmol)
*
               do j = 1, nrhox
                  t(j)    = t(j)    * factrp
                  tk(j)   = tk(j)   * factrp
                  tkev(j) = tkev(j) * factrp
                  call pfsaha((j), (id), (ncomp), c5, pfp)
*
*.... THE FOLLOWING LOOP IS NEEDED TO AVOID AN EQUIVALENCE IN THE
*.... ORIGINAL.  THE VALUES FOR EION ARE RETURNED FROM PFSAHA IN
*.... THE VARIABLE PFP(7+).
*
                  do ion = 1, 7
                     eion(ion) = pfp(ion + 6)
                  end do
*
                  t(j)    = t(j)    * factrm 
                  tk(j)   = tk(j)   * factrm
                  tkev(j) = tkev(j) * factrm
                  call pfsaha((j), (id), (ncomp), c5, pfm)
*
                  t(j)    = t(j)    * factre 
                  tk(j)   = tk(j)   * factre
                  tkev(j) = tkev(j) * factre 
                  ion = ncomp
                  edens(j) = edens(j) + xnmol(j, jmol) * tk(j) *
     &               (eion(ion) / tkev(j) + (pfp(ion) - pfm(ion)) / 
     &               (pfp(ion) + pfp(ion)) * 1000.0d0)
               end do
*
            else
*
               do j = 1, nrhox
                  tplus     = t(j) * factrp
                  tminus    = t(j) * 0.999d0
                  pfplus(j) = exp(-equil(2,jmol) + (equil(3,jmol) + 
     &               (-equil(4,jmol) + (equil(5,jmol) - equil(6,jmol) *
     &               tplus) * tplus) * tplus) * tplus) + 1.0d-30
                  pfmin(j)  = exp(-equil(2,jmol) + (equil(3,jmol) + 
     &               (-equil(4,jmol) + (equil(5,jmol) - equil(6,jmol) *
     &               tminus) * tminus) * tminus) * tminus) + 1.0d-30
               end do
*
               locj1  = locj(jmol)
               locj2  = locj(jmol+1) - 1
               lock   = locj1 - 1
               klimit = .false.
*
               do while(lock .lt. locj2 .and. .not. klimit)
                  lock = lock + 1
                  k    = kcomps(lock)
*
                  if(k .gt. nequa) then
                     klimit = .true.
*
                  else
                     id = idequa(k)
*
                     do j = 1, nrhox
                        t(j)    = t(j)    * factrp
                        tk(j)   = tk(j)   * factrp
                        tkev(j) = tkev(j) * factrp
                        call pfsaha((j), (id), 1, c3, frac)
                        pfplus(j) = pfplus(j) * frac(j, 1)
*
                        t(j)    = t(j)    * factrm 
                        tk(j)   = tk(j)   * factrm 
                        tkev(j) = tkev(j) * factrm 
                        call pfsaha((j),(id), 1, c3, frac)
                        pfmin(j) = pfmin(j) * frac(j, 1)
*
                        t(j)    = t(j)    * factre
                        tk(j)   = tk(j)   * factre 
                        tkev(j) = tkev(j) * factre 
                     end do
*
                  end if
*
               end do
*
               if (.not. klimit) then
*
                  do j = 1, nrhox
                     edens(j) = edens(j) + xnmol(j, jmol) * tk(j) *
     &                  (-equil(1, jmol) / tkev(j) + (pfplus(j) - 
     &                   pfmin(j)) / (pfplus(j) + pfmin(j)) * 1000.0d0)
                  end do
*
               end if
*
            end if
*
         end do
*
         do j = 1, nrhox
            edens(j) = edens(j) / rho(j)
         end do
! --- if ifedens was ture finished*
      else
*
         do k = 1, nequa
*
            do j = 1, nrhox
               xnsave(j, k) = xnz(j, k)
            end do
*
         end do
*
         write (16, 112) 'rhox', 't', 'p', 'xne', 'xnatom', 'rho',
     &       (j, rhox(j), t(j), p(j), xne(j), xnatom(j), rho(j), 
     &        j = 1, nrhox )
112      format ( /////
     &           11x, a, 9x, a, 11x, a, 10x, a, 8x, a, 8x, a /
     &           (i5, 1p6e12.3) )
         jmol1 = 1
*
         do while (jmol1 .le. nummol)
            write (16, 113) 'molecular number densities'
113         format ( /// 35x, a /)
            jmol10 = jmol1
            id     = code(jmol1)
*
            do while (int(code(jmol10 + 1)) .eq. id)
               jmol10 = jmol10 + 1
            end do
*
 
            write (16, 114) (code(jmol), jmol = jmol1, jmol10)
114         format (5x, 6f12.2)
*
            do j = 1, nrhox
               write(16, 115) j, 
     &                       (xnmol(j, jmol), jmol = jmol1, jmol10)
115            format ( i5, 1p6e12.3)
            end do
*
            jmol1 = jmol10 + 1
         end do
*
!------ if mode .ne. 2 
         if (mode .ne. 2 .and. mode .ne. 12) then
*
            do k = 2, nequa
               id = idequa(k)
*
               if (id .eq. 100) then
*
                  do j = 1, nrhox
                     xnz(j, k) = xnz(j, k) * 0.5d0 / 2.4148d15 / t(j) /
     &                           sqrt(t(j) )
                  end do
*
               else
*
*.... CALCULATE PARTITION FUNCTIONS
*
                  do j = 1, nrhox
                     call pfsaha ((j), (id), 1, c3, frac)
                     xnz(j, k) = xnz(j, k) / frac(j, 1) / 1.8786d20 /
     &                           sqrt((atmass(id) * t(j)) ** 3)
                  end do
*
               end if
*
            end do
*
            do jmol = 1, nummol
               ncomp = locj(jmol + 1) - locj(jmol)
*
               if (equil(1, jmol) .eq. 0.0d0) then
!--- equil(1, ..) is zero only for atoms 
                   id = code(jmol)
*
                  do j = 1, nrhox
                     call pfsaha ((j), (id), (ncomp), c3, frac)
                     xnmol(j, jmol) = xnmol(j, jmol) / frac(j, 1)
                  end do
!--- else is for molecules
               else
*
                  do j = 1, nrhox
                     xnmol(j, jmol) = exp(equil(1, jmol) / tkev(j))
                  end do
*
                  amass = 0.0d0
                  locj1 = locj(jmol)
                  locj2 = locj(jmol + 1) - 1
*
                  do lock = locj1, locj2
                     k = kcomps(lock)
 
                     if (k .eq. nequa1) then
*
                        do j = 1, nrhox
                           xnmol(j, jmol) = xnmol(j, jmol) /
     &                                      xnz(j,nequa)
                        end do
*
                     else
                        id = idequa(k)
                        if (id .lt. 100) amass = amass + atmass(id)
                        do j = 1, nrhox
                           xnmol(j, jmol) = xnmol(j, jmol) * xnz(j, k)
                        end do
*
                     end if
*
                  end do
*
                  do j = 1, nrhox
                     xnmol(j, jmol) = xnmol(j, jmol) * 1.8786d20 *
     &                                sqrt((amass * t(j)) ** 3)
                  end do
*
               end if
*
            end do
*
         end if
!------ end if for mode .ne. 2  ------*
         if (ifpnch(iterat) .eq. 5) then
            write (16, 150) 'number densities / partition functions'
150         format ( / 20x, a)
            write (16, 151) nummol, ' molecules'
            write (17, 151) nummol, ' molecules'
151         format (i5, a )
*
            do jmol = 1, nummol
               write (16, 152) code(jmol), (xnmol(j, jmol), j=1, nrhox)
               write (17, 152) code(jmol), (xnmol(j, jmol), j=1, nrhox)
152            format (f20.2 / (1p8e10.3))
            end do
*
            write (16, 158) ' xnatom', (xnatom(j), j = 1, nrhox)
            write (16, 158) ' rho',    (rho(j),    j = 1, nrhox)
            write (16, 158) ' xne',    (xne(j),    j = 1, nrhox)
*
            write (17, 158) ' xnatom', (xnatom(j), j = 1, nrhox)
            write (17, 158) ' rho',    (rho(j),    j = 1, nrhox)
            write (17, 158) ' xne',    (xne(j),    j = 1, nrhox)
158         format (a / (1p8e10.3))
         end if
*
      end if
*
      end subroutine 
*
*********** E N D   O F   S U B R O U T I N E   N M O L E C ************
*
      subroutine pfsaha (j, iz, nion, mode, answer)


      implicit none
*
*.... 1988 MAR - TO AVOID OVER/UNDERFLOW AT LOW DENSITY
*.... 1993 JAN - TO USE THE SAME CUTOFF AS BOB TO GET THE SAME RESULT
*....    MODE 1 RETURNS IONIZATION FRACTION /PARTITION FUNCTION
*....    MODE 2 RETURNS IONIZATION FRACTION
*....    MODE 3 RETURNS PARTITION FUNCTION
*....    MODE 4 RETURNS NUMBER OF ELECTRONS PRODUCED
*....    MODE 5 RETURNS ANSWER(ION) = PF   ANSWER(ION+7) = IP
*....    MODE + 10 RETURN ALL IONS TO NION.   MODE ALONE RETURN NION ONLY.
*
*----------------------------- COMMONS ---------------------------------
*
!      include 'common.sizebl'
      include 'common.depart'
      include 'common.stateb'
      include 'common.tempbl'
*
*--------------------------- DUMMY ARGUMENTS ---------------------------
*
      double precision, intent(out) ::  answer (maxd, *)
      integer, intent(in) ::  j, iz, nion, mode
*
*--------------------------- LOCAL VARIABLES --------------------------
*
      double precision  cf, charge, d1, d2, debye, dt, excess, f(6), g,
     &                  ip(6), p1, p2, part(6), pmin, potlo(6), potlow,
     &                  scale(4), t2000, tv, tvi, z
      integer  i, ion, it, k1, k2, k3, kp1, kscale, l, locz(29),
     &         mode1, n, nion2, nions, nnn(6, 365), nnn100
      logical  skip

! energy levels for more accurate partitions function for He, C, O, ...
      real(kind=8)  EHYD(6),EHE1(29),EHE2(6),EC1(14),EMG1(11),EAL1(9),
     & ESI1(11),EC2(6),EMG2(6),EAL2(1),ESI2(6),ECA1(8),ECA2(5),ENA1(8),
     & EO1(13),EB1(7),EK1(8)

      real(kind=8) GHYD(6),GHE1(29),GHE2(6),GC1(14),GMG1(11),GAL1(9),
     &GSI1(11),GC2(6),GMG2(6),GAL2(1),GSI2(6),GCA1(8),GCA2(5),GNA1(8),
     &GO1(13),GB1(7),GK1(8)

      real(kind=8) b

      save     locz, nnn, scale
!*
!*---------------------------- INITIALIZATION --------------------------
!*
      data locz /  1,   3,   6,  10,  14, 18,  22,  27,  33,  39,  45,
     &            51,  57,  63,  69,  75,  81, 86,  91,  96, 101, 106,
     &           111, 116, 121, 126, 131, 136, 141 /
!*
!*....  ( 1)( 2)   ( 3)( 4)   ( 5)( 6)   ( 7)( 8)   ( 9)(10)   ( IP ) G  REF
      data ((nnn(i, n), i = 1, 6), n = 1, 17) /
     & 200020001, 200020011, 201620881, 231228281, 378953411,  1359502, d+f 1.00
     & 100010001, 100010001, 100010001, 100010001, 100010001,  1359500, g   1.01
     & 100010001, 100010011, 102111241, 145022061, 363059451,  2458104, d+f 2.00
     & 200020001, 200020071, 208524971, 382669341, 128222452,  5440302, d+f 2.01
     & 100010001, 100010001, 100010001, 100010001, 100010001,  5440300, g   2.02
     & 200020011, 201220481, 212922881, 258731081, 394251691,   538901, d+f 3.00
     & 100010001, 100010201, 126225521,  67216512, 351165562,  7561907, d+f 3.01
     & 200020001, 200020211, 227936571,  69610342, 137217102, 12241800, d+f 3.02
     & 100010001, 100010001, 100010001, 100010001, 100010001, 12241800, g   3.03
     & 100010051, 104311441, 131615641, 190623681, 298037691,   931900, ael 4.00
     & 200120231, 211422771, 249627631, 309034911, 398545051,  1820600, ael 4.01
     & 100010001, 100010201, 126225521,  67216512, 351165562, 15385000, ael 4.02
     & 200020001, 200020011, 201220661, 223426161, 332644691, 21765700, ael 4.03
     & 600060001, 600560281, 608761991, 637466191, 693973361,   829500, ael 5.00
     & 100310831, 132016901, 214226411, 315736741, 419147071,  2514900, ael 5.01
     & 200721061, 233526401, 297533311, 369040481, 440747651,  3792000, ael 5.02
     & 100010001, 100010001, 100010001, 100010001, 100010001, 25929800/ g   5.03
!*
      data ((nnn(i,n), i = 1, 6), n =  18,  32) /
     & 893292271,  96110042, 105311262, 126315202, 196126432,  1125508, d+f 6.00
     & 595060251, 620865751, 713280191,  95712292, 167623542,  2437501, d+f 6.01
     & 105513201, 180324851, 341851341,  88416332, 296550722,  4787101, d+f 6.02
     & 204922771, 262630421, 350941931, 494556971, 644872001,  6447600, d+f 6.03
     & 403141851, 457051681, 594071181,  92913362, 203331152,  1452915, d+f 7.00
     & 919899541, 107211512, 124914302, 182526232, 403762662,  2959202, d+f 7.01
     & 596862721, 684177081,  88110342, 128317062, 239334312,  4742501, d+f 7.02
     & 112816481, 240733751, 462068491, 116419932, 283736822,  7744900, d+f 7.03
     & 210124681, 293634211, 391145791, 539862151, 703178471,  9786200, d+f 7.04
     & 874789691, 924795711,  99410492, 115213492, 169022242,  1361307, d+f 8.00
     & 424151091, 622874781,  91312832, 221842502,  79914013,  3510711, d+f 8.01
     &  95610702, 118113032, 149619922, 329761642, 101914173,  5488500, d+f 8.02
     & 603567171, 775391141, 106612482, 143716252, 181420032,  7739300, d+f 8.03
     & 124420321, 306943181, 606281181, 101712232, 142916342, 11387300, d+f 8.04
     & 215026541, 323137551, 421546491, 508255151, 594863811, 13807900/ ael 8.05
*
      data ((nnn(i,n), i = 1, 6), n =  33,  44) /
     & 575958511, 589859231, 595860671, 636470031, 815199581,  1741802, d+f 9.00
     & 900296401, 102610802, 113912542, 152921152, 318348952,  3498003, d+f 9.01
     & 469162651, 791295541, 121419552, 402686872, 154822203,  6264500, d+f 9.02
     &  99511422, 129214572, 170523002, 320140922, 498458762,  8713900, d+f 9.03
     & 615472711,  87710602, 127215002, 172919582, 218624152, 11421300, d+f 9.04
     & 135324181, 377252001, 661580261,  94410852, 122613672, 15711700, ael 9.05
     & 100010001, 100010051, 105313051, 210239461,  74013022,  2155808, d+f10.00
     & 580158751, 591759741, 642687101, 159332652,  64111533,  4106907, d+f10.01
     &  93510272, 110411662, 127116062, 257647882,  75110223,  6350000, d+f10.02
     & 529774371,  94611322, 135816202, 188221442, 240626682,  9701900, d+f10.03
     & 103312152, 140616092, 181320182, 222224262, 263128352, 12630000, ael10.04
     & 629178711,  98311802, 136715512, 173619202, 210422892, 15790900/ ael10.05
*
      data ((nnn(i,n), i = 1, 6), n =  45,  62) /
     & 200020001, 200320211, 207322131, 253031421, 417657451,   513802, d+f11.00
     & 100010001, 100010161, 119621261,  50711872, 246445382,  4728901, d+f11.01
     & 580158751, 591860351,  71813142, 321968812, 106014333,  7165000, d+f11.02
     &  96910772, 116012242, 130714232, 153916552, 177118872,  9888000, d+f11.03
     & 601386081, 108812932, 148916832, 187820722, 226624612, 13836900, ael11.04
     & 105712442, 144616652, 189221182, 234425702, 279630222, 17209000, ael11.05
     & 100010011, 101410621, 118414581, 204831781, 509479731,   764404, d+f12.00
     & 200120051, 202921001, 226926901, 368457091,  92814872,  1503101, d+f12.01
     & 100010001, 100110611, 177455431, 176546012,  99718753,  8011905, d+f12.02
     & 579758751, 591459501, 600560591, 611461681, 622362781, 10928900, ael12.03
     & 100611232, 120612752, 134214102, 147815462, 161416822, 14122900, ael12.04
     & 674896701, 121814462, 167018942, 211723412, 256527892, 18648900, ael12.05
     & 558857701, 583558761, 593260591, 635969541, 796790971,   598400, d+f13.00
     & 100310211, 110313021, 172828201,  55311252, 215637942,  1882203, d+f13.01
     & 200320201, 208622331, 250530971, 410251081, 611571211,  2844000, d+f13.02
     & 100010001, 100210881, 207436531, 523168101, 838999681, 11996000, d+f13.03
     & 577758651, 591259631, 604461351, 622563161, 640764981, 15377000, ael13.04
     & 103511582, 124713242, 140014772, 155316292, 170517812, 19042000/ ael13.05
*
      data ((nnn(i,n), i = 1, 6), n =  63,  74) /
     & 825189211,  95210052, 106211532, 134317202, 237934082,   814913, d+f14.00
     & 563057761, 588160311, 631768671, 791097651, 127817282,  1634000, d+f14.01
     & 101110771, 126716471, 232438081,  71914052, 262045302,  3346001, d+f14.02
     & 200720521, 217224081, 284439171, 551370951,  86810262,  4513000, d+f14.03
     & 100010001, 100210881, 207436531, 523168101, 838999681, 16672900, fak14.04
     & 575458521, 591459851, 610063201, 672674071, 843698661, 20510900, ael14.05
     & 402643441, 496757481, 658274401, 833492941, 103511532,  1048300, ael15.00
     & 874497931, 106011282, 119812802, 138415142, 164717802,  1972000, ael15.01
     & 564058061, 604164611, 709579551,  90410172, 112912422,  3015500, ael15.02
     & 100811411, 149720221, 280936121, 441552181, 602168241,  5135400, ael15.03
     & 200420781, 227025361, 281430911, 336936471, 392542021,  6500700, ael15.04
     & 100010001, 100010001, 100010001, 100010001, 100010001, 22041300/ g  15.05
*
      data ((nnn(i,n), i = 1, 6), n =  75,  90) /
     & 822887891, 930697831, 102610932, 121614492, 185124742,  1035708, d+f16.00
     & 443056011, 694982961,  96911522, 144218572, 227326892,  2339900, d+f16.01
     &  91610392, 113512242, 136416942, 233429882, 364242962,  3500000, d+f16.02
     & 560058861, 633871081,  82410062, 123314602, 168619132,  4728900, d+f16.03
     & 104512901, 177025421, 375163021, 122420462, 286036742,  7250000, d+f16.04
     & 202321571, 241428261, 358355061,  78310152, 124814802,  8802800, d+f16.05
     & 538155931, 571657911, 598067191,  89013782, 227737172,  1300916, d+f17.00
     & 873396771, 104411072, 118513532, 175525872, 406763932,  2379903, d+f17.01
     & 506569571,  87610522, 134421682, 439092662, 182132573,  3990006, d+f17.02
     &  95110872, 120013232, 154921252, 345149322, 641378942,  5350000, d+f17.03
     & 558960371, 677779341,  95311692, 138816082, 182720472,  6780000, d+f17.04
     & 100010001, 100010051, 106913911, 240147261,  90716112,  1575411, d+f18.00
     & 550256831, 578158781, 636585461, 151530162,  58010303,  2762007, d+f18.01
     &  92110362, 112412002, 133216772, 254443722,  76512833,  4090003, d+f18.02
     & 582082081, 103112292, 149920212, 309750502, 720793642,  5978900, d+f18.03
     &  97111072, 123213982, 172625622, 463976582, 106413633,  7500000/ d+f18.04
*
      data ((nnn(i,n), i = 1, 6), n =  91, 105) /
     & 200020011, 200720361, 211923291, 280137141, 525575741,   433803, d+f19.00
     & 100010001, 100110341, 135929551,  79119282, 405274892,  3180905, d+f19.01
     & 554657081, 581260301,  73012702, 285363872, 129023363,  4600005, d+f19.02
     &  96010862, 118413212, 180836632,  90321023, 416863253,  6090000, d+f19.03
     & 657793361, 119515082, 195826322, 352944302, 533162332,  8259900, d+f19.04
     & 100110061, 104311741, 145919971, 294345051,  69010322,   611003, d+f20.00
     & 205822781, 279234761, 427553061, 688994901, 136319772,  1186701, d+f20.01
     & 100010001, 100510821, 168744821, 130232522,  69012813,  5121003, d+f20.02
     & 555157161, 585662471,  82816862,  42510013, 168423663,  6700000, d+f20.03
     &  99411262, 123814062, 182930402, 484766392,  84310223,  8438900, d+f20.04
     & 924696691, 105212282, 151219062, 240530032, 368944512,   653900, ael21.00
     & 190424662, 297634542, 391743752, 482952832, 573761912,  1280000, ael21.01
     & 976799291, 101110322, 105810882, 111911502, 118112122,  2475000, ael21.02
     & 100010001, 100510821, 168744821, 130232522,  69012813,  7390000, fak21.03
     & 555157161, 585662471,  82816862,  42510013, 168423663,  9200000/ fak21.04
*
      data ((nnn(i,n), i = 1, 6), n = 106, 120) /
     & 181021172, 260333222, 430155582, 710089242, 110213293,   681900, d+f22.00
     & 474659872, 721284672,  98211413, 134515623, 177919963,  1356900, d+f22.01
     & 228327012, 308134272, 381143862, 534563472, 734983512,  2747000, d+f22.02
     & 971498311,  99210032, 102610572, 108711172, 114711782,  4324000, d+f22.03
     & 100010001, 100510821, 168744821, 130232522,  69012813,  9980000, fak22.04
     & 272835172, 425851532, 632278322,  97212013, 146817723,   674000, ael23.00
     & 373954132, 743597002, 121414713, 173920143, 229225713,  1464900, ael23.01
     & 323142642, 519660272, 679975352, 824789522,  96610363,  2930900, ael23.02
     & 248329302, 324234952, 373439752, 421744582, 469949412,  4800000, ael23.03
     & 970698231, 990699881, 100710152, 102410322, 104010482,  6500000, ael23.04
     & 717277611,  92911652, 152620872, 295141952, 550468122,   676400, d+f24.00
     &  71611552, 205635512, 558281952, 115315823, 205625293,  1649000, d+f24.01
     & 280639822, 538369722,  87610823, 129115003, 170919183,  3095000, d+f24.02
     & 377150952, 616070292, 791788382,  97610683, 116012523,  5000000, d+f24.03
     & 264730962, 341436462, 394042872, 463549832, 533056782,  7300000/ d+f24.04
*
      data ((nnn(i,n), i = 1, 6), n = 121, 135) /
     & 600060321, 629270891,  86911302, 151020222, 267534752,   743100, ael25.00
     & 739594821, 139921212, 309342852, 567372412,  97112553,  1563600, ael25.01
     &  98417472, 265535782, 454754842, 641973532, 828792212,  3369000, ael25.02
     & 328847052, 586668342, 771785912,  94710343, 112112093,  5300000, ael25.03
     & 422055132, 636770792, 779285062, 921999322, 106411363,  7600000, ael25.04
     & 197023222, 274433302, 416753952, 723799822, 139419053,   787038, d+f26.00
     & 409453722, 686687452, 110213823, 174322233, 286437043,  1617902, d+f26.01
     & 262136422, 501167232,  87911303, 138916483, 190721673,  3064300, d+f26.02
     &  98723522, 420363072,  87011423, 145117913, 215925463,  5700000, ael26.03
     & 388854482, 666275742, 846693572, 102511143, 120312923,  7900000, d+f26.04
     & 199427202, 335740022, 474957182, 708090462, 118315403,   786000, d+f27.00
     & 279739202, 490858232, 684582472, 104713233, 159818733,  1704900, d+f27.01
     & 279836622, 461857562, 720693022, 124915873, 192522633,  3349000, d+f27.02
     & 262136422, 501167232,  87911303, 138916483, 190821673,  5300000, fak27.03
     &  98723522, 420363072,  87011423, 145117913, 215925463,  8300000/ fak27.04
*
      data ((nnn(i,n), i = 1, 6), n = 136, 149) /
     & 227027622, 306233052, 356839222, 446052912, 652382292,   763314, d+f28.00
     & 108416342, 222428472, 353944332, 577378932, 110314303,  1814900, d+f28.01
     & 198724282, 293236452, 468362702,  86511123, 136016073,  3516000, d+f28.02
     & 279836622, 461857562, 720693022, 124915873, 192522633,  5600000, fak28.03
     & 262136422, 501167232,  87911303, 138916483, 190721673,  7900000, fak28.04
     & 201620781, 231026761, 314737361, 450555381, 692386911,   772301, d+f29.00
     & 109415761, 247938311,  58910042, 190937022,  68311693,  2028903, d+f29.01
     & 897195961, 107212972, 165021182, 260230862, 356940532,  3682900, d+f29.02
     & 100010001, 100410231, 108712611, 167124841, 388460411,   939102, d+f30.00
     & 200020021, 201620761, 223726341, 351352061,  80812472,  1796001, d+f30.01
     & 100610471, 122617301, 300566361, 149924112, 332342352,  3970000, d+f30.02
     & 403245601, 493151431, 529654331, 559358091, 611065171,   600000, ael31.00
     &  99710051, 104511541, 135016501, 208226431, 321837921,  2050900, ael31.01
     & 199820071, 204521391, 229124761, 266028451, 302932131,  3070000/ ael31.02
*
      data ((nnn(i,n), i = 1, 6), n = 150, 164) /
     & 502665261, 755183501, 901496201, 102410942, 117912812,   787900, ael32.00
     & 422848161, 512153401, 557458941, 636270361, 794489061,  1593000, ael32.01
     & 100010261, 114613921, 175221251, 249828711, 324436181,  3421000, ael32.02
     & 403143241, 491856701, 649173781, 840396751, 113013392,   981000, ael33.00
     & 593676641, 884697521, 105911572, 129515012, 180322212,  1858700, ael33.01
     & 484470541,  91510972, 125614082, 157017612, 199722912,  2829900, ael33.02
     & 630172361, 799686381, 919797221, 102810942, 117712832,   975000, ael34.00
     & 438055511, 691582151,  94510732, 121413672, 152016732,  2150000, ael34.01
     & 651982921,  94610382, 113212492, 139515462, 169718482,  3200000, ael34.02
     & 437347431, 498951671, 538559501,  74710812, 169126672,  1183910, d+f35.00
     & 705183611,  93510092, 111614162, 222932532, 427652992,  2160000, d+f35.01
     & 510869921,  87410312, 123116552, 236530712, 377744832,  3590000, d+f35.02
     & 100010001, 100010051, 105012781, 198535971,  65911422,  1399507, d+f36.00
     & 461049811, 522254261, 609088131, 168935052,  68612253,  2455908, d+f36.01
     & 759990901, 101911142, 129017782, 302856642,  99414333,  3690000/ d+f36.02
*
      data ((nnn(i,n), i = 1, 6), n = 165, 179) /
     & 200020011, 200720361, 211523021, 269434141, 459163351,   417502, d+f37.00
     & 100010001, 100110321, 129524961,  61014202, 291753192,  2750004, d+f37.01
     & 473650891, 533156051,  66810932, 232950852,  99915303,  4000000, d+f37.02
     & 100110041, 104111741, 146019721, 281941411, 607785251,   569202, d+f38.00
     & 202621931, 255331271, 384347931, 624085761, 122417632,  1102600, d+f38.01
     & 100010001, 100110321, 129524961,  61014202, 291753192,  4300000, fak38.02
     & 791587851, 100012192, 155119942, 254031782, 389946932,   637900, ael39.00
     & 118217102, 220827002, 319036792, 416646512, 513256072,  1223000, ael39.01
     &  92510012, 104710862, 112311612, 120212472, 132814282,  2050000, ael39.02
     & 141320802, 291439702, 531170262,  92712273, 162521053,   684000, d+f40.00
     & 354454352, 724689652, 107212643, 148517093, 193321573,  1312900, d+f40.01
     & 209727032, 324537052, 415446282, 510255752, 604965222,  2298000, d+f40.02
     & 256636022, 465759302, 749693962, 116514243, 171520333,   687900, ael41.00
     & 335157222,  84511463, 147718363, 221826083, 299933893,  1431900, ael41.01
     & 223725352, 280830972, 340937362, 406844002, 473150632,  2503900/ ael41.02
*
      data ((nnn(i,n), i = 1, 6), n = 180, 194) /
     & 703972941,  82610822, 154822682, 327244912, 571469372,   709900, d+f42.00
     &  75714552, 274347322, 718897632, 123414913, 174920063,  1614900, d+f42.01
     & 267645462, 669890262, 115514323, 173620673, 242528083,  2714900, ael42.02
     &  90613732, 184823562, 291735332, 419949102, 565764332,   728000, ael43.00
     & 131318312, 227126932, 311735452, 397644072, 483852692,  1525900, ael43.01
     & 204721673, 234725733, 284031463, 348738613, 426546943,  3000000, ael43.02
     & 176824122, 318941082, 515263202, 761790472, 106112303,   736400, ael44.00
     & 221934642, 501968372,  88911173, 136316243, 189221613,  1675900, ael44.01
     & 210622722, 241025422, 267928262, 297731272, 327834282,  2846000, ael44.02
     & 148520202, 255230902, 364942462, 489656082, 638872352,   746000, ael45.00
     & 153421292, 288137912, 484660322, 720187062, 101011483,  1807000, ael45.01
     & 254537212, 492362292, 770592182, 107312243, 137615273,  3104900, ael45.02
     & 115919651, 320746011, 607576761,  95011642, 141817172,   832900, ael46.00
     & 755087211, 105913442, 173122222, 282034722, 412247732,  1941900, ael46.01
     & 180223462, 289735212, 414247632, 538460052, 662672472,  3292000/ ael46.02
*
      data ((nnn(i,n), i = 1, 6), n = 195, 209) /
     & 200020001, 200220141, 206422141, 257633021, 455164681,   757403, d+f47.00
     & 100810581, 125817401, 260641031,  66210072, 135316982,  2148000, d+f47.01
     & 795887491,  97711762, 156620252, 248329422, 340038582,  3481900, d+f47.02
     & 100010001, 100410241, 109212891, 176827421, 444268771,   899003, d+f48.00
     & 200020021, 201720921, 233329881, 451475371, 127520782,  1690301, d+f48.01
     & 100310281, 114815371, 246138311, 519265531, 791492761,  3747000, d+f48.02
     & 252431921, 368440461, 433746521, 512259221, 723389021,   578400, d+f49.00
     & 100110071, 104611651, 146118581, 225426511, 304734431,  1886000, d+f49.01
     & 200120111, 205021611, 243628031, 317035371, 390442701,  2802900, d+f49.02
     & 232637101, 488058571, 669074381, 816189091,  97210632,   734200, ael50.00
     & 286335941, 408144471, 479351961, 571862901, 686274341,  1462700, ael50.01
     & 100010251, 114013811, 175321601, 256829751, 338337901,  3049000, ael50.02
     & 404043481, 494656811, 646772781, 813490751, 101411372,   863900, ael51.00
     & 303147981, 618472951, 827392621, 103711702, 131214532,  1650000, ael51.01
     & 313037601, 429347901, 536260591, 689477591, 862494881,  2529900/ ael51.02
*
      data ((nnn(i,n), i = 1, 6), n = 210, 224) /
     & 526258801, 657372351, 784284071, 897095741, 102711082,   900900, ael52.00
     & 440855541, 686481251,  93810792, 125414792, 176321132,  1860000, ael52.01
     & 349054751, 699883081,  96611302, 134216202, 197724212,  2800000, ael52.02
     & 405342041, 438645621, 475751071, 587974491, 102214572,  1045404, d+f53.00
     & 568567471, 773485861,  94510362, 112712182, 130914002,  1909000, d+f53.01
     & 514269581,  86910562, 130716652, 215327742, 351843662,  3200000, ael53.02
     & 100010001, 100010091, 109515351, 291060661, 119621482,  1212716, d+f54.00
     & 414844131, 465649111, 538464651,  87112232, 158019362,  2120000, d+f54.01
     & 615475101, 867797531, 112213462, 157618062, 203622662,  3209900, d+f54.02
     & 200020001, 201020501, 215623871, 283536181, 462756261,   389300, d+f55.00
     & 100010001, 100310371, 119016501, 269146361,  77912412,  2510000, d+f55.01
     & 424445601, 481750061, 516953311, 549356551, 581759791,  3500000, d+f55.02
     & 101210791, 135119351, 282340571, 574580391, 111015062,   521002, d+f56.00
     & 262638611, 504160621, 698579371,  91010692, 129115952,  1000000, d+f56.01
     & 100010001, 100310351, 118416321, 264945521,  76512182,  3700000/ fak56.02
*
      data ((nnn(i,n), i = 1, 6), n = 225, 239) /
     &  71111992, 172323592, 312540402, 510763182, 765791012,   558000, ael57.00
     & 204529582, 383647882, 582469262, 807992692, 104911723,  1106000, ael57.01
     &  94712552, 148416582, 179819212, 203621522, 227424042,  1916900, ael57.02
     & 295959132, 103515693, 215527593, 335939413, 449650223,   565000, ael58.00
     &  79718153, 289639443, 495159253, 686877533, 863794813,  1085000, ael58.01
     & 298640242, 475053692, 596965912, 725379692, 872094692,  2008000, ael58.02
     & 460693672, 158523823, 327242303, 519661563, 709379783,   541900, fak59.00
     & 455480232, 114014653, 178521013, 240927073, 299232633,  1055000, ael59.01
     &  46410533, 183826893, 354443773, 518459633, 674375243,  2320000, ael59.02
     & 139623042, 364860002,  96114603, 209828633, 373446973,   549000, ael60.00
     & 460493692, 158523823, 327142303, 519661563, 709279783,  1073000, ael60.01
     & 455480232, 114014653, 178521013, 240927073, 299232633,  2000000, fak60.02
     & 131720482, 280535692, 441254492, 676583972, 103412583,   555000, ael61.00
     & 139623042, 364860002,  96114603, 209828633, 373446973,  1089900, fak61.01
     & 460493682, 158523823, 327142303, 519661563, 709279783,  2000000/ fak61.02
*
      data ((nnn(i,n), i = 1, 6), n = 240, 254) /
     &  92915672, 222431062, 444763802,  89612173, 159520253,   562900, ael62.00
     & 315059662,  97114563, 204627093, 342541693, 490556383,  1106900, ael62.01
     & 269037812, 520270372,  91111273, 133915483, 172719093,  2000000, ael62.02
     & 800080571, 851699301, 127617362, 240433032, 444958442,   568000, ael63.00
     & 125416052, 211828182, 375549622, 644381732, 101112213,  1125000, ael63.01
     & 800080571, 851699301, 127617362, 240433032, 444958442,  2000000, fak63.02
     & 240432982, 427555202, 708489962, 112613853, 167319843,   615900, ael64.00
     & 534793262, 139219123, 247730843, 371043333, 495055893,  1210000, ael64.01
     & 364145232, 514756362, 604864112, 673870372, 732276072,  2000000, ael64.02
     & 480767202,  89011393, 144118243, 230028753, 354142883,   584900, ael65.00
     & 480767192,  89011393, 144118243, 230028753, 354142883,  1151900, fak65.01
     & 480767202,  89011393, 144118243, 230028753, 354142883,  2000000, fak65.02
     & 343147532, 645887152, 115314793, 183322063, 257729373,   593000, fak66.00
     & 343147532, 645887142, 115314793, 183322063, 257729373,  1167000, ael66.01
     & 343147532, 645887142, 115314793, 183322063, 257729373,  2000000/ fak66.02
*
      data ((nnn(i,n), i = 1, 6), n = 255, 269) /
     & 222635002, 542276772, 100312353, 145716713, 187020703,   602000, fak67.00
     & 222635002, 542276772, 100312353, 145716713, 187020703,  1180000, fak67.01
     & 222635002, 542276772, 100312353, 145716713, 187020703,  2000000, ael67.02
     & 133715382, 209130152, 429859382,  79410293, 129815983,   609900, ael68.00
     & 265934782, 497877532, 120517733, 245032063, 400448073,  1193000, ael68.01
     & 265934782, 497877532, 120517733, 245032063, 400448073,  2000000, fak68.02
     & 800381111,  87510702, 147621462, 310343462, 585475982,   618000, ael69.00
     & 156718872, 279244452, 678196342, 128316243, 197823443,  1205000, ael69.01
     &  93517192, 364666132, 103414613, 192624193, 293334613,  2370000, ael69.02
     & 100010011, 101310651, 118613951, 169120661, 250629971,   625000, ael70.00
     & 200120901, 270345231,  81714042, 223533112, 461959862,  1217000, ael70.01
     & 100312561, 250851931,  91914182, 198626022, 323638692,  2000000, ael70.02
     & 514664441, 759086851,  99211442, 133315612, 182721252,   609900, ael71.00
     & 125924831, 438667801,  98714112, 199727872, 380850742,  1389900, ael71.01
     & 323948621, 661297271, 158626482, 426865032,  93712843,  1900000/ ael71.02
*
      data ((nnn(i,n), i = 1, 6), n = 270, 284) /
     & 659294081, 128016962, 222528952, 372047062, 585171462,   700000, ael72.00
     &  99117882, 274638812, 520867322,  84410313, 123314453,  1489900, ael72.01
     & 187427702, 343739872, 448049452, 539358282, 625266642,  2329900, ael72.02
     &  65210892, 171325762, 373552252, 705192012, 116414343,   787900, ael73.00
     & 192837842, 600784802, 111113823, 165419233, 218524383,  1620000, ael73.01
     &  99117872, 274638812, 520867312,  84410313, 123314453,  2400000, fak73.02
     & 398981651, 130019172, 273438022, 516168382,  88411163,   797900, ael74.00
     & 131429482, 523279952, 111414623, 183422233, 262130233,  1770000, ael74.01
     & 192837842, 600784792, 111113823, 165419233, 218524383,  2500000, fak74.02
     & 600963001,  75910412, 150121572, 301940972, 539168952,   787000, ael75.00
     &  73710852, 190731262, 464964142,  83810503, 127315053,  1660000, ael75.01
     & 131429482, 523279952, 111414623, 183422233, 262130233,  2600000, fak75.02
     & 110815502, 216829732, 398752322, 672484682, 104612673,   850000, ael76.00
     & 168225972, 362046562, 566766422, 757484612,  93010103,  1700000, ael76.01
     &  73710852, 190731262, 464964142,  83810503, 127315053,  2700000/ fak76.02
*
      data ((nnn(i,n), i = 1, 6), n = 285, 299) /
     & 129117892, 239430882, 388748292, 596173252,  89510843,   910000, ael77.00
     & 110815502, 216829732, 398752322, 672484682, 104612673,  2000000, fak77.01
     & 168225972, 362046562, 566766422, 757484612,  93010103,  2800000, fak77.02
     & 158918512, 207523002, 254328242, 316335762, 407246582,   900000, ael78.00
     &  98115462, 224930742, 401150612, 623475412,  89910583,  1855900, ael78.01
     & 110815502, 216829732, 398752322, 672484682, 104612673,  2900000, fak78.02
     & 203222611, 265731251, 364042301, 494958601, 702084731,   922000, ael79.00
     & 120521331, 357753801,  75310062, 130516572, 206925452,  2050000, ael79.01
     & 651780821, 108814772, 195925252, 316338622, 460853882,  3000000, ael79.02
     & 100010001, 100110111, 105211851, 152122101, 341552811,  1043002, d+f80.00
     & 200320211, 210023021, 268834231, 480472341, 111416912,  1875000, d+f80.01
     & 104012871, 186129471, 458664151,  82410072, 119013732,  3420000, d+f80.02
     & 200420711, 222424271, 265429161, 325637371, 442853911,   610500, ael81.00
     & 100010021, 101910801, 121414641, 189525811, 358949721,  2041900, ael81.01
     & 200020311, 216624611, 296337451, 489064791,  85711212,  2979900/ ael81.02
*
      data ((nnn(i,n), i = 1, 6), n = 300, 314) /
     & 103411711, 147819101, 244331781, 434862751,  93113762,   741404, d+f82.00
     & 204122231, 248227841, 311535621, 429153941, 651976431,  1502800, d+f82.01
     & 100210131, 106812201, 154522671, 381665951,  95512512,  3192900, d+f82.02
     & 400140351, 416944121, 474851591, 564362181, 690477231,   728700, ael83.00
     & 106814451, 204427341, 350744811, 586879131, 108314772,  1667900, ael83.01
     & 205523051, 264830231, 345439921, 469156001, 675281671,  2555900, ael83.02
     & 500950661, 518153561, 559058941, 628968071, 748483501,   843000, ael84.00
     & 443756241, 696282451,  95411012, 128615262, 182922012,  1900000, fak84.01
     & 336953201, 682481011,  93810882, 127915272, 184622442,  2700000, fak84.02
     & 402841621, 431544771, 463148311, 520059491, 734896851,   930000, fak85.00
     & 576168741, 788387631,  96910642, 116012552, 135014462,  2000000, fak85.01
     & 490265341, 812797201, 116614322, 179622692, 285035302,  2900000, fak85.02
     & 100010001, 100010031, 102311051, 133018071, 264539391,  1074500, ael86.00
     & 402841621, 431544771, 463148311, 520059491, 734996851,  2000000, fak86.01
     & 576168741, 788387631,  96910642, 116012552, 135014462,  3000000/ fak86.02
*
      data ((nnn(i,n), i = 1, 6), n = 315, 329) /
     & 200020011, 201220591, 218124481, 296538611, 488859141,   400000, fak87.00
     & 100010001, 100010031, 102311051, 133018071, 264539401,  2200000, fak87.01
     & 421645151, 477449611, 511852711, 542455761, 572958821,  3300000, fak87.02
     & 100010041, 105212131, 153220271, 270435641, 460258111,   527600, ael88.00
     & 201221791, 258131471, 381645781, 546365131, 777592781,  1014400, ael88.01
     & 100010001, 100010031, 102311051, 133018071, 264539391,  3400000, fak88.02
     & 510064491,  82710872, 142718412, 232328712, 348341572,   690000, ael89.00
     & 228951571,  88513232, 183324132, 305537492, 448152402,  1210000, ael89.01
     & 723989131, 103511752, 130814352, 155416652, 177018682,  2000000, ael89.02
     & 620099241, 162725772, 391457072,  80110833, 141818023,   600000, ael90.00
     & 620099241, 162725772, 391457072,  80110833, 141818023,  1200000, fak90.01
     & 620099251, 162725772, 391457072,  80110833, 141818023,  2000000, fak90.02
     & 347877992, 129318323, 240730533, 380546863, 570368573,   600000, ael91.00
     & 347877992, 129318323, 240730533, 380546863, 570368573,  1200000, fak91.01
     & 347777992, 129318323, 240730533, 380546863, 570368573,  2000000/ fak91.02
*
      data ((nnn(i,n), i = 1, 6), n = 330, 344) /
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, ael92.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak92.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, fak92.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak93.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak93.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, fak93.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak94.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak94.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, fak94.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak95.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak95.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, fak95.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak96.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak96.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000/ fak96.02
*
      data ((nnn(i,n), i = 1, 6), n = 345, 353) /
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak97.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak97.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, fak97.02
     & 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak98.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak98.01
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000, fak98.02
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000, fak99.00
     & 209530092, 450866762,  96613623, 186524763, 318839893,  1200000, fak99.01
     & 209530092, 450866762,  96613623, 186524763, 318839893,  2000000/ fak99.02
*
      data ((nnn(i,n), i = 1, 6), n = 354, 365) /
     & 893292271,  96110042, 105311262, 126315202, 196126432,  1125508, d+f 6.00
     & 595060251, 620865751, 713280191,  95712292, 167623542,  2437501, d+f 6.01
     & 105513201, 180324851, 341851341,  88416332, 296550722,  4787101, d+f 6.02
     & 204922771, 262630421, 350941931, 494556971, 644872001,  6447600, d+f 6.03
     & 100010001, 100010001, 100010001, 100010001, 100010001, 39207700, g   6.04
     & 200020001, 200020001, 200020001, 200020001, 200020001, 48998100, g   6.05
     & 403141851, 457051681, 594071181,  92913362, 203331152,  1452915, d+f 7.00
     & 919899541, 107211512, 124914302, 182526232, 403762662,  2959202, d+f 7.01
     & 596862721, 684177081,  88110342, 128317062, 239334312,  4742501, d+f 7.02
     & 112816481, 240733751, 462068491, 116419932, 283736822,  7744900, d+f 7.03
     & 210124681, 293634211, 391145791, 539862151, 703178471,  9786200, d+f 7.04
     & 100010001, 100010001, 100010001, 100010001, 100010001, 55205700/ g   7.05
*
! more accurate partition functions for He, B, C, O, Na, Mg, AL, Si, K, Ca 
      DATA EHE1/0.,159856.069,166277.546,169087.007,171135.000,
     & 183236.892,184864.936,185564.694,186101.654,186105.065,
     & 186209.471,190298.210,190940.331,191217.14,191444.588,
     & 191446.559,191451.80,191452.08,191492.817,
     & 193347.089,193663.627,193800.78,193917.245,193918.391,
     & 193921.31,193921.37,193922.5,193922.5,193942.57/

      DATA GHE1/1.,3.,1.,9.,3.,3.,1.,9.,15.,5.,3.,3.,1.,9.,15.,5.,21.,
     & 7.,3.,3.,1.,9.,15.,5.,21.,7.,27.,9.,3./
      DATA EHE2/0.,329182.321,390142.359,411477.925,421353.135,
     & 426717.413/
      DATA GHE2/2.,8.,18.,32.,50.,72./
!
      DATA EB1/10.17,28810.,40039.65,47856.99,48613.01,54767.74,
     & 55010.08/
      DATA GB1/6.,12.,2.,10.,6.,10.,2./
!
      DATA EC1/29.60,10192.66,21648.02,33735.20,60373.00,61981.82,
     & 64088.85,68856.33,69722.00,70743.95,71374.90,72610.72,73975.91,
     & 75254.93/
      DATA GC1/9.,5.,1.,5.,9.,3.,15.,3.,15.,3.,9.,5.,1.,9./
      DATA EC2/42.48,43035.8,74931.11,96493.74,110652.10,116537.65/
      DATA GC2/6.,12.,10.,2.,6.,2./
!
      DATA EO1/77.975,15867.862,33792.583,73768.200,76794.978,86629.089,
     & 88630.977,95476.728,96225.049,97420.748,97488.476,99094.065,
     & 99681.051/
      DATA GO1/9.,5.,1.,5.,3.,15.,9.,5.,3.,25.,15.,15.,9./
!
      DATA ENA1/0.,16956.172,16973.368,25739.991,29172.889,29172.839,
     & 30266.99,30272.58/
      DATA GNA1/2.,2.,4.,2.,6.,4.,2.,4./
!
      DATA EMG1/0.,21890.854,35051.264,41197.403,43503.333,
     & 46403.065,47847.797,47957.034,49346.729,51872.526,52556.206/
      DATA GMG1/1.,9.,3.,3.,1.,5.,9.,15.,3.,3.,1./
      DATA EMG2/0.,35730.36,69804.95,71490.54,80639.85,92790.51/
      DATA GMG2/2.,6.,2.,10.,6.,2./
!
      DATA EAL1/74.707,25347.756,29097.11,32436.241,32960.363,
     & 37689.413,38932.139,40275.903,41319.377/
      DATA GAL1/6.,2.,12.,10.,6.,2.,10.,6.,14./
!
      DATA ESI1/149.681,6298.850,15394.370,33326.053,39859.920,
     & 40991.884,45303.310,47284.061,47351.554,48161.459,49128.131/
      DATA GSI1/9.,5.,1.,5.,9.,3.,15.,3.,5.,15.,9./
      DATA ESI2/191.55,43002.27,55319.11,65500.73,76665.61,79348.67/
      DATA GSI2/6.,12.,10.,2.,2.,10./
      DATA EK1/0.,12985.170,13042.876,21026.551,21534.680,21536.988,
     & 24701.382,24720.139/
      DATA GK1/2.,2.,4.,2.,6.,4.,2.,4./
      DATA ECA1/0.,15263.089,20356.265,21849.634,23652.304,31539.495,
     & 33317.264,35831.203/
      DATA GCA1/1.,9.,15.,5.,3.,3.,1.,21./
      DATA ECA2/0.,13686.60,25340.10,52166.93,56850.78/
      DATA GCA2/2.,10.,6.,2.,10./


      data scale / 0.001, 0.01, 0.1, 1.0 /
*
*------------------------------- EXECUTION -----------------------------
*
      mode1 = mode
      if(mode1 .gt. 10) mode1 = mode1 - 10
*
*.... LOWERING OF THE IONIZATION POTENTIAL IN VOLTS FOR UNIT ZEFF
*
      charge = 2.0d0 * xne(j)
      excess = charge - p(j) / tk(j)
*
*.... ALLOWANCE FOR DOUBLY IONIZED HELIUM
*
      if(excess .gt. 0.0d0) charge = charge + 2.0 * excess
      debye  = sqrt(tk(j) / 2.8965d-18 / charge)
      potlow = min(1.0d0, 1.44d-7 / debye)
      tv     = tkev(j)
      tvi    = 1.0d0 / tv
*
      if (iz .gt. 28) then
         n     = 3 * iz + 54
         nions = 3
*
      else if (iz .eq. 6) then
         n     = 354
         nions = 6
*
      else if (iz .eq. 7) then
         n     = 360
         nions = 6
      else
         n     = locz(iz)
         nions = locz(iz + 1) - n
      end if
*
      nion2 = min(nion + 2, nions)
      n     = n - 1
*
      do ion = 1, nion2
         z          = ion
         potlo(ion) = potlow * z
         n          = n + 1
         nnn100     = nnn(6, n) / 100
         ip(ion)    = dble(nnn100) * 1.0d-3
         g          = nnn(6, n) - nnn100 * 100
*
         if(n .eq. 1) then
            skip    = .false.
            part(1) = 2.0d0 * bhyd(j, 1)
            part(1) = part(1) +  
     &                 8.0d0 * bhyd(j, 2) * exp( -10.196d0 * tvi) + 
     &                18.0d0 * bhyd(j, 3) * exp( -12.084d0 * tvi) + 
     &                32.0d0 * bhyd(j, 4) * exp( -12.745d0 * tvi) + 
     &                50.0d0 * bhyd(j, 5) * exp( -13.051d0 * tvi) + 
     &                72.0d0 * bhyd(j, 6) * exp( -13.217d0 * tvi)
            d1      = 13.595d0 / 6.5d0 / 6.5d0 * tvi
            d2      = potlo(1) * tvi
*
! Helium 
         else if (n .eq. 3) then 
             b=1.0d0 
             part(1)=b
          
             if(t(j).lt.15000.) then 
               skip = .true.
             else 
               do  i=2,29
                b= 1.0d0 
                part(1)=part(1)+ghe1(i)*b  *exp(-ehe1(i)*hckt(j))
               end do 

               d1=109677.576/5.5/5.5*hckt(j)
               skip = .false.
             end if 
! Helium
         else if (n .eq. 4) then 
            b=1.
            part(2)=2.*b

            do  i=2,6
              b=1.
              part(2)=part(2)+ghe2(i)*b  *exp(-ehe2(i)*hckt(j))
            end do  
            d1=4.*109722.267/6.5/6.5*hckt(j)
            skip = .false.  

! Carbon 
         else if (n .eq. 354) then 
           b=1.
           part(1)=b*(1.+3.*exp(-16.42*hckt(j))+5.*exp(-43.42*hckt(j)))

           do  i=2,14
             b=1.
             part(1)=part(1)+gc1(i)*b  *exp(-ec1(i)*hckt(j))
           end do  

           part(1)=part(1)+108.*exp(-80000.*hckt(j))+
     &             189*exp(-84000.*hckt(j))+247.*exp(-87000.*hckt(j))+
     &             231.*exp(-88000.*hckt(j))+190.*exp(-89000.*hckt(j))+
     &             300.*exp(-90000.*hckt(j))

           skip = .true. 
! C 
         else if (n .eq. 355) then
            b=1.
            part(2)=b*(2.+4.*exp(-63.42*hckt(j)))
    
            do  i=2,6 
              b=1.
              part(2)=part(2)+gc2(i)*b*exp(-ec2(i)*hckt(j))
            end do 
 
            part(2)=part(2)+6.*exp(-131731.80*hckt(j))+
     & 4.*exp(-142027.1*hckt(j))+10.*exp(-145550.13*hckt(j))+
     & 10.*exp(-150463.62*hckt(j))+2.*exp(-157234.07*hckt(j))+
     & 6.*exp(-162500.*hckt(j))+42.*exp(-168000.*hckt(j))+
     & 56.*exp(-178000.*hckt(j))+102.*exp(-183000.*hckt(j))+
     & 400.*exp(-188000.*hckt(j))
     
            skip = .true. 

! magnesium
         else if (n .eq. 51) then

            b=1.
            part(1)=b

            do  i=2,11
               b=1.
               part(1)=part(1)+gmg1(i)*b   *exp(-emg1(i)*hckt(j))
            end do

            part(1)=part(1)+5.*exp(-53134.*hckt(j))+15.*
     &              exp(-54192.*hckt(j))+
     &              28.*exp(-54676.*hckt(j))+9.*exp(-57853.*hckt(j))
            g=4.
            d1=109734.83/4.5/4.5*hckt(j)
            skip = .false.
! Magnesium

        else if (n .eq. 52) then

            b=1.
            part(2)=b*2.

            do i=2,6
               b=1.
               part(2)=part(2)+gmg2(i)*b        *exp(-emg2(i)*hckt(j))
            end do

            part(2)=part(2)+10.*exp(-93310.80*hckt(j))+
     &          14.*exp(-93799.70*hckt(j))+6.*exp(-97464.32*hckt(j))+
     &          10.*exp(-103419.82*hckt(j))+14.*exp(-103689.89*hckt(j))+
     &          18.*exp(-103705.66*hckt(j))
            g=2.
            d1=4.*109734.83/5.5/5.5*hckt(j)
            skip = . true.

! Aluminium

         else if (n .eq. 57) then 
            b=1.
            part(1)=b *(2.+4.*exp(-112.061*hckt(j)))

            do  i=2,9
              b=1.
              part(1)=part(1)+gal1(i)*b *exp(-eal1(i)*hckt(j))
            end do

            part(1)=part(1)+10.*exp(-42235.*hckt(j))+
     &              14.*exp(-43831.*hckt(j))
            g=2.
            d1=109735.08/5.5/5.5*hckt(j)
            skip = .false.
! Si 
         else if (n .eq. 63) then 

            b=1.
            part(1)=b  *(1.+3.*exp(-77.115*hckt(j))+
     &              5.*exp(-223.157*hckt(j)))

            do  i=2,11
               b=1.
               part(1)=part(1)+gsi1(i)*b *exp(-esi1(i)*hckt(j))
            end do 

            part(1)=part(1)+76.*exp(-53000.*hckt(j))+ 
     &              71.*exp(-57000.*hckt(j))+
     &              191.*exp(-60000.*hckt(j))+240.*exp(-62000.*hckt(j))+
     &              251.*exp(-63000.*hckt(j))+300.*exp(-65000.*hckt(j))
            skip = .true. 

! Si 
         else if (n .eq. 64) then
            b=1.
            part(2)=b*(2.+4.*exp(-287.32*hckt(j)))

            do  i=2,6
               b=1.
               part(2)=part(2)+gsi2(i)*b *exp(-esi2(i)*hckt(j))
            end do 

            part(2)=part(2)+6.*exp(-81231.59*hckt(j))+
     &           6.*exp(-83937.08*hckt(j))+10.*exp(-101024.09*hckt(j))+
     &           14.*exp(-103556.35*hckt(j))+10.*exp(-108800.*hckt(j))+
     &           42.*exp(-115000.*hckt(j))+6.*exp(-121000.*hckt(j))+
     &           38.*exp(-125000.*hckt(j))+34.*exp(-132000.*hckt(j))
            g=2.
            d1=4.*109734.83/4.5/4.5*hckt(j)
            skip = .false. 
! Oxygen
         else if (n .eq. 367) then 
            b=1.
           part(1)=b*(5.+3.*exp(-158.265*hckt(j))+exp(-226.977*hckt(j)))

            do i=2,13
              b=1.
              part(1)=part(1)+go1(i)*b*exp(-eo1(i)*hckt(j))
            end do 
 
            part(1)=part(1)+15.*exp(-101140.*hckt(j))+
     &           131.*exp(-103000.*hckt(j))+128.*exp(-105000.*hckt(j))+
     &           600.*exp(-107000.*hckt(j))
            
            skip = .true. 
! Natrium
         else if (n .eq. 45 ) then 

            b=1.
            part(1)=b*2.

            do  i=2,8
             b=1.
             part(1)=part(1)+gna1(i)*b*exp(-ena1(i)*hckt(j))
            end do  

            part(1)=part(1)+10.*exp(-34548.745*hckt(j))+
     &               14.*exp(-34586.96*hckt(j))
            g=2.
            d1=109734.83/4.5/4.5*hckt(j)
            skip = .false. 
! B - Bor
         else if (n .eq. 14 ) then 

            b=1.
            part(1)=b*(2.+4.*exp(-15.25*hckt(j)))

            do  i=2,7
             b=1.
             part(1)=part(1)+gb1(i)*b*exp(-eb1(i)*hckt(j))
            end do 

            part(1)=part(1)+6.*exp(-57786.80*hckt(j))+
     &           10.*exp(-59989.*hckt(j))+14.*exp(-60031.03*hckt(j))+
     &           2.*exp(-63561.*hckt(j))
            g=2.
            d1=109734.83/4.5/4.5*hckt(j)
            skip = .false. 


! Kalium 
         else if (n .eq. 91 ) then 
   
            b=1.
            part(1)=b*2.

            do i=2,8
              b=1.
             part(1)=part(1)+gk1(i)*b*exp(-ek1(i)*hckt(j))
            end do 
  
            part(1)=part(1)+10.*exp(-27397.077*hckt(j))+
     &              14.*exp(-28127.85*hckt(j))
            g=2.
            d1=109734.83/5.5/5.5*hckt(j)
            skip = .false. 

         else
            t2000  = ip(ion) * 2000.0d0 / 11.0d0
            it     = max(1, min(9,  int(t(j) / t2000 - 0.5d0)))
            dt     = t(j) / t2000 - dble(it) - 0.5d0
            pmin   = 1.0d0
            i      = (it + 1) / 2
            k1     = nnn(i, n) / 100000
            k2     = nnn(i, n) - k1 * 100000
            k3     = k2 / 10
            kscale = k2 - k3 * 10
*
            if(mod(it, 2) .eq. 0) then
               p1     = dble(k3) * scale(kscale)
               k1     = nnn(i + 1, n) / 100000
               kscale = mod(nnn(i + 1, n), 10)
               p2     = dble(k1) * scale(kscale)
*
            else
               p1 = dble(k1) * scale(kscale)
               p2 = dble(k3) * scale(kscale)
*
               if(dt .lt. 0.0d0 .and. kscale .le. 1)then
                  kp1 = p1
                  if(kp1 .eq. int(p2 + 0.5d0)) pmin = kp1
               end if
*
            end if
*
            part(ion) = max(pmin, p1 + (p2 - p1) * dt)
*
            if(   g       .eq. 0.0d0 .or. 
     &         potlo(ion) .lt. 0.1d0 .or. 
     &           t(j)     .lt. t2000 * 4.0d0) then
               skip = .true.
*
            else
               skip = .false.
               if(t(j) .gt. (t2000 * 11.0d0)) tv = (t2000 * 11.0d0) * 
     &                                              8.6171d-5
               tvi = 1.0d0 / tv
               d1 = 0.1d0 * tvi
            end if
*
         end if
*
         if(.not. skip) then
            d2        = potlo(ion) * tvi
            part(ion) = part(ion) + g * exp( - ip(ion) * tvi) * 
     &                  (sqrt(13.595d0 * z * z * tvi / d2) ** 3 * 
     &                  (1.0d0 / 3.0d0 + (1.0d0 - (0.5d0 + (1.0d0 /
     &                  18.0d0 + d2 / 120.0d0) * d2) * d2) * d2) - 
     &                   sqrt(13.595d0 * z * z * tvi / d1) ** 3 * 
     &                   (1.0d0 / 3.0d0 + (1.0d0 - (0.5d0 + (1.0d0 / 
     &                  18.0d0 + d1 / 120.0d0) * d1) * d1) * d1))
            tv  = tkev(j)
            tvi = 1.0d0 / tv
         end if
*
      end do
*
      if(mode1 .eq. 5) then
         answer(7, 1) = 0.0d0
*
         do ion = 1, nion
            answer(ion, 1)     = part(ion)
            answer(ion + 7, 1) = ip(ion) + answer(ion + 6, 1)
         end do
*
      else
*
         if(mode1 .ne. 3) then
            n  = n - nion2
            cf = 2.0d0 * 2.4148d15 * t(j) * sqrt(t(j)) / xne(j)
*
            do ion = 2, nion2
               n      = n + 1
               f(ion) = cf * part(ion) / part(ion - 1) * 
     &            exp(-min((ip(ion-1) - potlo(ion-1))/tv,75.0d0))
*
*.... THE CUTOFF AT EXP(-75.0) IS TO BE CONSISTENT WITH BOB'S VERSION.
*.... IT IS NOT NEEDED ON SUN IN DOUBLE PRECISION, WHERE UNDERFLOW 
*.... STARTS AT EXP(-709).
*
            end do
*
            f(1) = 1.0d0
            l    = nion2 + 1
*
            do ion = 2, nion2
               l    = l - 1
               f(1) = 1.0d0 + f(l) * f(1)
            end do
*
            f(1) = 1.0d0 / f(1)
*
            do ion = 2, nion2
               f(ion) = f(ion - 1) * f(ion)
            end do
*
         end if
*
         if(mode .ge. 10) then
*
            if(mode1 .eq. 1) then
*
               do ion = 1, nion
                  answer(j, ion) = f(ion) / part(ion)
               end do
*
            else if(mode1 .eq. 2) then
*
               do ion = 1, nion
                  answer(j, ion) = f(ion)
               end do
*
            else if(mode1 .eq. 3) then
*
               do ion = 1, nion
                  answer(j, ion) = part(ion)
               end do
*
            else if(mode1 .eq. 4) then
               answer(j, 1) = 0.0d0
*
               do ion = 2, nion2
                  answer(j, 1) = answer(j, 1) + f(ion) * dble(ion - 1)
               end do
*
            end if
*
         else
*
            if(mode1 .eq. 1) then
               answer(j, 1) = f(nion) / part(nion)
*
            else if(mode1 .eq. 2) then
               answer(j, 1) = f(nion)
*
            else if(mode1 .eq. 3) then
               answer(j, 1) = part(nion)
*
            else if(mode1 .eq. 4) then
               answer(j, 1) = 0.0d0
*
               do ion = 2, nion2
                  answer(j, 1) = answer(j, 1) + f(ion) * dble(ion - 1)
               end do
*
            end if
*
         end if
*
      end if
*
      end
*
*********** E N D   O F   S U B R O U T I N E   P F S A H A ************

      end module population

