 subroutine calcross

    use types
    use atlcomm
    use kappa_cal
    use population

    implicit none

      include 'common.odfnlte'
      include 'common.constb'
      include 'common.abross'

      include 'common.abtotb'
      include 'common.contbl'
      include 'common.convbl'
      include 'common.depart'
      include 'common.edenbl'
      include 'common.elembl'
      include 'common.fluxbl'
      include 'common.height'
      include 'common.ifblkk'
      include 'common.ifopbl'
      include 'common.ionsbl'
      include 'common.iterbl'
      include 'common.junkbl'
      include 'common.musblk'
      include 'common.opsblk'
      include 'common.optotb'
      include 'common.putblk'
      include 'common.stateb'
      include 'common.steplg'
      include 'common.taushj'
      include 'common.tcorrb'
      include 'common.teffbl'
      include 'common.tempbl'
      include 'common.turbpr'
      include 'common.waveyb'
      include 'common.xabund'
      include 'common.xnfblk'
      include 'common.xnfpbl'
      include 'common.xnmolb'
      include 'common.freebl'

 
      integer i, l, n, nu, nmod, savenrhox, savei 
      integer j, kk
      integer nsteps

      integer numppoints 
      integer saveitemp 
      double precision savet(maxd), savep(maxd), srhox(maxd), sprad(maxd), sxne(maxd), stauros(maxd) 
      double precision sabross(maxd), spradk(maxd)   
      double precision tl, plog 

      integer    ip, ipj(maxd), it, itj(maxd)
      integer    countpairs
      integer    pbot(maxd), ptop(maxd), tarr(maxd), tcount
      integer    matrixtabs(maxd, 2) 
      integer    startj, stopj, imod

      logical  finish, cast_version

!--- local from main:       
      double precision  contin, freq15, rco, rcowt, stepwt, sumwt, wave, x
      double precision part(maxd, 6)


!----
      external  ross, ross2 , ross3

      saveitemp = itemp
      itemp = 101
      cast_version = .true. 
!--- set allowance of pressure points to number of pressure in ODF

      numppoints = numpres

      allocate(tabrosskap(numt, numpres))


      savenrhox = nrhox
 
      do j = 1, nrhox
       savet(j) = t(j)
       savep(j) = p(j)
       srhox(j) = rhox(j)
       sxne(j) = xne(j) 
       sprad(j) = prad(j) 
       spradk(j) = pradk(j)
       tlog(j) = log(t(j))
       stauros(j) = tauros(j) 
       sabross(j) = abross(j) 
      end do 




! first find all tj's and pj's :

      do j = 1, nrhox 
          tl = min(max(tlog(j) / tenlog, tabt(1)), tabt(numt) )
          it = 2
!
          do while (tabt(it) .le. tl .and. it .lt. numt)
             it = it + 1
          end do
          itj(j+1) = it

           plog = min(max(log10(p(j)), tabp(1)), tabp(numpres))
           ip   = 2
!
           do while (tabp(ip) .le. plog .and. ip .lt. numpres)
              ip = ip + 1
           end do
!
           ipj(j+1) = ip
           pbot(j+1) = max(1,ip-12)
           ptop(j+1) = min(numpres, ip+12)
!
       end do

!--- take also a temperature that is below the smallest T-value



         countpairs = 0 
         tcount = 1
         
         itj(1) = max(itj(2)-5, 1)  
         tarr(1) = itj(1)
         pbot(1) = pbot(2)        
         ptop(1) = ptop(2) 

! and one T values larger
         itj(nrhox+2) = itj(nrhox+1)+1
         pbot(nrhox+2) = pbot(nrhox+1)
         ptop(nrhox+2) = ptop(nrhox+1)

!  second get tarr with all j's and count them! 

         do j = 2, nrhox+2
           
           if (itj(j) .eq. tarr(tcount) ) then 
             pbot(tcount) = min(pbot(tcount), pbot(j))
             ptop(tcount) = max(ptop(tcount), ptop(j)) 
           else
! if not equal it also could be that the itj is not only one point away so 
             startj = tarr(tcount) +1
             stopj = itj(j)
             do i = startj, stopj
 
              tcount = tcount +1
              tarr(tcount) = i 
              pbot(tcount) = min(pbot(tcount-1), pbot(j) )
              ptop(tcount) = max(ptop(tcount-1), ptop(j))

              countpairs = countpairs + (ptop(tcount-1) -pbot(tcount-1)+1) 
             end do 

           end if 

         end do 
        
         countpairs = countpairs +(ptop(tcount) -pbot(tcount)+1) 

         nmod = ceiling( dble(countpairs)/dble(maxd-numppoints))

         do j = 1, tcount 
   
         end do  
!    get the T, P values from the array

         j = 1 
         savei = 1 


       do imod = 1, nmod          
         itemp = itemp +1
         nrhox = maxd 
         if (savei .le. tcount) then    
!    get the T, P values from the array
         
         i = savei 
         j = 1
         do while ((j.le. (nrhox-numppoints)) .and. i .le. tcount)

           kk = ptop(i)-pbot(i) +1
           if (kk .gt. numppoints) print*, 'we got not enough p-points allowance,', kk 

           do l = 0, kk-1
             t(j) = 10.0d0**(tabt(tarr(i)))
             p(j) = 10.0d0**(tabp(pbot(i)+l))
             matrixtabs(j,1) = tarr(i) 
             matrixtabs(j,2) = pbot(i) + l 
             j = j +1 
           end do
           i = i +1

         end do
         savei = i 
         nrhox = j-1



       do j = 1, nrhox
          tk(j)     = k * t(j)
          hckt(j)   = hc / tk(j)
          hkt(j)    = h / tk(j)
          tkev(j)   = k_ev * t(j)
          tlog(j)   = log(t(j))
          xnatom(j) = p(j) / tk(j) - xne(j)
          rho(j)    = xnatom(j) * wtmole * 1.660d-24
          if (ifturb) pturb(j) = 0.5 * rho(j) * vturb(j) ** 2

       end do


            call calc_se

! initialise freq intervals
            call ross


! start frequency interval 

!            
            do nu = nulo, nuhi
!------------ if wave open ---------------------------------------
!----------------------------------------------------------------*
                  if (ifwave) then
!
                     if (wbegin .le. 1.0d10) then
                        wave = wbegin + dble(nu-nulo) * deltaw
                        freq = 2.997925d17 / wave
                        rco  = abs (deltaw / wave * freq)
!*
                     else
!
!.... EQUALLY SPACED FREQUENCIES
!
                        freq = wbegin + dble(nu-nulo) * deltaw
                        rco  = deltaw
                     end if
!
                  else
                     freq = freset (nu)
                     rco  = rcoset (nu)


                  end if
!-------------------------------------------------------------------*
! --------  close if wave ------------------------------------------


                  freqlg = log10(freq)
                  freqln = log(freq)
                  freq15 = freq * 1.0d-15
                  waveno = freq / c_cm
                  wave =  2.997925d17 / freq
!
                  do j = 1, nrhox
                     ehvkt(j)  = exp(-freq * hkt(j))
                     stim(j)   = 1.0d0 - ehvkt(j)
                     bnu(j)    = 1.47439d-2 * freq15**3 * ehvkt(j) / stim(j)
                     dbnudt(j) = bnu(j) * freq * hkt(j) / t(j) / stim(j)
                     if(numnu .eq. 1) dbnudt(j) = 4.0d0 * sigma / pi *t(j) ** 3
                  end do

!  first calculate continuum at the frequency 

                  call kapp

!.... THIS NEXT LOOP IS NEEDED TO INITIALIZE ALINE, ... EACH NU
!
                     do j = 1, nrhox
                        aline(j)  = 0.0d0 

                        sline(j)  = bnu(j)
                         if (cast_version ) then
                         if (aline(j) .gt. 0.0d0) sline(j) = & 
     &                       (ahline(j) * shline(j) + & 
     &                        alines(j) * bnu(j) + & 
     &                        axline(j) * sxline(j)) / aline(j)
                         end if


                        sigmal(j) = siglin(j) + sigxl(j)
                     end do
!

                     do j = 1, nrhox
                       abtot(j) = acont(j)  + aline(j) + sigmac(j) + sigmal(j)
                       alpha(j) = (sigmac(j) + sigmal(j)) / abtot(j)
                     end do

!
!  ifsurf = 0 : 
!
                     do j = 1, nrhox
                           abtotc(j)  = abtot(j)
                           alphac(j)  = alpha(j)
                           residc(j)  = 0.0d0
                     end do
!
!
                     sumwt = 0.0d0
!
!.... METALLIC LINE BLANKETING
                    finish = .false.
                     n = 0
!
                     do while (.not. finish)
                        n = n + 1
!------ --  here are the ODFs read (three different possibilities:
!          a) once and, b) from memory c) several files if microtrubulence varies
!           throughout atmosphere
!- --------------------------------------------------------------------------
!
                        if (ifop(15) .and. lodf .eq. 'constant') then

                           if (ifkbin ) then
                             ! old routines used for mpsa.odf in binary
                             ! for that ifkbin has to be set ture

                             if (freqid .eq. 'little') then
                              call linopl ((n), nsteps, stepwt)
                             else
                              call linopb ((n), nsteps, stepwt)
                             end if
                           else
                             ! .nc files with defined grid were already read in mpsa.read
                             ! here we only interpolate on the T,P grid

                             call interPT(n, stepwt)
                             nsteps = isubbin

                           endif

                        else if (ifop(15) .and. lodf .eq. 'memory  ') then
                           call linopm ((n), nsteps, stepwt)
!
                        else if (ifop(15) .and. lodf .eq. 'variable') then
                           call linopv ((n), nsteps, stepwt)
                        end if

                        do j = 1, nrhox
                           aline(j)  = alines(j)
                           sline(j)  = bnu(j)

                           if (cast_version ) then
                             if (aline(j) .gt. 0.0d0) sline(j) = & 
     &                       (ahline(j) * shline(j) +  & 
     &                        alines(j) * bnu(j) + & 
     &                        axline(j) * sxline(j)) / aline(j)
                           end if

                           sigmal(j) = siglin(j) + sigxl(j)
                        end do
!

                       do j = 1, nrhox
                         abtot(j) = acont(j)  + aline(j) + sigmac(j) + sigmal(j)
                         alpha(j) = (sigmac(j) + sigmal(j)) / abtot(j)
                       end do

                        sumwt = sumwt + stepwt
                        rcowt = rco * stepwt
                        residc(1) = 1.0d0
!
!
!
                        if (stepwt .gt. 0.0d0 ) then 
                           call ross2 ((rcowt))
                        end if


                        put  = stepwt
                        iput = nsteps
!
!.... CHECK TO SEE IF THIS FREQUENCY IS FINISHED
!

                        if (n .eq. nsteps) then
                           finish = .true.
                           n = 0
                        end if


                    end do
!*
!*.... FINISH OFF THIS FREQUENCY INTERVAL
!*
                       residc(1) = 1.0d0
                       stepwt    = 1.0d0 - sumwt
                       if (stepwt .lt. 0.0001d0) stepwt = 0.0d0
!
! ifsurf .eq. 0 : *
                        do j = 1, nrhox
                           abtot(j)  = abtotc(j)
                           alpha(j)  = alphac(j)
                        end do
!
                     sumwt = sumwt + stepwt
                     rcowt = rco * stepwt
!
                     if (stepwt .gt. 0.0d0 ) then
                        call ross2 ((rcowt))
                     end if
                     put  = stepwt
                     iput = nsteps




            end do ! close frequency loop 


!*.... FINISH OFF THIS ITERATION
                  call ross3
            do j = 1, nrhox
              kk = matrixtabs(j,1)
              l = matrixtabs(j,2) 
              tabrosskap(kk,l) = nint(log10(abross(j))*1000.)
            end do 
              

         end if ! if we still have t,p values to do
       end do ! over all full atmospheres with T,P points





! give the saved values back : 
      nrhox = savenrhox 
      t = 0.0d0
      p = 0.0d0 
      rhox = 0.0d0
      xne = 0.0d0 
      prad = 0.0d0

      do j = 1, nrhox
       t(j) = savet(j)
       p(j) = savep(j)
       rhox(j) = srhox(j)
       xne(j) = sxne(j)
       prad(j) = sprad(j)
       pradk(j) = spradk(j)
       tauros(j) = stauros(j) 
       abross(j) = sabross(j)
      end do
      itemp =saveitemp

      do j = 1, nrhox
          tk(j)     = k * t(j)
          hckt(j)   = hc / tk(j)
          hkt(j)    = h / tk(j)
          tkev(j)   = k_ev * t(j)
          tlog(j)   = log(t(j))
          xnatom(j) = p(j) / tk(j) - xne(j)
          rho(j)    = xnatom(j) * wtmole * 1.660d-24

      end do



 end subroutine calcross


 subroutine inter_ross(rosst) 
     use types
     use atlcomm

      implicit none

!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE Rosstab has been calculated 
!.... ONLY FOR THAT VTURB
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision, intent(out) ::  rosst(maxd) 
      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd)
      double precision  plog, tl,  x, y

      integer    j, ip, ipj(maxd), it, itemp1, itj(maxd)
      save       co1, co2, co3, co4,  ipj, itj
!
!---------------------------------- EXECUTION --------------------------
!
         do j = 1, nrhox
            tl = min(max(tlog(j) / tenlog, tabt(1)), tabt(numt) )
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. numt)
               it = it + 1
            end do
!
            plog = min(max(log10(p(j)), tabp(1)), tabp(numpres))
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. numpres)
               ip = ip + 1
            end do
!
            ipj(j) = ip
            itj(j) = it
            x      = (tl   - tabt(it-1)) / (tabt(it) - tabt(it-1))
            y      = (plog - tabp(ip-1)) / (tabp(ip) - tabp(ip-1))
            co1(j) = (1.0d0 - x) * (1.0d0 - y) * ttenlg
            co2(j) = (1.0d0 - x) * y * ttenlg
            co3(j) = x * (1.0d0 - y) * ttenlg
            co4(j) = x * y * ttenlg
!.... THE STEPS HAVE BEEN SCALED BY 1000
         end do
!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         a  = exp( co1(j) * dble(tabrosskap( it-1, ip-1)) + & 
     &             co2(j) * dble(tabrosskap( it-1  , ip)) + & 
     &             co3(j) * dble(tabrosskap(it, ip-1)) + & 
     &             co4(j) * dble(tabrosskap( it  , ip)) )
         rosst(j) = a
      end do
 end subroutine


  double precision function introssk( tv, pv)

     use types
     use atlcomm

      implicit none

!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE Rosstab has been calculated
!.... ONLY FOR THAT VTURB
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision tv, pv

      double precision  a, co1, co2, co3, co4
      double precision  plog, tl,  x, y

      integer    j, ip, it
!
!---------------------------------- EXECUTION --------------------------
!
            tl = min(max(log10(tv) , tabt(1)), tabt(numt) )
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. numt)
               it = it + 1
            end do
!
            plog = min(max(log10(pv), tabp(1)), tabp(numpres))
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. numpres)
               ip = ip + 1
            end do
!
            x      = (tl   - tabt(it-1)) / (tabt(it) - tabt(it-1))
            y      = (plog - tabp(ip-1)) / (tabp(ip) - tabp(ip-1))
            co1 = (1.0d0 - x) * (1.0d0 - y) * ttenlg
            co2 = (1.0d0 - x) * y * ttenlg
            co3 = x * (1.0d0 - y) * ttenlg
            co4 = x * y * ttenlg
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         a  = exp( co1 * dble(tabrosskap( it-1, ip-1)) + &
     &             co2 * dble(tabrosskap( it-1  , ip)) + &
     &             co3 * dble(tabrosskap(it, ip-1)) + &
     &             co4 * dble(tabrosskap( it  , ip)) )
         introssk= a
 end function 


