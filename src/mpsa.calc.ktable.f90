 subroutine calc_ktab

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





!    get the T, P values from the array



  do i = 1, numt           
         itemp = itemp +1
         nrhox = numpres 
 
!    get the T, P values from the array

       do j = 1, nrhox 
          t(j) = 10.0d0**(tabt(i))
          p(j) = 10.0d0**(tabp(j))

!----------------------------------------
          tk(j)     = k * t(j)
          hckt(j)   = hc / tk(j)
          hkt(j)    = h / tk(j)
          tkev(j)   = k_ev * t(j)
          tlog(j)   = log(t(j))
          xnatom(j) = p(j) / tk(j) - xne(j)
          rho(j)    = xnatom(j) * wtmole * 1.660d-24
          if (ifturb) pturb(j) = 0.5 * rho(j) * vturb(j) ** 2

       end do

!------- solve the equilibrium numbers,...
            call calc_se


!    --- find the frequency that you need ! if not already found in the mpsa.read routine 
         
            nu =  nulo
            freq = freset (nu)
            rco  = rcoset (nu)


            freqlg = log10(freq)
            freqln = log(freq)
            freq15 = freq * 1.0d-15
            waveno = freq / c_cm
            wave =  2.997925d17 / freq
            print*,' wave = ', wave !
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
               abtot(j)  = acont(j) + aline(j) + sigmac(j) +sigmal(j) 
               sigmal(j) = siglin(j) + sigxl(j)
           end do
!


!
           sumwt = 0.0d0
!
!.... METALLIC LINE BLANKETING ----- Juts the highest sub - bin 
            n = 1
!
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
                           sigmal(j) = siglin(j) + sigxl(j)

                         abtot(j) = acont(j)  + aline(j) + sigmac(j) + sigmal(j)

                       end do


       
            do j = 1, nrhox
              tabrosskap(i,j) = nint(log10(abtot(j))*1000.0)
            end do 



       end do ! over all full atmospheres with T,P points


      open(unit = 9 , file='kappa_table.dat', form = 'formatted'  ) 
! ---- at this point the table is calculated .
      write(9,*) numt, numpres
      write(9,*) (tabt(i), i = 1 , numt)
      write(9,*) (tabp(i), i = 1, numpres) 
      do i = 1, numt
         write(9,*) (tabrosskap(i,j), j = 1, numpres )
      end do
          
      close(unit=9)


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



 end subroutine calc_ktab 



