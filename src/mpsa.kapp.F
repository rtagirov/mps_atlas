      module kappa_cal 
      use types
      use atlcomm

      implicit none

      include 'common.rhoxbl'
      include 'common.freqbl'


      Contains


      subroutine kapp

      implicit none
!!
!!.... 1985 APR - MODIFIED TO INCLUDE A NEW HLINOP AND TO REORGANIZE
!                KAPP BY REMOVING ALL THE LINE OPACITY CALLS
!.... 1988 MAR - MODIFIED LINOP TO HANDLE EXTRAPOLATION OFF THE ODF
!.... 1991 JUN - MODIFIED TO BRING INTO LINE WITH ATLAS9.  MADE THREE
!                VERSIONS OF LINOP FOR THE THREE POSSIBILITIES:
!                LINOPC = A SINGLE, CONSTANT TURBULENT VELOCITY
!                LINOPM = A SINGLE, CONSTANT TURBULENT VELOCITY IN MEMORY
!                LINOPV = VARIABLE TURBULENT VELOCITY
!.... 1992 APR - CHANGE TO DOUBLE PRECISION
!.... 1993 JUN - CHANGED HMINOP TO LINEAR INTERPOLATION
!.... 1994 JUL - FIXED ERROR WITH INITIALIZATION USING DUMMY - 
!                CHANGED START OF EQUIVALENCE FROM AHYD(1) TO AAL1(1)
!                CHANGED SIZE OF DUMMY FROM 25 TO 37
!                CHANGED DO LIMIT FROM 25 TO 37
!.... 1995 MAR - REMOVED waveno CALCULATIONS FROM ALL ROUTINES.
!                waveno IS CALCULATED IN THE MAIN PROGRAM AND PASSED 
!                THROUGH common.freqbl.  THIS FORCED SOME CHANGES
!                IN THE USE OF freq AS A DUMMY VARIABLE.
!.... 1995 JUL - ADDING common.freqbl AS DESCRIBED ABOVE MADE THE
!                VARIABLE freqlg VISIBLE IN ROUTINES he12s1s, he12s3s,
!                he12p1p, he12p3p.  IN EACH OF THESE ROUTINES, freqlg
!                IS REDEFINED TO BE LOG10, BUT THIS IS PASSED OUT IN
!                common.freqbl.  
!              - SOLUTION: IN MAIN, CHANGED freqlg TO log10(freq)
!                AND CREATED freqln = log(freq), WHICH IS ADDED TO 
!                common.freqbl.  CHANGED ALL USES OF freqlg IN OPACITY
!                ROUTINES TO BE CONSISTENT.
!.... 1995 OCT - REVERSED TESTS IN SEVERAL ROUTINES TO TEST THE INDEX
!                FIRST TO AVOID ARRAY OUT OF BOUNDS
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'


      include 'common.ifopbl'
      include 'common.ionsbl'
      include 'common.opsblk'
      include 'common.optotb'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!-------------------------------- EXTERNALS ---------------------------
!
!      external coolop, elecop, h2plop, h2raop, he1op,  he2op,  hemiop, 
!     &         heraop, hlinop, hminop, hop,    hotop,  hrayop, linsop, 
!     &         lukeop, pops,   xconop, xlinop, xlisop, xsop
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  a, dummy(maxd, 38)
      integer           i, j
!
!----------------------------- EQUIVALENCE ----------------------------
!
      equivalence (aal1(1), dummy(1, 1))
!
!------------------------------- EXECUTION ----------------------------
!
!.... EACH CALL TO KAPP USES DUMMY TO INITIALIZE ALL OPACITIES TO 0.0
!
      do i = 1, 38
!
         do  j = 1, nrhox
            dummy(j, i) = 0.0d0
         end do
!
      end do
!
      if (ifop(1)) call hop
      if (ifop(2) .and. freq .le. 3.28805d15) call h2plop
      if (ifop(3)) call hminop
      if (ifop(4)) call hrayop
      if (ifop(5)) call he1op
      if (ifop(6)) call he2op
      if (ifop(7)) call hemiop
      if (ifop(8)) call heraop
      if (ifop(9) .and. freq .le. 3.28805d15) call coolop
#ifdef NHOP
      if (ifop(9) .and. freq .le. 5.9958d15) call nhop
#endif
      if (ifop(10)) call lukeop
      if (ifop(11)) call hotop
      if (ifop(12)) call elecop
      if (ifop(13)) call h2raop
      if (ifop(19)) call xconop
      if (ifop(20)) call xsop
!
      do j = 1, nrhox
         a         = ah2p(j) + ahe1(j) + ahe2(j) + ahemin(j) + 
     &               acool(j) + aluke(j) + ahot(j) 
#ifdef NHOP
     & + anhop(j)
#endif 
 
         acont(j)  = a + ahyd(j) + ahmin(j) + axcont(j)
         scont(j)  = bnu(j)
         if (acont(j) .gt. 0.0d0) scont(j) = 
     &               (a * bnu(j) + ahyd(j) * shyd(j) + 
     &                ahmin(j)* shmin(j) + axcont(j) * sxcont(j)) / 
     &                acont(j)
         sigmac(j) = sigh(j) + sighe(j) + sigel(j) + sigh2(j) +sigx(j)


      end do            





!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   K A P P !!!!!!!!!!!!!!!!
!
#ifdef DFSY

      SUBROUTINE HOP

      IMPLICIT REAL(kind=8) (A-H,O-Z)
      integer kw, j 
      PARAMETER (kw=99)
!     REQUIRES FUNCTION COULFF
!      include 'common.sizebl'
      include 'common.constb'
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'


      FREQ3=2.815D29/FREQ/FREQ/FREQ
      DO 100 J=1,NRHOX
!     LIMIT 109678.764   RYDBERG 109677.576
!  N=16 TO INFINITY
      Htemp=FREQ3*2./2./(109677.576*HCKT(J))*
     1(EXP(-MAX(109250.336D0,109678.764D0-WAVENO)*HCKT(J))-
     2EXP(-109678.764*HCKT(J)))*STIM(J)
      S=Htemp*BNU(J)
!  N=15
      IF(WAVENO.LT.    487.456)GO TO 30
      X=XKARSAS(1.D0,15,15)
      A=X*450.*EXP(-   109191.313*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=14
      IF(WAVENO.LT.    559.579)GO TO 30
      X=XKARSAS(1.D0,14,14)
      A=X*392.*EXP(-   109119.188*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=13
      IF(WAVENO.LT.    648.980)GO TO 30
      X=XKARSAS(1.D0,13,13)
      A=X*338.*EXP(-   109029.789*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=12
      IF(WAVENO.LT.    761.649)GO TO 30
      X=XKARSAS(1.D0,12,12)
      A=X*288.*EXP(-   108917.117*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=11
      IF(WAVENO.LT.    906.426)GO TO 30
      X=XKARSAS(1.D0,11,11)
      A=X*242.*EXP(-   108772.336*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=10
      IF(WAVENO.LT.   1096.776)GO TO 30
      X=XKARSAS(1.D0,10,10)
      A=X*200.*EXP(-   108581.992*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=9
      IF(WAVENO.LT.  1354.044)GO TO 30
      X=XKARSAS(1.D0,9,9)
      A=X*162.*EXP(-  108324.719*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=8
      IF(WAVENO.LT.1713.713)GO TO 30
      X=XKARSAS(1.D0,8,8)
      A=X*128.*EXP(-107965.051*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=7
      IF(WAVENO.LT.2238.320)GO TO 30
      X=XKARSAS(1.D0,7,7)
      A=X*98.*EXP(-107440.444*HCKT(J))*STIM(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=6
      IF(WAVENO.LT.3046.604)GO TO 30
      X=XKARSAS(1.D0,6,6)
      A=X*72.*EXP(-106632.160*HCKT(J))*(BHYD(J,6)-EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHYD(J,6)-EHVKT(J))
!  N=5
      IF(WAVENO.LT.4387.113)GO TO 30
      X=XKARSAS(1.D0,5,5)
      A=X*50.*EXP(-105291.651*HCKT(J))*(BHYD(J,5)-EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHYD(J,5)-EHVKT(J))
!  N=4
      IF(WAVENO.LT.6854.871)GO TO 30
      X=XKARSAS(1.D0,4,4)
      A=X*32.*EXP(-102823.893*HCKT(J))*(BHYD(J,4)-EHVKT(J))
      Htmep=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHYD(J,4)-EHVKT(J))
!  N=3
      IF(WAVENO.LT.12186.462)GO TO 30
      X=XKARSAS(1.D0,3,3)
      A=X*18.*EXP(- 97492.302*HCKT(J))*(BHYD(J,3)-EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHYD(J,3)-EHVKT(J))
!  N=2
      IF(WAVENO.LT.27419.659)GO TO 30
      X=XKARSAS(1.D0,2,2)
      A=X* 8.*EXP(- 82259.105*HCKT(J))*(BHYD(J,2)-EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHYD(J,2)-EHVKT(J))
!  N=1
      IF(WAVENO.LT.109678.764)GO TO 30
      X=XKARSAS(1.D0,1,1)
      A=X* 2.*1. *(BHYD(J,1)-EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHYD(J,1)-EHVKT(J))
!C
   30 Htemp=Htemp*XNFPH(J,1)/RHO(J)
      S=S*XNFPH(J,1)/RHO(J)
!     FREE-FREE
      A=3.6919D8/SQRT(T(J))*COULFF(J,1)/FREQ*XNE(J)/FREQ*XNFPH(J,2)/
     1FREQ*STIM(J)/RHO(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
      AHYD(J)=Htemp
      IF(Htemp.GT.0.)SHYD(J)=S/Htemp
  100 CONTINUE
      RETURN
      END SUBROUTINE 



#else
      subroutine hop
      implicit none
!
!.... REQUIRES COULFF AND XKARSAS
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  bolt(maxd, 8), boltex(maxd), c, cfree, cont(8), 
     &                   ex, exlim(maxd), freet(maxd), freqi, 
     &                  hyd, s, xr
      save bolt, boltex, exlim, freet, xr
      integer itemp1, j, n
!
!---------------------------- INITIALIZATION --------------------------
!
      save itemp1
      data itemp1 / 0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do j = 1, nrhox
!
            do n = 1, 8
               bolt(j, n) = exp( -(hydip - hydip / dble(n*n) ) /
     &           tkev(j) ) * 2.0d0 * dble(n*n) * xnfph(j,1) / rho(j)
            end do
!
            do n = 1, 8
               bolt(j, n) = bolt(j, n) * bhyd(j, n)
            end do
!
!.... TAKE 1.0D10 OUT OF FREET AND PUT IT WITH CFREE
            freet(j) = 1.0d-10 * xne(j) * xnfh(j,2) /
     &                 (rho(j) * sqrt(t(j)))
            xr = xnfph(j, 1) * (2.0d0 / 2.0d0 / hydip) * tkev(j) / 
     &           rho(j)
            boltex(j) = exp(-13.427d0 / tkev(j) ) * xr
            exlim(j) = exp( -hydip / tkev(j) ) * xr
         end do
!
      end if
!
      do n = 1, 8
         cont(n) = xkarsas (1.0d0, n, n)
      end do
!
      freqi = 1.0d0 / freq
!
!.... PUT IN THE FACTOR OF 1.0D10 HERE FROM FREET(J)
!
      cfree = freqi * (freqi * (freqi * 3.6919d18 ) )
      c     = freqi * (freqi * (freqi * 2.815d29 ) )
!
      do j = 1, nrhox
         ex = boltex(j)
         if (freq .lt. 4.05933d13) ex = exlim(j) / ehvkt(j)
         hyd = (cont(7) * bolt(j, 7) + cont(8) * bolt(j, 8) +
     &       (ex-exlim(j)) * c + coulff(j, 1) * cfree * freet(j)) *
     &        stim(j)
         s = hyd * bnu(j)
!
         do n = 1, 8
            hyd = hyd + cont(n) * bolt(j, n) * (1.0d0 - ehvkt(j) / 
     &              bhyd(j, n))
            s = s + cont(n) * bolt(j, n) * bnu(j) * stim(j) / 
     &              bhyd(j, n)
         end do
!
         ahyd(j) = hyd
         shyd(j) = s / hyd
      end do
!
      end subroutine 

#endif

!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H O P !!!!!!!!!!!!!!!!!!!!! 
!
      double precision function coulx (n, freq, z)
      implicit none
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision freq, z
      integer n
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision a(6), b(6), c(6)
!
!
!------------------------------ INITIALIZATION ------------------------
!
      data a / 0.9916, 1.105, 1.101, 1.101,  1.102, 1.0986 /
      data b /  2.719d13, -2.375d14, -9.863d13, -5.765d13, -3.909d13,
     &         -2.704d13 /
      data c / -2.268d30, 4.077d28, 1.035d28, 4.593d27, 2.371d27,
     &          1.229d27 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (freq .lt. z * z * 3.28805d15 / dble(n * n)) then
         coulx = 0.0d0
!
      else
         coulx = 2.815d29 / freq/ freq/ freq/ dble(n ** 5) * z ** 4
!
         if (n .eq. 1) then
            coulx = coulx * coulbf1s (freq, z)
!
         else if (n .le. 6) then
            coulx = coulx * (a(n) + (b(n) + c(n) * (z * z / freq)) * 
     &              (z * z / freq))
!
         end if
!
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   C O U L X !!!!!!!!!!!!!!
!
      double precision function coulbf1s (freq, z)
      implicit none
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision freq, z
!
!---------------------------- LOCAL VARIABLES --------------------------
!
      double precision  elog, gaunt1s(151)
      integer i
!
!---------------------------- INITIALIZATION ---------------------------
!
      data gaunt1s /
     &   0.7973, 0.8094, 0.8212, 0.8328, 0.8439, 0.8548, 0.8653, 0.8754, 
     &   0.8852, 0.8946, 0.9035, 0.9120, 0.9201, 0.9278, 0.9351, 0.9420, 
     &   0.9484, 0.9544, 0.9601, 0.9653, 0.9702, 0.9745, 0.9785, 0.9820, 
     &   0.9852, 0.9879, 0.9903, 0.9922, 0.9938, 0.9949, 0.9957, 0.9960, 
     &   0.9960, 0.9957, 0.9949, 0.9938, 0.9923, 0.9905, 0.9884, 0.9859, 
     &   0.9832, 0.9801, 0.9767, 0.9730, 0.9688, 0.9645, 0.9598, 0.9550, 
     &   0.9499, 0.9445, 0.9389, 0.9330, 0.9269, 0.9206, 0.9140, 0.9071, 
     &   0.9001, 0.8930, 0.8856, 0.8781, 0.8705, 0.8627, 0.8546, 0.8464, 
     &   0.8381, 0.8298, 0.8213, 0.8128, 0.8042, 0.7954, 0.7866, 0.7777, 
     &   0.7685, 0.7593, 0.7502, 0.7410, 0.7318, 0.7226, 0.7134, 0.7042, 
     &   0.6951, 0.6859, 0.6767, 0.6675, 0.6584, 0.6492, 0.6401, 0.6310, 
     &   0.6219, 0.6129, 0.6039, 0.5948, 0.5859, 0.5769, 0.5680, 0.5590, 
     &   0.5502, 0.5413, 0.5324, 0.5236, 0.5148, 0.5063, 0.4979, 0.4896, 
     &   0.4814, 0.4733, 0.4652, 0.4572, 0.4493, 0.4415, 0.4337, 0.4261, 
     &   0.4185, 0.4110, 0.4035, 0.3962, 0.3889, 0.3818, 0.3749, 0.3680, 
     &   0.3611, 0.3544, 0.3478, 0.3413, 0.3348, 0.3285, 0.3222, 0.3160, 
     &   0.3099, 0.3039, 0.2980, 0.2923, 0.2866, 0.2810, 0.2755, 0.2701, 
     &   0.2648, 0.2595, 0.2544, 0.2493, 0.2443, 0.2394, 0.2345, 0.2298, 
     &   0.2251, 0.2205, 0.2160, 0.2115, 0.2072, 0.2029, 0.1987/
!
!-------------------------------- EXECUTION ----------------------------
!
      if(freq / z ** 2 .ge. 3.28805d15) then
         elog = log10(freq / z ** 2 /3.28805d15)
         i = elog / 0.02
         i = max(min(i+1, 150), 1)
         coulbf1s = gaunt1s(i) + (gaunt1s(i+1) - gaunt1s(i)) / 
     &              0.02d0 * (elog - (i-1) * 0.02d0)
      else
         coulbf1s = 0.0d0
      end if
!
      end function 
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   C O U L B F 1 S !!!!!!!!
!
      double precision function coulff(j, nz)
      implicit none
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j, nz
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'

      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision a(11, 12), gamlog, hvktlg, p, q, z4log(6)
      integer igam, ihvkt
!
!------------------------------ INITIALIZATION ------------------------
!
      data z4log / 0.0, 1.20412, 1.90849, 2.40824, 2.79588, 3.11216 /
!
      data a / 5.53, 5.49, 5.46, 5.43, 5.40, 5.25, 5.00, 4.69, 4.48,
     &         4.16, 3.85, 4.91, 4.87, 4.84, 4.80, 4.77, 4.63, 4.40,
     &         4.13, 3.87, 3.52, 3.27, 4.29, 4.25, 4.22, 4.18, 4.15,
     &         4.02, 3.80, 3.57, 3.27, 2.98, 2.70, 3.64, 3.61, 3.59,
     &         3.56, 3.54, 3.41, 3.22, 2.97, 2.70, 2.45, 2.20, 3.00,
     &         2.98, 2.97, 2.95, 2.94, 2.81, 2.65, 2.44, 2.21, 2.01,
     &         1.81, 2.41, 2.41, 2.41, 2.41, 2.41, 2.32, 2.19, 2.02,
     &         1.84, 1.67, 1.50, 1.87, 1.89, 1.91, 1.93, 1.95, 1.90,
     &         1.80, 1.68, 1.52, 1.41, 1.30, 1.33, 1.39, 1.44, 1.49,
     &         1.55, 1.56, 1.51, 1.42, 1.33, 1.25, 1.17, 0.90, 0.95,
     &         1.00, 1.08, 1.17, 1.30, 1.32, 1.30, 1.20, 1.15, 1.11,
!err fixed by Bob 13 april 1988
!!!! &         0.45, 0.48, 0.52, 0.60, 0.75, 0.91, 1.15, 1.18, 1.15,
     &         0.55, 0.58, 0.62, 0.70, 0.85, 1.01, 1.15, 1.18, 1.15,
     &         1.11, 1.08, 0.33, 0.36, 0.39, 0.46, 0.59, 0.76, 0.97,
     &         1.09, 1.13, 1.10, 1.08, 0.19, 0.21, 0.24, 0.28, 0.38,
     &         0.53, 0.76, 0.96, 1.08, 1.09, 1.09 /
!
!------------------------------- EXECUTION ----------------------------
!
      gamlog = 10.39638 - tlog(j) / 1.15129 + z4log(nz)
      igam   = max0 ( min0 ( int(gamlog + 7.0) , 10) , 1)
      hvktlg = (freqln - tlog(j) ) / 1.15129 - 20.63764
      ihvkt  = max0 ( min0 ( int(hvktlg + 9.0) , 11) , 1)
      p      = gamlog - dble(igam - 7)
      q      = hvktlg - dble(ihvkt - 9)
      coulff = (1.0 - p) * ((1.0 - q) * a(igam, ihvkt) + q * 
     &         a(igam, ihvkt + 1)) + p * ( (1.0 - q ) * 
     &         a(igam + 1, ihvkt) + q * a(igam + 1, ihvkt + 1) )
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   C O U L F F !!!!!!!!!!!!
!
      subroutine h2plop

      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  es, fr, fr15, xp2
      integer j
!
!------------------------------- EXECUTION ----------------------------
!
      if ( freq .le. 3.28805d15) then 

      fr15 = freq * 1.0d-15
      fr = -3.0233d3 + 
     &      (3.7797d2 + (-1.82496d1 + (3.9207d-1 - 3.1672d-3 *
     &       freqln) * freqln) * freqln) * freqln
      es = -7.342d-3 + 
     &      (-2.409 + (1.028 + (-4.230d-1 + (1.224d-1 - 1.351d-2 *
     &       fr15) * fr15) * fr15) * fr15) * fr15
!
      do j = 1, nrhox
         xp2 = 0.5 * (fr - es / tkev(j))
         ah2p(j) = exp(xp2) * xnfph(j, 1) * 2.0 * bhyd(j, 1) *
     &             xnfph(j, 2) * exp(xp2) / rho(j) * stim(j)
      end do

      else

         ah2p = 0.0d0

      end if

!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H 2 P L O P !!!!!!!!!!!!
!
      subroutine hminop
      implicit none
!-------------------------------- COMMONS -----------------------------

!      include 'common.sizebl'
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision bf(85), ff(11, 22), fflog(22, 11), fftheta(maxd), 
     &                 fftlog(1), fftt(11), h, hminbf, hminff, 
     &                 theta(maxd), 
     &                 thetaff(11), wave, wavek(22), wavelog(1), 
     &                 wbf(85), 
     &                 wfflog(22), xhmin(maxd)
      integer          i, itemp1, itheta, iwave, j, map1, maxwave
      logical          start
      save itemp1, start, fflog, theta, wfflog, xhmin
!
!------------------------------ EXTERNAL FUNCTION ----------------------
!
      external  map1
!
!------------------------------ INITIALIZATION ------------------------
!
      data itemp1, start / 0, .true. /
!
!.... FROM MATHISEN (1984), AFTER WISHART(1979) AND BROAD AND REINHARDT (1976)
!
      data bf/   0.067,   0.088,    0.117,    0.155,    0.206,  
     &           0.283,   0.414,    0.703,    1.24,     2.33,  
     &          11.60,   13.90,    24.30,    66.70,    95.00,
     &          56.60,   20.00,    14.60,     8.50,     7.10,   
     &           5.43,    5.91,     7.29,     7.918,    9.453,  
     &          11.08,   12.75,    14.46,    16.19,    17.92,  
     &          19.65,   21.35,    23.02,    24.65,    26.24,  
     &          27.77,   29.23,    30.62,    31.94,    33.17,  
     &          34.32,    35.37,    36.32,    37.17,    37.91,  
     &          38.54,   39.07,    39.48,    39.77,    39.95,  
     &          40.01,   39.95,    39.77,    39.48,    39.06,
     &          38.53,   37.89,    37.13,    36.25,    35.28,  
     &          34.19,   33.01,    31.72,    30.34,    28.87,  
     &          27.33,   25.71,    24.02,    22.26,    20.46,  
     &          18.62,   16.74,    14.85,    12.95,    11.07,  
     &           9.211,   7.407,    5.677,    4.052,    2.575,  
     &           1.302,   0.8697,   0.4974,   0.1989,   0.00 /
!
      data ((ff(i, j), i = 1, 11), j = 1, 11) /
     &.0178,.0222,.0308,.0402,.0498,.0596,.0695,.0795,.0896, .131, .172,   1823
     &.0228,.0280,.0388,.0499,.0614,.0732,.0851,.0972, .110, .160, .211,   2278
     &.0277,.0342,.0476,.0615,.0760,.0908, .105, .121, .136, .199, .262,   2604
     &.0364,.0447,.0616,.0789,.0966, .114, .132, .150, .169, .243, .318,   3038
     &.0520,.0633,.0859, .108, .131, .154, .178, .201, .225, .321, .418,   3645
     &.0791,.0959, .129, .161, .194, .227, .260, .293, .327, .463, .602,   4557
     &.0965, .117, .157, .195, .234, .272, .311, .351, .390, .549, .711,   5063
     & .121, .146, .195, .241, .288, .334, .381, .428, .475, .667, .861,   5696
     & .154, .188, .249, .309, .367, .424, .482, .539, .597, .830, 1.07,   6510
     & .208, .250, .332, .409, .484, .557, .630, .702, .774, 1.06, 1.36,   7595
     & .293, .354, .468, .576, .677, .777, .874, .969, 1.06, 1.45, 1.83/   9113
!
      data ((ff(i, j), i = 1, 11), j = 12, 22) /
     & .358, .432, .572, .702, .825, .943, 1.06, 1.17, 1.28, 1.73, 2.17,  10126
     & .448, .539, .711, .871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.60,  11392
     & .579, .699, .924, 1.13, 1.33, 1.51, 1.69, 1.86, 2.02, 2.67, 3.31,  13019
     & .781, .940, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48, 2.69, 3.52, 4.31,  15189
     & 1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.20, 3.51, 3.80, 4.92, 5.97,  18227
     & 1.73, 2.08, 2.74, 3.37, 3.90, 4.50, 5.01, 5.50, 5.95, 7.59, 9.06,  22784
     & 3.04, 3.65, 4.80, 5.86, 6.86, 7.79, 8.67, 9.50, 10.3, 13.2, 15.6,  30378
     & 6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2, 23.0, 29.5, 35.0,  45567
     & 27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117., 140.,  91134
     & 42.3, 50.6, 66.4, 80.8, 94.5, 107., 120., 131., 142., 183., 219., 113918
     & 75.1, 90.0, 118., 144., 168., 191., 212., 234., 253., 325., 388./ 151890
!
      data thetaff / 
     &     0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6 /
!
!.... BELL AND BERRINGTON J.PHYS.B,VOL. 20, 801-806,1987.
!
      data wavek /
     &   0.50,  0.40, 0.35, 0.30, 0.25, 0.20, 0.18, 0.16, 0.14, 0.12, 
     &   0.10,  0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 
     &   0.008, 0.006 /
!
      data wbf/  18.00,   19.60,   21.40,   23.60,   26.40,  
     &           29.80,   34.30,   40.40,   49.10,   62.60, 
     &          111.30,  112.10,  112.67,  112.95,  113.05,
     &          113.10,  113.20,  113.23,  113.50,  114.40, 
     &          121.00,  139.00,  164.00,  175.00,  200.00, 
     &          225.00,  250.00,  275.00,  300.00,  325.00, 
     &          350.00,  375.00,  400.00,  425.00,  450.00, 
     &          475.00,  500.00,  525.00,  550.00,  575.00, 
     &          600.00,  625.00,  650.00,  675.00,  700.00, 
     &          725.00,  750.00,  775.00,  800.00,  825.00, 
     &          850.00,  875.00,  900.00,  925.00,  950.00,
     &          975.00, 1000.00, 1025.00, 1050.00, 1075.00,
     &         1100.00, 1125.00, 1150.00, 1175.00, 1200.00,
     &         1225.00, 1250.00, 1275.00, 1300.00, 1325.00,
     &         1350.00, 1375.00, 1400.00, 1425.00, 1450.00,
     &         1475.00, 1500.00, 1525.00, 1550.00, 1575.00,
     &         1600.00, 1610.00, 1620.00, 1630.00, 1643.91/
!
!------------------------------- EXECUTION ----------------------------
!
      if (start) then
         start = .false.
!
         do iwave = 1, 22
!
!.... THE NUMBER 91.134 IS TAKEN FROM BELL AND BERRINGTON
!
            wfflog(iwave) = log (91.134d0 / wavek(iwave))
!
            do itheta = 1, 11
               fflog(iwave, itheta) = log (ff(itheta, iwave) * 1.0d-26)
            end do
!
         end do
!
      end if
!
      if (itemp .ne. itemp1 ) then
         itemp1 = itemp
!
         do j = 1, nrhox
            theta(j) = 5040.0d0 / t(j)
!
!....  0.754209 HOTOP AND LINEBERGER J.PHYS.CHEM.REF.DATA VOL 14,731-752,1985.
!
            xhmin(j) = exp(0.754209d0 / tkev(j)) / 
     &                 (2.0d0 * 2.4148d15 * t(j) * sqrt(t(j))) * 
     &                 bmin(j) * bhyd(j, 1) * xnfph(j, 1) * xne(j)
         end do
!
      end if
!
      wave    = 2.99792458d17 / freq
      wavelog(1) = log(wave)
!
      do itheta = 1, 11
         call linter(wfflog, fflog(1,itheta), 22, wavelog, fftlog, 1)
         fftt(itheta) = exp(fftlog(1)) / thetaff(itheta) * 5040.0d0 * 
     &                  1.380658d-16
      end do
!
      hminbf = 0.0d0
      if(freq .gt. 1.82365d14) maxwave = map1(wbf,bf,85,wave,hminbf,1)
! interpolation does not work for wave < 12.8
! to avoid h < 0
      if(hminbf .lt. 0.0d0) hminbf = 0.0d0
!
      do j = 1, nrhox
         call linter(thetaff, fftt, 11, theta(j), fftheta(j), 1)
         hminff = fftheta(j) * xnfph(j,1) * 2.0d0 * bhyd(j,1) * xne(j) /
     &            rho(j)
         h = hminbf * 1.0d-18 * (1.0d0 - ehvkt(j) / bmin(j)) * xhmin(j)/
     &       rho(j)
         ahmin(j) = h + hminff
         shmin(j) = (h * bnu(j) * stim(j) / (bmin(j) - ehvkt(j)) + 
     &               hminff * bnu(j)) / ahmin(j)

      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H M I N O P !!!!!!!!!!!!
!
      subroutine linter(xold, yold, nold, xnew, ynew, nnew)
      implicit none
!
!---------------------------- DUMMY VARIABLES  ------------------------
!.... xold and xnew increasing
!
      double precision xnew(*), xold(*), ynew(*), yold(*)
      integer          nnew, nold
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      integer inew, iold
!
!------------------------------- EXECUTION ----------------------------
!
      iold = 2
!
      do inew = 1, nnew
!
         do while(iold .lt. nold .and. xold(iold) .lt. xnew(inew))
            iold = iold + 1
         end do
!
         ynew(inew) = yold(iold-1) + (yold(iold) - yold(iold-1)) /
     &                               (xold(iold) - xold(iold-1)) * 
     &                               (xnew(inew) - xold(iold-1))
      end do
!
      end subroutine 
!
!!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L I N T E R !!!!!!!!!!!!!
!
!      subroutine hrayop
!      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.depart'
!      include 'common.freqbl'
!      include 'common.ionsbl'
!      include 'common.opsblk'
!      !include 'common.rhoxbl'
!      include 'common.stateb'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
!      double precision sig, wave, ww
!      integer j
!
!-------------------------------- EXECUTION ---------------------------
!
!.... SET THE CUTOFF AT THE LYMAN ALPHA, WAVELENGTH IN ANGSTROMS
!
!      wave = 2.997925d18 / min (freq, 2.463d15)
!
!.... CHANGE FROM DIVISION TO MULTIPLICATION
!
!      ww = 1.0d0 / (wave !! 2)
!      sig = (5.799d-13 + 1.422d-6 * ww + 2.784d0 * (ww * ww)) * 
!     &      (ww * ww)
!
!      do j = 1, nrhox
!         sigh(j) = sig * xnfph(j, 1) * 2.0d0 * bhyd(j, 1) / rho(j)
!      end do
!
!      end
!
!-------------------------------(above it is an old routine) ----
!
!   VW: implemented routine from Castelli, Jan. 2018:

      subroutine hrayop

!     gavrila, m. 1966, coherent scattering of light by atomic hydrogen,
!     jila report no. 86
      implicit none

!      include 'common.sizebl'
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
      double precision xsect,  sig, wave, ww, g 
      integer idum, map1, i, j
!
      logical cast_version 

      double precision gavrilam(74),gavrilamab(27),gavrilambc(24)
      double precision gavrilamcd(22), gavrilalymancont(64)
      double precision fgavrilalymancont(64)
     
      cast_version = .true.
 
      if (cast_version ) then 

!     fgavrilam .01 to .74 by .01
      data gavrilam/ -0.000113,  -0.000450,  -0.001014,  -0.001804,
     &  -0.002823,
     &  -0.004072,  -0.005553,  -0.007269,  -0.009223,  -0.011419,
     &  -0.013861,  -0.016553,  -0.019500,  -0.022709,  -0.026185,
     &  -0.029936,  -0.033968,  -0.038291,  -0.042913,  -0.047843,
     &  -0.053093,  -0.058674,  -0.064599,  -0.070882,  -0.077537,
     &  -0.084581,  -0.092031,  -0.099907,  -0.108230,  -0.117022,
     &  -0.126308,  -0.136117,  -0.146477,  -0.157422,  -0.168987,
     &  -0.181213,  -0.194143,  -0.207825,  -0.222313,  -0.237667,
     &  -0.253953,  -0.271245,  -0.289626,  -0.309189,  -0.330041,
     &  -0.352300,  -0.376103,  -0.401605,  -0.428985,  -0.458448,
     &  -0.490235,  -0.524625,  -0.561947,  -0.602591,  -0.647023,
     &  -0.695805,  -0.749619,  -0.809306,  -0.875910,  -0.950750,
     &  -1.035515,  -1.132403,  -1.244337,  -1.375285,  -1.530787,
     &  -1.718821,  -1.951320,  -2.246993,  -2.636960,  -3.177142,
     &  -3.979234,  -5.303624,  -7.930999, -15.763602/
!     fgavrilamab  .755 to .885 by .005
!
      data gavrilamab/ 31.008832, 15.382871, 10.160646,  7.538338,
     &  5.955062,
     &   4.890397,   4.121176,   3.535672,   3.071659,   2.691623,
     &   2.371483,   2.094936,   1.850395,   1.629203,   1.424526,
     &   1.230596,   1.042127,   0.853766,   0.659460,   0.451533,
     &   0.219115,  -0.054939,  -0.400868,  -0.879559,  -1.637857,
     &  -3.150374,  -8.326078/
!     fgavrilambc  .890 to .936 by .002

      data gavrilambc/  32.260389, 11.880702,7.418436, 5.442077,
     &   4.313409,
     &   3.573504,   3.043218,   2.637983,   2.312466,   2.039959,
     &   1.803441,   1.591244,   1.394717,   1.206823,   1.021148,
     &   0.831020,   0.628449,   0.402484,   0.136127,  -0.200462,
     &  -0.667435,  -1.410661,  -2.906862,  -8.169314/
!    5  -0.667435,  -1.410661,  -2.906862,  -8.169314, -26.332319/
!     fgavrilamcd  .938 to ,959 by .001
      data gavrilamcd/ 27.981406, 9.816495, 6.145775, 4.544224,
     &   3.630968,
     &   3.029081,   2.593248,   2.255265,   1.978565,   1.741426,
     &   1.529699,   1.333240,   1.143898,   0.954154,   0.755875,
     &   0.538760,   0.287687,  -0.022759,  -0.441666,  -1.081712,
     &  -2.278530,  -5.705843/
!
      data gavrilalymancont/ 2.667783, 2.526696, 2.408970, 2.308970,
     & 2.222736,
     &   2.147415,   2.080913,   2.021653,   1.968431,   1.920304,
     &   1.876527,   1.799739,   1.734455,   1.678180,   1.629118,
     &   1.585943,   1.547643,   1.513435,   1.482700,   1.454941,
     &   1.429751,   1.406798,   1.385804,   1.366536,   1.348797,
     &   1.332419,   1.317257,   1.303187,   1.290100,   1.277901,
     &   1.266509,   1.255848,   1.245856,   1.236474,   1.227652,
     &   1.219344,   1.190492,   1.167227,   1.148153,   1.132293,
     &   1.118945,   1.107593,   1.097848,   1.089413,   1.082059,
     &   1.075606,   1.069908,   1.064850,   1.060338,   1.056294,
     &   1.052655,   1.038936,   1.030042,   1.023928,   1.019536,
     &   1.016269,   1.011814,   1.008986,   1.007074,   1.005720,
     &   1.004724,   1.003970,   1.003385,   1.003140/
!
      data fgavrilalymancont/ 1.00,1.05,1.10,1.15,1.20,
     & 1.25,1.30,1.35,1.40,1.45,
     & 1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
     & 3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.4,4.8,5.2,5.6,
     & 6.0,6.4,6.8,7.2,7.6,8.0,8.4,8.8,9.2,9.6,10.,12.,14,16.,18.,
     & 20.,24.,28.,32.,36.,40.,44.,48.,50./

      xsect=0.

      if(freq .lt. 3.289841e13) then
        xsect=6.65246e-25*gavrilam(1)**2*(freq/3.289841e13)**4

      else if (freq .le. 2.4190611e15) then
        i=int(freq/3.289841e13) 
        i=min(i+1,74)
        g=gavrilam(i-1)+(gavrilam(i)-gavrilam(i-1))/3.289841e13*
     & (freq-(i-1)*3.289841e13)
        xsect=6.65246e-25*g**2

      else if (freq .lt. 0.77*3.289841e15) then 
!        nothing goes straight to end
      else if (freq .le. 0.885*3.289841e15) then
        i=int((freq-.755*3.289841e15)/1.644921e13) 
        i=i+1
        i=min(i+1,27)
        g=gavrilamab(i-1)+(gavrilamab(i)-gavrilamab(i-1))/1.644921e13*
     & (freq-(.755*3.289841e15+((i-1)-1)*1.6644921e13))
        xsect=6.65246e-25*g**2

      else if (freq .lt. 0.890*3.289841e15) then 
!         nothing goes straight to end 

      else if (freq .le. 0.936*3.289841e15) then
        i=int((freq-.890*3.289841e15)/0.657968e13) 
        i=i+1
        i=min(i+1,24)
        g=gavrilambc(i-1)+(gavrilambc(i)-gavrilambc(i-1))/0.657968e13*
     & (freq-(.890*3.289841e15+((i-1)-1)*0.657968e13))
        xsect=6.65246e-25*g**2

      else if (freq .lt. 0.938*3.289841e15) then 
!         goes to end go to 1

      else if (freq .le. 0.959*3.289841e15) then
 
        i=int((freq-.938*3.289841e15)/0.3289841e13) 
        i=i+1
        i=min(i+1,22)
        g=gavrilamcd(i-1)+(gavrilamcd(i)-gavrilamcd(i-1))/0.3289841e13*
     & (freq-(.938*3.289841e15+((i-1)-1)*0.3289841e13))
        xsect=6.65246e-25*g**2

      else if (freq .lt. 0.961*3.289841e15)  then 


      else 
        xsect=6.65246e-25*gavrilalymancont(1)

        if(freq .gt. 1.000*3.289841e15) then 

          idum=map1(fgavrilalymancont,gavrilalymancont,64,
     &          dble(freq/3.2898041d15),g,1)
          xsect=6.65246e-25*g
        end if 

      end if

      do j=1,nrhox
       sigh(j)=xsect*xnfph(j,1)*2.*bhyd(j,1)/rho(j)
      end do   


!---- cast_version end 
      else

!.... SET THE CUTOFF AT THE LYMAN ALPHA, WAVELENGTH IN ANGSTROMS
!
      wave = 2.997925d18 / min (freq, 2.463d15)
!
!.... CHANGE FROM DIVISION TO MULTIPLICATION
!
      ww = 1.0d0 / (wave ** 2)
      sig = (5.799d-13 + 1.422d-6 * ww + 2.784d0 * (ww * ww)) * 
     &      (ww * ww)

      do j = 1, nrhox
         sigh(j) = sig * xnfph(j, 1) * 2.0d0 * bhyd(j, 1) / rho(j)
      end do

      end if 

      end subroutine 
!---------------------------------------------------------------------
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H R A Y O P !!!!!!!!!!!!
!

#ifdef DFSY

      SUBROUTINE HE1OP

      IMPLICIT REAL*8 (A-H,O-Z)
      integer kw, j
      PARAMETER (kw=99)
!     REQUIRES FUNCTION COULFF
!      include 'common.sizebl'
      include 'common.constb'
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'




      FREQ3=2.815D29/FREQ/FREQ/FREQ
      DO 100 J=1,NRHOX
!     LIMIT 198310.76    RYDBERG 109722.267
!  N=6 TO INFINITY       BHE1=BHE2
      Htemp=FREQ3*4./2./(109722.267*HCKT(J))*
     1(EXP(-MAX(195262.919D0,198310.76D0-WAVENO)*HCKT(J))-
     2EXP(-198310.76*HCKT(J)))*STIM(J)*BHE2(J,1)
      S=Htemp*BNU(J)
!  5P 1P
      IF(WAVENO.LT.4368.190)GO TO 30
      X=FREQ3/3125.
      A=X* 3.*EXP(-193942.57 *HCKT(J))*(BHE1(J,29)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,29)/BHE2(J,1)-EHVKT(J))
!  5G 1G
      IF(WAVENO.LT.4388.260)GO TO 30
      X=FREQ3/3125.
      A=X* 9.*EXP(-193922.5  *HCKT(J))*(BHE1(J,28)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,28)/BHE2(J,1)-EHVKT(J))
!  5G 3G
      IF(WAVENO.LT.4388.260)GO TO 30
      X=FREQ3/3125.
      A=X*27.*EXP(-193922.5  *HCKT(J))*(BHE1(J,27)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,27)/BHE2(J,1)-EHVKT(J))
!  5F 1F
      IF(WAVENO.LT.4389.390)GO TO 30
      X=FREQ3/3125.
      A=X* 7.*EXP(-193921.37 *HCKT(J))*(BHE1(J,26)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,26)/BHE2(J,1)-EHVKT(J))
!  5F 3F
      IF(WAVENO.LT.4389.450)GO TO 30
      X=FREQ3/3125.
      A=X*15.*EXP(-193921.31 *HCKT(J))*(BHE1(J,25)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,25)/BHE2(J,1)-EHVKT(J))
!  5D 1D
      IF(WAVENO.LT.4392.369)GO TO 30
      X=FREQ3/3125.
      A=X* 5.*EXP(-193918.391*HCKT(J))*(BHE1(J,24)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,24)/BHE2(J,1)-EHVKT(J))
!  5D 3D
      IF(WAVENO.LT.4393.515)GO TO 30
      X=FREQ3/3125.
      A=X*15.*EXP(-193917.245*HCKT(J))*(BHE1(J,23)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,23)/BHE2(J,1)-EHVKT(J))
!  5P 3P
      IF(WAVENO.LT.4509.980)GO TO 30
      X=FREQ3/3125.
      A=X* 9.*EXP(-193800.78 *HCKT(J))*(BHE1(J,22)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,22)/BHE2(J,1)-EHVKT(J))

!  5S 1S
      IF(WAVENO.LT.4647.133)GO TO 30
      X=FREQ3/3125.
      A=X* 1.*EXP(-193663.627*HCKT(J))*(BHE1(J,21)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,21)/BHE2(J,1)-EHVKT(J))
!  5S 3S
      IF(WAVENO.LT.4963.671)GO TO 30
      X=FREQ3/3125.
      A=X* 3.*EXP(-193347.089*HCKT(J))*(BHE1(J,20)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,20)/BHE2(J,1)-EHVKT(J))
!  4P 1P
      IF(WAVENO.LT.6817.943)GO TO 30
      X=FREQ3/1024.
      A=X* 3.*EXP(-191492.817*HCKT(J))*(BHE1(J,19)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,19)/BHE2(J,1)-EHVKT(J))
!  4F 1F
      IF(WAVENO.LT.6858.680)GO TO 30
      X=FREQ3/1024.
      A=X* 7.*EXP(-191452.08 *HCKT(J))*(BHE1(J,18)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,18)/BHE2(J,1)-EHVKT(J))
!  4F 3F
      IF(WAVENO.LT.6858.960)GO TO 30
      X=FREQ3/1024.
      A=X*21.*EXP(-191451.80 *HCKT(J))*(BHE1(J,17)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,17)/BHE2(J,1)-EHVKT(J))
!  4D 1D
      IF(WAVENO.LT.6864.201)GO TO 30
      X=FREQ3/1024.
      A=X* 5.*EXP(-191446.559*HCKT(J))*(BHE1(J,16)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,16)/BHE2(J,1)-EHVKT(J))
!  4D 3D
      IF(WAVENO.LT.6866.172)GO TO 30
      X=FREQ3/1024.
      A=X*15.*EXP(-191444.588*HCKT(J))*(BHE1(J,15)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,15)/BHE2(J,1)-EHVKT(J))
!  4P 3P
      IF(WAVENO.LT.7093.620)GO TO 30
      X=FREQ3/1024.
      A=X* 9.*EXP(-191217.14 *HCKT(J))*(BHE1(J,14)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,14)/BHE2(J,1)-EHVKT(J))
!  4S 1S
      IF(WAVENO.LT.7370.429)GO TO 30
      X=FREQ3/1024.
      A=X* 1.*EXP(-190940.331*HCKT(J))*(BHE1(J,13)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,13)/BHE2(J,1)-EHVKT(J))
!  4S 3S
      IF(WAVENO.LT.8012.550)GO TO 30
      X=FREQ3/1024.
      A=X* 3.*EXP(-190298.210*HCKT(J))*(BHE1(J,12)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,12)/BHE2(J,1)-EHVKT(J))
!  3P 1P
      IF(WAVENO.LT.12101.289)GO TO 30
      X=EXP(58.81-2.89*FREQLN)
      A=X* 3.*EXP(-186209.471*HCKT(J))*(BHE1(J,11)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,11)/BHE2(J,1)-EHVKT(J))

!  3D 1D
      IF(WAVENO.LT.12205.695)GO TO 30
      X=EXP(85.20-3.69*FREQLN)
      A=X* 5.*EXP(-186105.065*HCKT(J))*(BHE1(J,10)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,10)/BHE2(J,1)-EHVKT(J))
!  3D 3D
      IF(WAVENO.LT.12209.106)GO TO 30
      X=EXP(85.20-3.69*FREQLN)
      A=X*15.*EXP(-186101.654*HCKT(J))*(BHE1(J, 9)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,9)/BHE2(J,1)-EHVKT(J))
!  3P 3P
      IF(WAVENO.LT.12746.066)GO TO 30
      X=EXP(49.30-2.60*FREQLN)
      A=X* 9.*EXP(-185564.694*HCKT(J))*(BHE1(J, 8)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,8)/BHE2(J,1)-EHVKT(J))
!  3S 1S
      IF(WAVENO.LT.13445.824)GO TO 30
      X=EXP(23.85-1.86*FREQLN)
      A=X* 1.*EXP(-184864.936*HCKT(J))*(BHE1(J, 7)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,7)/BHE2(J,1)-EHVKT(J))
!  3S 3S
      IF(WAVENO.LT.15073.868)GO TO 30
      X=EXP(12.69-1.54*FREQLN)
      A=X* 3.*EXP(-183236.892*HCKT(J))*(BHE1(J, 6)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,6)/BHE2(J,1)-EHVKT(J))
!  2P 1P
      IF(WAVENO.LT.27175.760)GO TO 30
      X=EXP(81.35-3.5*FREQLN)
      A=X* 3.*EXP(-171135.000*HCKT(J))*(BHE1(J, 5)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,5)/BHE2(J,1)-EHVKT(J))
!  2P 3P
      IF(WAVENO.LT.29223.753)GO TO 30
      X=EXP(61.21-2.9*FREQLN)
      A=X* 9.*EXP(-169087.007*HCKT(J))*(BHE1(J, 4)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,4)/BHE2(J,1)-EHVKT(J))
!  2S 1S
      IF(WAVENO.LT.32033.214)GO TO 30
      X=EXP(26.83-1.91*FREQLN)
      A=X* 1.*EXP(-166277.546*HCKT(J))*(BHE1(J, 3)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,3)/BHE2(J,1)-EHVKT(J))
!  2S 3S
      IF(WAVENO.LT.38454.691)GO TO 30
      X=EXP(-390.026+(21.035-.318*FREQLN)*FREQLN)
      A=X* 3.*EXP(-159856.069*HCKT(J))*(BHE1(J, 2)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,2)/BHE2(J,1)-EHVKT(J))
!  1S 1S
      IF(WAVENO.LT.198310.760)GO TO 30
      X=EXP(33.32-2.*FREQLN)
      A=X* 1.*1.                      *(BHE1(J, 1)-BHE2(J,1)*EHVKT(J))
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE1(J,1)/BHE2(J,1)-EHVKT(J))
   30 Htemp=Htemp*XNFPHE(J,1)/RHO(J)
      S=S*XNFPHE(J,1)/RHO(J)
!     FREE-FREE
      A=3.619D8/SQRT(T(J))*COULFF(J,1)/FREQ*XNE(J)/FREQ*
     1XNFHE(J,2)/FREQ*STIM(J)/RHO(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
      AHE1(J)=Htemp
      SHE1(J)=BNU(J)
      IF(Htemp.GT.0.)SHE1(J)=S/Htemp 
  100 CONTINUE
      RETURN
      END SUBROUTINE 



#else

!
      subroutine he1op

      implicit none
!
!.... 1995 OCT - TEST imin .le. 10 BEFORE hefreq(imin) TO AVOID ARRAY
!     OUT OF BOUNDS
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.constb'

!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  ahe1bound(maxd), ahe1free(maxd), bolt(maxd, 10),
     &                  boltn(maxd, 27), boltex(maxd), c, cfree, 
     &                  chi(10),  ex, exlim(maxd), 
     &                  freet(maxd), freqhen(10), freqi, g(10), he1, 
     &                  hefreq(11), 
     &                  rhoinv, rn2, trans(10), transn(27),  
     &                  xr, zeff2
      integer  imin, itemp1, j, n
      logical  first
!
!---------------------------- INITIALIZATION --------------------------
!
      data chi /  0.000d0, 19.819d0, 20.615d0, 20.964d0, 21.217d0, 
     &           22.718d0, 22.920d0, 23.006d0, 23.073d0, 23.086d0 /
!
      data first / .true. /
!
      data g / 1.0d0,  3.0d0, 1.0d0, 9.0d0, 3.0d0, 3.0d0, 1.0d0, 9.0d0, 
     &         20.0d0, 3.0d0 /
!
      data hefreq / 5.945209d15,  1.152844d15,  0.9603331d15, 
     &              0.8761076d15, 0.8147104d15, 0.4519048d15, 
     &              0.4030971d15, 0.3821191d15, 0.3660215d15,
     &              0.3627891d15, 0.0  /
!
      data itemp1 / 0 /
!
      save bolt, boltex, boltn, chi, exlim, first, freet, freqhen, g, 
     &     hefreq, itemp1
!
!------------------------------- EXECUTION ----------------------------
!
      if (first) then
         first = .false.
!
!.... IONIZATION FREQUENCY FROM HEI TO HEII N = 2
!
         freqhen(2)  = (he2lim - 159856.069d0) * c_cm
         freqhen(3)  = (he2lim - 166277.546d0) * c_cm
         freqhen(4)  = (he2lim - 169087.000d0) * c_cm
         freqhen(5)  = (he2lim - 171135.000d0) * c_cm
!
!.... IONIZATION FREQUENCY FROM HEI TO HEII N = 3
!
         freqhen(6)  = (he3lim - 183236.000d0) * c_cm
         freqhen(7)  = (he3lim - 184864.000d0) * c_cm
         freqhen(8)  = (he3lim - 185564.000d0) * c_cm
         freqhen(9)  = (he3lim - 186101.000d0) * c_cm
         freqhen(10) = (he3lim - 186209.471d0) * c_cm
      end if
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do j = 1, nrhox
            rhoinv = 1.0d0 / rho(j)
!
            do n = 1, 10
               bolt(j,n) = exp(-chi(n) / tkev(j)) * g(n) * 
     &                     xnfphe(j,1) * rhoinv
            end do
!
            do n = 4, 27
               rn2 = n * n
               boltn(j,n) = exp(-24.587d0 * (1.0d0 - 1.0d0 / rn2) / 
     &                      tkev(j)) * 4.0d0 * rn2 * xnfphe(j,1) * 
     &                      rhoinv
            end do
!
!.... THE FACTOR OF E10 HAS BEEN SWITCHED TO CFREE TO AVOID UNDERFLOW.
!
            freet(j)  = 1.0d-10 * xne(j) * rhoinv * xnfhe(j,2) / 
     &                  sqrt(t(j))
            xr        = xnfphe(j,1) * (4.0d0 / 2.0d0 / 13.595d0) * 
     &                  tkev(j) * rhoinv
            boltex(j) = exp(-23.730d0 / tkev(j)) * xr
            exlim(j)  = exp(-24.587d0 / tkev(j)) * xr
         end do
!
      end if
!
!.... FACTOR OF E10 IN CFREE IS FROM FREET
!
      freqi = 1.0d0 / freq
      cfree = ((1.0d10 * 3.6919d8 * freqi) * freqi) * freqi
      c     = ((2.815d29 * freqi) * freqi) * freqi
      imin  = 1
!
      do while (imin .le. 10 .and. hefreq(imin) .gt. freq)
         trans(imin) = 0.0d0
         imin = imin + 1
      end do
!
      if (imin .gt. 10) then
         imin = 0
!
      else
         if (imin .eq. 1 ) trans(1) = crosshe(freq)
         if (imin .le. 2 ) trans(2) = he12s3s(freq)
         if (imin .le. 3 ) trans(3) = he12s1s(freq)
         if (imin .le. 4 ) trans(4) = he12p3p(freq)
         if (imin .le. 5 ) trans(5) = he12p1p(freq)
!
!.... 1S3S 3S
!
         if (imin .le. 6 ) trans(6) = xkarsas(1.236439d0, 3, 0)
!
!.... 1S3S 1S
!
         if (imin .le. 7 ) trans(7) = xkarsas(1.102898d0, 3, 0)
!
!.... 1S3P 3P
!
         if (imin .le. 8 ) trans(8) = xkarsas(1.045499d0, 3, 1)
!
!.... 1S3D 3D+1D
!
         if (imin .le. 9 ) trans(9) = xkarsas(1.001427d0, 3, 2)
!
!.... 1S3P 1P.  NO TEST NEEDED HERE
!
         trans(10) = xkarsas(0.9926d0, 3, 1)
!
!..... TO HEII N=2
!
         n = 5
!
         do while (freq .ge. freqhen(n) .and. n .ge. 2)
            zeff2    = freqhen(n) / ryd
            trans(n) = trans(n) + xkarsas(zeff2, 1, 0)
            n = n - 1
         end do
!
!.... TO HEII N = 3
!
         n = 10
!
         do while (freq .ge. freqhen(n) .and. n .ge. 6)
            zeff2    = freqhen(n) / ryd
            trans(n) = trans(n) + xkarsas(zeff2, 1, 0)
            n = n - 1
         end do
!
         if(freq .ge. 1.25408d+16) then
!
            do n = 4, 27
               rn2 = n * n
               zeff2 = 4.0d0 - 3.0d0 / rn2
               transn(n) = xkarsas(zeff2, 1, 0)
            end do
!
         end if
!
      end if
!
      do j = 1, nrhox
         ex = boltex(j)
         if (freq .lt. 2.055d14) ex = exlim(j) / ehvkt(j)
         he1 = (ex - exlim(j)) * c
!
         if (imin .gt. 0) then
!
            do n = imin, 10
               he1 = he1 + trans(n) * bolt(j,n)
            end do
!
            if (freq .ge. 1.25408d+16) then
!
               do n = 4, 27
                  he1 = he1 + transn(n) * boltn(j,n)
               end do
!
            end if
!
         end if
!
         ahe1bound(j) = he1 * stim(j)
         ahe1free(j)  = (coulff(j,1) * freet(j) * cfree) * stim(j)
         ahe1(j) = (he1 + coulff(j,1) * freet(j) * cfree) * stim(j)
      end do
!
      end subroutine 

#endif

!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H E 1 O P !!!!!!!!!!!!!!
!
      double precision function crosshe(freq)

      implicit none
!.... MARR, G.V. AND WEST, J.B. ATOMIC DATA AND NUCLEAR DATA TABLES,
!.... VOL 18, 497-508, 1976.
!
!------------------------------ DUMMY VARIABLE -------------------------
!
      double precision freq
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!
!------------------------------ LOCAL VARIABLES ------------------------
!
      double precision  wave, x505(92), x50(16), x20(11), x10(21)
      integer i
!
!------------------------------ LOCAL VARIABLES ------------------------
!
      data x505 /
     &     7.58,   7.46,   7.33,   7.19,  7.06,  6.94,  6.81,  6.68,
     &     6.55,   6.43,   6.30,   6.18,  6.05,  5.93,  5.81,  5.69,
     &     5.57,   5.45,   5.33,   5.21,  5.10,  4.98,  4.87,  4.76,
     &     4.64,   4.53,   4.42,   4.31,  4.20,  4.09,  4.00,  3.88,
     &     3.78,   3.68,   3.57,   3.47,  3.37,  3.27,  3.18,  3.08,
     &     2.98,   2.89,   2.80,   2.70,  2.61,  2.52,  2.44,  2.35,
     &     2.26,   2.18,   2.10,   2.02,  1.94,  1.86,  1.78,  1.70,
     &     1.63,   1.55,   1.48,   1.41,  1.34,  1.28,  1.21,  1.14,
     &     1.08,   1.02,   0.961,  0.903, 0.847, 0.792, 0.738, 0.687,
     &     0.637,  0.588,  0.542,  0.497, 0.454, 0.412, 0.373, 0.335,
     &     0.299,  0.265,  0.233,  0.202, 0.174, 0.147, 0.123, 0.100,
     &     0.0795, 0.0609, 0.0443, 0.0315/
!
      data x50 /
     &     0.0315,   0.0282,   0.0250,   0.0220,   0.0193,
     &     0.0168,   0.0145,   0.0124,   0.0105,   0.00885,
     &     0.00736,  0.00604,  0.00489,  0.00389,  0.00303,
     &     0.00231/
!
      data x20 /
     &   0.00231,  0.00199,  0.00171,  0.00145,  0.00122, 0.00101,
     &   0.000832, 0.000673, 0.000535, 0.000417, 0.000318/
!
      data x10/
     &   0.000318, 0.000274, 0.000235, 0.000200,  0.000168,
     &   0.000139, 0.000115, 0.000093, 0.000074,  0.000057,
     &   0.000044, 0.000032, 0.000023, 0.000016,  0.000010,
     &   0.000006, 0.000003, 0.000001, 0.0000006, 0.0000003,
     &   0.000000 /
!
!------------------------------ EXECUTION ------------------------------
!
      if (freq .ge. 5.945209d15) then
         wave = 2.99792458d18 / freq
!
         if (wave .gt. 50.0d0) then
            i = 93.0 - (wave - 50.0d0) /  5.0d0 
            i = min(92, max(2, i))
            crosshe = ((wave - 50.0d0 - dble((92 - i) * 5) ) /  5.0d0  * 
     &                (x505(i-1) - x505(i)) + x505(i)) * 1.0d-18
!
         else if (wave .gt. 20.0d0) then
            i = 17.0 - (wave - 20.0d0) / 2.0
            i = min(16, max(2, i))
            crosshe = ((wave - 20.0d0 - dble((16 - i) * 2) ) / 2.0 *
     &                (x50(i-1) - x50(i)) + x50(i)) * 1.0d-18
!
         else if (wave .gt. 10.0d0) then
            i = 12.0 - (wave - 10.0d0) / 1.0
            i = min(11, max(2, i))
            crosshe = ((wave - 10.0d0 - dble((11 - i) * 1) ) / 1.0 * 
     &                (x20(i-1) - x20(i)) + x20(i)) * 1.0d-18
!
         else
            i = 22.0 - wave / 0.5
            i = min(21, max(2, i))
            crosshe = ((wave - dble(21 - i) * 0.5) / 0.5 * 
     &                (x10(i-1) - x10(i)) + x10(i)) * 1.0d-18
         end if
!
      else
         crosshe = 0.0d0
      end if
!
      end function
!
*!!!!!!!!!! E N D   O F   S U B R O U T I N E   C R O S S H E !!!!!!!!!!
!
      double precision function he111s(freq)
      implicit none
!
!.... FOLLOWING MATHISEN
!
!-------------------------------- DUMMY VARIABLE -----------------------
!
      double precision freq
!
!-------------------------------- LOCAL VARIABLE -----------------------
!
      double precision  w(64), wave, x(64)
      integer i
!
!-------------------------------- INITIALIZATION -----------------------
!
      data w /
     &   504.3,  501.5, 498.7, 493.3, 488.1, 480.3, 477.8, 454.0, 
     &   443.0,  395.0, 356.4, 348.2, 324.6, 302.0, 298.1, 275.6, 
     &   260.6,  256.2, 239.4, 224.6, 220.0, 215.0, 210.0, 205.0,  
     &   200.0,  195.0, 190.0, 185.0, 180.0, 175.0, 170.0, 165.0,  
     &   160.0,  155.0, 150.0, 145.0, 135.0, 130.0, 125.0, 120.0,  
     &   115.0,  110.0, 105.0, 100.0,  95.0,  90.0,  85.0,  80.0,   
     &    75.0,  70.0,   65.0,  60.0,  55.0,  50.0,  45.0,  40.0,   
     &    35.0,  30.0,   25.0,  20.0,  15.0,  10.0,   5.0,   0.0 /
!
      data x /
     &   7.346,   7.317,   7.259,    7.143,   7.030,   
     &   6.857,   6.800,   6.284,    6.041,   4.977,  
     &   4.138,   3.961,   3.474,    3.025,   2.945, 
     &   2.522,   2.259,   2.179,    1.901,   1.684,   
     &   1.61,    1.53,    1.45,     1.38,    1.30,  
     &   1.22,    1.14,    1.08,     1.02,    0.961,   
     &   0.903,   0.847,   0.792,    0.738,   0.687,   
     &   0.637,   0.542,   0.497,    0.454,   0.412, 
     &   0.373,   0.335,   0.299,    0.265,   0.233,   
     &   0.202,   0.174,   0.147,    0.124,   0.103,  
     &   0.0840,  0.0676,  0.0535,   0.0414,  0.0311,
     &   0.0266,  0.0158,  0.0104,   0.00637, 0.00349, 
     &   0.00161, 0.00054, 0.000083, 0.0 /
!
!------------------------------- EXECUTON -----------------------------
!
      if (freq .ge. 5.945209d15) then
         wave = 2.99792458d18 / freq
         i = 2
!
         do while (w(i) .ge. wave .and. i .lt. 64)
            i = i + 1
         end do
!
         he111s = ((wave - w(i)) / (w(i-1) - w(i)) * 
     &             (x(i-1) - x(i)) + x(i)) * 1.0d-18
!
      else
         he111s = 0.0d0
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   H E 1 1 1 S !!!!!!!!!!!!
!
      double precision function he12s1s(f)

      implicit none
!
!.... 1995 JUL - ADDING common.freqbl MADE THE VARIABLE freqlg VISIBLE 
!                HERE, WHERE IT WAS REDEFINED TO BE LOG10.  TO SOLVE THIS
!                CORRUPTION, CHANGE THE DEFINITION OF freqlg = log10(freq),
!                AND CREATE freqln = log(freq)
!
!---------------------------- DUMMY VARIABLE ---------------------------
!
      double precision f
!
!------------------------------- COMMON --------------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!      include 'common.freqbl'
!
!------------------------------- CONSTANT ------------------------------
!
!.... FRQLIM1 = 32033.214 * 2.99792458D10
!.... FRQLIM2 = 2.4 * 109722.267 * 2.99792458D10
!
      double precision frqlim1, frqlim2
      parameter ( frqlim1 = 9.6033159627d14,
     &            frqlim2 = 7.89453794911d15)
!
!--------------------------- LOCAL VARIABLES ---------------------------
!
      double precision  ek, eps, freq1s(16), x, x1s(16)
      integer i
!
!---------------------------- INITIALIZATION ---------------------------
!
      data freq1s /
     &    15.947182,   15.913654,   15.877320,   15.837666,   15.794025,
     &    15.745503,   15.690869,   15.628361,   15.555317,   15.467455,
     &    15.357189,   15.289399,   15.251073,   15.209035,   15.162487,
     &    14.982421/
!
      data x1s/
     &   -19.635557,  -19.159345,  -18.958474,  -18.809535,  -18.676481,
     &   -18.546006,  -18.410962,  -18.264821,  -18.100205,  -17.909165,
     &   -17.684370,  -17.557867,  -17.490360,  -17.417876,  -17.349386,
     &   -17.084441/
!
!------------------------------ EXECUTION  ----------------------------
!
      if (freq .ge. frqlim1) then
!
         if(freq .gt. frqlim2) then
            ek = (waveno - 32033.214d0) / 109722.267d0
            eps = 2.0 * (ek - 2.612316d0) / 0.00322d0
            he12s1s = 0.008175d0 * (484940.0d0 / waveno) ** 2.71d0 * 
     &                8.067d-18 * (eps + 76.21d0) ** 2 / 
     &               (1.0d0 + eps**2)
!
         else
            i = 2
!
            do while (freq1s(i) .ge. freqlg .and. i .lt. 16)
               i = i + 1
            end do 
!
            x = (freqlg - freq1s(i)) / (freq1s(i-1) - freq1s(i)) *
     &          (x1s(i-1) - x1s(i)) + x1s(i)
            he12s1s = exp(x * tenlog) 
!
         end if
!
      else
         he12s1s = 0.0d0
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   H E 1 2 S 1 S !!!!!!!!!!
!
      double precision function he12s3s(f)

      implicit none
!
!.... 1995 JUL - ADDING common.freqbl MADE THE VARIABLE freqlg VISIBLE 
!                HERE, WHERE IT WAS REDEFINED TO BE LOG10.  TO SOLVE THIS
!                CORRUPTION, CHANGE THE DEFINITION OF freqlg = log10(freq),
!                AND CREATE freqln = log(freq)
!
!---------------------------- DUMMY VARIABLE ---------------------------
!
      double precision f
!
!------------------------------- COMMON --------------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!      include 'common.freqbl'
!
!------------------------------- CONSTANT ------------------------------
!
!.... FRQLIM1 = 38454.691 * 2.99792458D10
!.... FRQLIM2 = 2.4 * 109722.267 * 2.99792458D10
!
      double precision  frqlim1, frqlim2
      parameter ( frqlim1 = 1.15284263365d15,
     &            frqlim2 = 7.89453794911d15)
!
!--------------------------- LOCAL VARIABLES ---------------------------
!
      double precision  ek, eps, freq3s(16), x, x3s(16)
      integer i
!
!---------------------------- INITIALIZATION ---------------------------
!
      data freq3s /
     &    15.956523,   15.923736,   15.888271,   15.849649,   15.807255,
     &    15.760271,   15.707580,   15.647601,   15.577992,   15.495055,
     &    15.392451,   15.330345,   15.295609,   15.257851,   15.216496,
     &    15.061770/
!
      data x3s /
     &   -18.426022,  -18.610700,  -18.593051,  -18.543304,  -18.465513,
     &   -18.378707,  -18.278574,  -18.164329,  -18.033346,  -17.882435,
     &   -17.705542,  -17.605584,  -17.553459,  -17.500667,  -17.451318,
     &   -17.266686/
!
!------------------------------ EXECUTION  ----------------------------
!
      if (freq .ge. frqlim1) then
!
         if(freq .gt. frqlim2) then
            ek = (waveno - 38454.691d0) / 109722.267d0
            eps = 2.0 * (ek - 2.47898d0) / 0.000780d0
            he12s3s = 0.01521d0 * (470310.0d0 / waveno) ** 3.12d0 * 
     &                8.067d-18 * (eps - 122.4d0) ** 2 / 
     &                (1.0d0 + eps ** 2)
!
         else
            i = 2
!
            do while (freq3s(i) .ge. freqlg .and. i .lt. 16)
               i = i + 1
            end do
!
            x = (freqlg - freq3s(i)) / (freq3s(i-1) - freq3s(i)) *
     &          (x3s(i-1) - x3s(i)) + x3s(i)
            he12s3s = exp(x * tenlog)
         end if
!
      else
         he12s3s = 0.0d0
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   H E 1 2 S 3 1 !!!!!!!!!!
!
      double precision function he12p1p(f)
      implicit none
!
!.... 1995 JUL - ADDING common.freqbl MADE THE VARIABLE freqlg VISIBLE 
!                HERE, WHERE IT WAS REDEFINED TO BE LOG10.  TO SOLVE THIS
!                CORRUPTION, CHANGE THE DEFINITION OF freqlg = log10(freq),
!                AND CREATE freqln = log(freq)
!
!---------------------------- DUMMY VARIABLE ---------------------------
!
      double precision f
!
!------------------------------- COMMON --------------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!      include 'common.freqbl'
!
!------------------------------- CONSTANT ------------------------------
!
!.... FRQLIM1 = 27175.76 * 2.99792458D10
!.... FRQLIM2 = 2.4 * 109722.267 * 2.99792458D10
!
      double precision  frqlim1, frqlim2
      parameter ( frqlim1 = 8.14708788842d14,
     &            frqlim2 = 7.89453794911d15)
!
!--------------------------- LOCAL VARIABLES ---------------------------
!
      double precision  ek, eps1d, eps1s, freq1p(16), x,
     &                  x1p(16)
      integer i
!
!---------------------------- INITIALIZATION ---------------------------
!
      data freq1p /
     &    15.939981,   15.905870,   15.868850,   15.828377,   15.783742,
     &    15.733988,   15.677787,   15.613218,   15.537343,   15.445346,
     &    15.328474,   15.255641,   15.214064,   15.168081,   15.116647,
     &    14.911002/
!
      data x1p /
     &   -18.798876,  -19.685922,  -20.011664,  -20.143030,  -20.091354,
     &   -19.908333,  -19.656788,  -19.367745,  -19.043016,  -18.674484,
     &   -18.240861,  -17.989700,  -17.852015,  -17.702677,  -17.525347,
     &   -16.816344/
!
!------------------------------ EXECUTION  ----------------------------
!
      if (freq .ge. frqlim1) then
!
         if(freq .gt. frqlim2) then
            ek = (waveno - 27175.76d0) / 109722.267d0
            eps1s = 2.0 * (ek - 2.446534d0) / 0.01037d0
            eps1d = 2.0 * (ek - 2.59427d0)  / 0.00538d0
            he12p1p = 0.0009487d0 * (466750.0d0 / waveno) ** 3.69d0 * 
     &                8.067d-18 * ((eps1s - 29.30d0) ** 2 / 
     &                (1.0d0 + eps1s ** 2) + (eps1d + 172.4d0) ** 2 / 
     &                (1.0d0 + eps1d ** 2))
!
         else
            i = 2
!
            do while (freq1p(i) .ge. freqlg .and. i .lt. 16)
               i = i + 1
            end do
!
            x = (freqlg - freq1p(i)) / (freq1p(i-1) - freq1p(i)) *
     &          (x1p(i-1) - x1p(i)) + x1p(i)
            he12p1p = exp(x * tenlog)
         end if
!
      else
         he12p1p = 0.0d0
      end if
!
      end function 
!
!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   H E 1 2 P 1 P !!!!!!!!!!
!
      double precision function he12p3p(f)
      implicit none
!
!.... 1995 JUL - ADDING common.freqbl MADE THE VARIABLE freqlg VISIBLE 
!                HERE, WHERE IT WAS REDEFINED TO BE LOG10.  TO SOLVE THIS
!                CORRUPTION, CHANGE THE DEFINITION OF freqlg = log10(freq),
!                AND CREATE freqln = log(freq)
!
!---------------------------- DUMMY VARIABLE ---------------------------
!
      double precision f
!
!------------------------------- COMMON --------------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!      include 'common.freqbl'
!
!------------------------------- CONSTANT ------------------------------
!
!.... FRQLIM = 29223.753 * 2.99792458D10
!
      double precision frqlim
      parameter ( frqlim = 8.76106074385d14)
!
!--------------------------- LOCAL VARIABLES ---------------------------
!
      double precision  freq3p(16), x, x3p(16)
      integer i
!
!---------------------------- INITIALIZATION ---------------------------
!
      data freq3p /
     &    15.943031,   15.909169,   15.872441,   15.832318,   15.788107,
     &    15.738880,   15.683351,   15.619667,   15.545012,   15.454805,
     &    15.340813,   15.270195,   15.230054,   15.185821,   15.136567,
     &    14.942557/
!
      data x3p /
     &   -19.791021,  -19.697886,  -19.591421,  -19.471855,  -19.337053,
     &   -19.183958,  -19.009750,  -18.807990,  -18.570571,  -18.288361,
     &   -17.943476,  -17.738737,  -17.624154,  -17.497163,  -17.403183,
     &   -17.032999/
!
!------------------------------ EXECUTION  ----------------------------
!
      if (freq .ge. frqlim) then
         i = 2
!
         do while (freq3p(i) .ge. freqlg .and. i .lt. 16)
            i = i + 1
         end do
!
         x = (freqlg - freq3p(i)) / (freq3p(i-1) - freq3p(i)) *
     &       (x3p(i-1) - x3p(i)) + x3p(i)
         he12p3p = exp(x * tenlog)
!
      else
         he12p3p = 0.0d0
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   H E 1 2 P 3 P !!!!!!!!!!
!
#ifdef DFSY
      SUBROUTINE HE2OP
      IMPLICIT REAL*8 (A-H,O-Z)
      integer kw, j
      PARAMETER (kw=99)
!     REQUIRES FUNCTION COULFF
!      include 'common.sizebl'
      include 'common.constb'
!
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'




      FREQ3=2.815D29/FREQ/FREQ/FREQ
      DO 100 J=1,NRHOX
!     LIMIT 438908.85    4*RYDBERG 438889.068
!  N=10 TO INFINITY
      XNFPRHO=XNFPHE(J,2)/RHO(J)
      Htemp=FREQ3*16. *2./2./(438889.068*HCKT(J))*
     1(EXP(-MAX(434519.959D0,438908.85D0-WAVENO)*HCKT(J))-
     2EXP(-438908.85*HCKT(J)))*STIM(J)*XNFPRHO
      S=Htemp*BNU(J)
!  N=9
      IF(WAVENO.LT.5418.390)GO TO 30
      X=FREQ3/59049.*16.
      A=X*162.*EXP(-433490.46*HCKT(J))*STIM(J)*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=8
      IF(WAVENO.LT.6857.660)GO TO 30
      X=FREQ3*16./32768.
      A=X*128.*EXP(-432051.19*HCKT(J))*STIM(J)*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=7
      IF(WAVENO.LT.8956.950)GO TO 30
      X=FREQ3*16./16807.
      A=X*98.*EXP(-429951.90*HCKT(J))*STIM(J)*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)
!  N=6
      IF(WAVENO.LT.12191.437)GO TO 30
      X=FREQ3*16./7776.*(1.0986+(-2.704D13+1.229D27/FREQ)/FREQ)
      A=X*72.*EXP(-426717.413*HCKT(J))*(BHE2(J,6)-EHVKT(J))*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE2(J,6)-EHVKT(J))
!  N=5
      IF(WAVENO.LT.17555.715)GO TO 30
      X=FREQ3*16./3125.*(1.102+(-3.909D13+2.371D27/FREQ)/FREQ)
      A=X*50.*EXP(-421353.135*HCKT(J))*(BHE2(J,5)-EHVKT(J))*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE2(J,5)-EHVKT(J))
!  N=4
      IF(WAVENO.LT.27430.925)GO TO 30
      X=FREQ3*16./1024.*(1.101+(-5.765D13+4.593D27/FREQ)/FREQ)
      A=X*32.*EXP(-411477.925*HCKT(J))*(BHE2(J,4)-EHVKT(J))*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE2(J,4)-EHVKT(J))
!  N=3
      IF(WAVENO.LT.48766.491)GO TO 30
      X=FREQ3*16./243.*(1.101+(-9.863D13+1.035D28/FREQ)/FREQ)
      A=X*18.*EXP(-390142.359*HCKT(J))*(BHE2(J,3)-EHVKT(J))*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE2(J,3)-EHVKT(J))
!  N=2
      IF(WAVENO.LT.109726.529)GO TO 30
      X=FREQ3*16./32.*(1.105+(-2.375D14+4.077D28/FREQ)/FREQ)
      A=X* 8.*EXP(-329182.321*HCKT(J))*(BHE2(J,2)-EHVKT(J))*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE2(J,2)-EHVKT(J))

!  N=1
      IF(WAVENO.LT.438908.850)GO TO 30
      X=FREQ3*16./1.*(.9916+(2.719D13-2.268D30/FREQ)/FREQ)
      A=X* 2.*1.                      *(BHE2(J,1)-EHVKT(J))*XNFPRHO
      Htemp=Htemp+A
      S=S+A*BNU(J)*STIM(J)/(BHE2(J,1)-EHVKT(J))
!
!     FREE-FREE
!   30 A=3.6919E8/SQRT(T(J))*COULFF(J,2)/FREQ*XNE(J)/FREQ*
   30 A=3.6919D8*4./SQRT(T(J))*COULFF(J,2)/FREQ*XNE(J)/FREQ*
     1XNFPHE(J,3)/FREQ*STIM(J)/RHO(J)
      Htemp=Htemp+A
      S=S+A*BNU(J)
      AHE2(J)=Htemp
      SHE2(J)=BNU(J)
      IF(Htemp.GT.0.)SHE2(J)=S/Htemp
  100 CONTINUE
! 100 SHE2(J)=S/H
      RETURN
      END SUBROUTINE 




#else
      subroutine he2op
      implicit none
!
!.... REQUIRES FUNCTIONS coulff AND xkarsas
!.... FREQUENCIES ARE 4X HYDROGEN, CHI ARE FOR ION POT = 54.403
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  ahe2bound(maxd), ahe2free(maxd), bolt(maxd, 9), 
     &                  boltex(maxd), c, cfree, cont(9),  ex, 
     &                  exlim(maxd), freet(maxd), freqi, he2, rhoinv, 
     &                  rn2, xr
      integer itemp1, j, n
      save bolt, boltex, exlim, freet, xr
!
!------------------------------ INITIALIZATION ------------------------
!
      save itemp1
      data itemp1 / 0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do j = 1, nrhox
            rhoinv = 1.0 / rho(j)
!
            do n = 1, 9
               rn2 = n * n
               bolt(j, n) = exp(-(54.403 - 54.403 / rn2) / tkev(j)) * 
     &                      2.0 * rn2 * xnfphe(j, 2) * rhoinv
            end do
!
!.... TAKE OUT 1.0E10 HERE AND APPLY TO CFREE
            freet(j) = 1.0d-10* xne(j) * xnfhe(j, 3) /
     &                 (rho(j) * sqrt(t(j)))
            xr = xnfphe(j, 2) * (2.0 / 2.0 / 13.595) * tkev(j) * rhoinv
            boltex(j) = exp(-53.859d0 / tkev(j)) * xr
            exlim(j)  = exp(-54.403d0 / tkev(j)) * xr
         end do
!
      end if
!
      do n = 1, 9
         cont(n) = xkarsas(4.0d0, n, n)
      end do
!
      freqi = 1.0d0 / freq
!
!.... HAVE INCREASED CONSTANT BY 1.0E10 FROM FREET
!
      cfree = freqi * (freqi * (freqi * 3.6919d18) ) * 4.0
      c = freqi * (freqi * (freqi * 2.815d29) ) * 4.0
!
      do j = 1, nrhox
         ex = boltex(j)
         if (freq .lt. 1.31522d14) ex = exlim(j) / ehvkt(j)
         he2 = (ex - exlim(j)) * c
!
         do n = 1, 9
            he2 = he2 + cont(n) * bolt(j, n)
         end do
!
         ahe2bound(j) = he2 * stim(j)
         ahe2free(j) = (coulff(j, 2) * cfree * freet(j)) * stim(j)
         ahe2(j) = (he2 + coulff(j, 2) * cfree * freet(j)) * stim(j)
      end do
!
      end SUBROUTINE 

#endif

!*
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H E 2 O P !!!!!!!!!!!!!!
!
      subroutine hemiop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision a, b, c
      integer j
!
!------------------------------- EXECUTION ----------------------------
!
!.... REDUCE A, B AND C BY 1.0D-30 TO AVOID VAX OVERFLOW. PUT IT BACK IN
!.... THE CALCULATION OF AHEMIN
!
      a = 3.397d-16 + (-5.216d-1 + 7.039d15 / freq) / freq
      b = -4.116d-12 + (1.067d4 + 8.135d19 / freq) / freq
      c = 5.081d-7 + (-8.724d7 - 5.659d22 / freq) / freq
!
      do j = 1, nrhox
         ahemin(j) = (a * t(j) + b + c / t(j) ) * xne(j) * 
     &               xnfphe(j, 1) / rho(j) * 1.0d-30
      end do
!
      end subroutine
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H E M I O P !!!!!!!!!!!!
!
      subroutine heraop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
!      !include 'common.rhoxbl'
      include 'common.stateb'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  sig, wave, ww
      integer j
!
!------------------------------- EXECUTION ----------------------------
!
      wave = 2.997925d18 / min (freq, 5.15d15)
      ww = wave ** 2
      sig = 5.484d-14 / ww / ww * (1.0 + (2.44d5 + 5.94d10 / 
     &      (ww - 2.90d5) ) / ww) ** 2
!
      do j = 1, nrhox
         sighe(j) = sig * xnfphe(j, 1) / rho(j)
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H E R A O P !!!!!!!!!!!!
!
      subroutine coolop
      implicit none
!
!.... SI1, MG1, AL1, C1, FE1, OH, CH
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.freqbl'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xnfpbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
!      double precision  al1op, c1op, chop, fe1op, mg1op, ohop, si1op
      integer j
      logical cast_version  
!
!-------------------------------- EXTERNALS ---------------------------
!
!      external al1op, c1op, chop, fe1op, mg1op, ohop, si1op
!
!------------------------------- EXECUTION ----------------------------
!
      cast_version = .true.  

      if (cast_version ) then 
        call si1opcastelli
 
      do j = 1, nrhox
         acool(j) = (c1op(j)  * xnfpc(j,1)  +
     &               al1op(j) * xnfpal(j,1) +
     &               fe1op(j) * xnfpfe(j,1) +
     &               mg1op(j) * xnfpmg(j,1) +
     &               chop(j)  * xnfpch(j)   +
     &               ohop(j)  * xnfpoh(j))  * stim(j) / rho(j) +    
     &               asi1(j)  * xnfpsi(j,1) / rho(j)

      end do

    
      else


      do j = 1, nrhox
         acool(j) = (c1op(j)  * xnfpc(j,1)  + 
     &               al1op(j) * xnfpal(j,1) + 
     &               fe1op(j) * xnfpfe(j,1) +
     &               mg1op(j) * xnfpmg(j,1) +
     &               si1op(j) * xnfpsi(j,1) +
     &               chop(j)  * xnfpch(j)   + 
     &               ohop(j)  * xnfpoh(j))  * stim(j) / rho(j)
      end do
     
      end if
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   C O O L O P !!!!!!!!!!!!
!
      double precision function al1op(j)
      implicit none
!
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
!
!--------------------------------- DUMMY VARIABLE ----------------------
!.... NOT USED BUT NEEDED TO KEEP THE COMPILER HAPPY
!
      integer j
!
!------------------------------- COMMON --------------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
!
!------------------------------- EXECUTION -----------------------------
!
      if (freq .ge. 1.443d15) then
!         al1op = 6.5d-17 * (1.443d15 / freq) ** 5 * 6.0 ! Lester version !
          al1op = 6.5d-17 * (1.443d15 / freq) ** 4.4 * 6.0
!
!
      else
         al1op = 0.0d0
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   A L 1 O P !!!!!!!!!!!!!!
!
      double precision function c1op (j)
      implicit none
!
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
!  needed for Castelli version
      include 'common.constb'
 
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  c1240(maxd), c1444(maxd), freq1,  x1100,
     &                  x1240, x1444
      integer itemp1, l
!
!------------------------------ INITIALIZATION ------------------------
!
      save c1240, c1444, freq1, itemp1, x1100, x1240, x1444
      data freq1, itemp1 / 0.0d0, 0 /
!----- VW change to current Castelli version : 
      logical cast_version
      double precision xs0, xs1, xd0, xd1, xd2, eps
!
!------------------------------- EXECUTION ----------------------------
!
      cast_version = .true.

      
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do l = 1, nrhox
            c1240(l) = 5.0 * exp(-1.264 / tkev(l))
            c1444(l) = exp(-2.683 / tkev(l))
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1 = freq
         x1100 = 0.0d0
         x1240 = 0.0d0
         x1444 = 0.0d0
        
        if ( cast_version ) then

          xs0 = 0.0d0
          xs1 = 0.0d0
          xd0 = 0.0d0
          xd1 = 0.0d0
          xd2 = 0.0d0

         if (freq .ge. 2.7254d15) x1100 = 10.d0**(-16.80-(waveno-
     &                                    90777.0d0)/3.0d0/ryd/c_cm)

         if (freq .ge. 2.4196e15) then 
!     Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.!
          xd0 = 10.0d0**(-16.8d0-(waveno - 80627.760)/3.0d0/ryd/c_cm)
!     Fit to Burke, P.G. and Taylor, K.T. 1979, J. Phys. B, 12, 2971-2984.          .
          eps = (waveno - 93917.0d0)*2.0d0/2743.0d0
          xd1 = ((22.0d-18)*eps+(26.0d-18))/(eps**2+1.0d0)
!
          eps = (waveno - 111130.0d0)*2.0d0/2743.0d0  
          xd2 = ((-10.5d-18)*eps +(46.0d-18))/(eps**2 + 1.0d0)
! 
          x1240 = xd0 +xd1 + xd2
         end if

         if (freq .ge. 2.0761d15) then 
!      Luo, D. and Pradhan, A.K. 1989, J.Phys. B, 22, 3377-3395.
          xs0 = 10.0d0**(-16.80-(waveno -69172.400)/3.d0/ryd/c_cm )
          eps = (waveno - 97700.0d0)*2.0d0/2743.0d0
!     Fit to Burke, P.G. and Taylor, K.T. 1979, J. Phys. B, 12, 2971-2984.
          xs1 = (68.0d-18*eps + 118.0d-18) /(eps**2 +1.d0) 
          x1444=xs0 + xs1 
         end if
 

        else
!       here the lester version : 
         if (freq .ge. 2.7254d15) x1100 = seaton (2.7254d15, 1.219d-17, 
     &                                            2.0d0, 3.317d0)

         if (freq .ge. 2.4196d15) x1240 = seaton (2.4196d15, 1.03d-17,  
     &                                            1.5d0, 2.789d0)
         if (freq .ge. 2.0761d15) x1444 = seaton (2.0761d15, 9.59d-18,  
     &                                            1.5d0, 3.501d0)

        end if ! ---> close if for different versions

!
      end if
!
      c1op = x1100 * 9.0 + x1240 * c1240(j) + x1444 * c1444(j)
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   C 1 O P !!!!!!!!!!!!!!!!
!
      double precision function chop(j)
      implicit none
!
!.... CROSS-SECTION TIMES PARTITION FUNCTION
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.freqbl'
      include 'common.tempbl'
      include 'common.stateb'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  crossch(15, 105), crosscht(15), en, evolt, 
     &                  freq1, part, partch(41), tn
      integer it, n
      save    crossch, crosscht, freq1, n, partch
!
!----------------------------- INITIALIZATION --------------------------
!
      data (crossch(it, 1), it = 1, 15) /
     &        -38.000, -38.000, -38.000, -38.000, -38.000,               0.1
     &        -38.000, -38.000, -38.000, -38.000, -38.000,               0.1
     &        -38.000, -38.000, -38.000, -38.000, -38.000/               0.1
!
      data (crossch(it, 2), it = 1, 15) /
     &        -32.727, -31.151, -30.133, -29.432, -28.925,               0.2
     &        -28.547, -28.257, -28.030, -27.848, -27.701,               0.2
     &        -27.580, -27.479, -27.395, -27.322, -27.261/               0.2
!
      data (crossch(it, 3), it = 1, 15) /
     &        -31.588, -30.011, -28.993, -28.290, -27.784,               0.3
     &        -27.405, -27.115, -26.887, -26.705, -26.558,               0.3
     &        -26.437, -26.336, -26.251, -26.179, -26.117/               0.3
!
      data (crossch(it, 4), it = 1, 15) /
     &        -30.407, -28.830, -27.811, -27.108, -26.601,               0.4
     &        -26.223, -25.932, -25.705, -25.523, -25.376,               0.4
     &        -25.255, -25.154, -25.069, -24.997, -24.935/               0.4
!
      data (crossch(it, 5), it = 1, 15) /
     &        -29.513, -27.937, -26.920, -26.218, -25.712,               0.5
     &        -25.334, -25.043, -24.816, -24.635, -24.487,               0.5
     &        -24.366, -24.266, -24.181, -24.109, -24.047/               0.5
!
      data (crossch(it, 6), it = 1, 15) /
     &        -28.910, -27.341, -26.327, -25.628, -25.123,               0.6
     &        -24.746, -24.457, -24.230, -24.049, -23.902,               0.6
     &        -23.782, -23.681, -23.597, -23.525, -23.464/               0.6
!
      data (crossch(it, 7), it = 1, 15) /
     &        -28.517, -26.961, -25.955, -25.261, -24.760,               0.7
     &        -24.385, -24.098, -23.873, -23.694, -23.548,               0.7
     &        -23.429, -23.329, -23.245, -23.174, -23.113/               0.7
!
      data (crossch(it, 8), it = 1, 15) /
     &        -28.213, -26.675, -25.680, -24.993, -24.497,               0.8
     &        -24.127, -23.843, -23.620, -23.443, -23.299,               0.8
     &        -23.181, -23.082, -22.999, -22.929, -22.869/               0.8
!
      data (crossch(it, 9), it = 1, 15) /
     &        -27.942, -26.427, -25.446, -24.769, -24.280,               0.9
     &        -23.915, -23.635, -23.416, -23.241, -23.100,               0.9
     &        -22.983, -22.887, -22.805, -22.736, -22.677/               0.9
!
      data (crossch(it, 10), it = 1, 15) /
     &        -27.706, -26.210, -25.241, -24.572, -24.088,               1.0
     &        -23.728, -23.451, -23.235, -23.063, -22.923,               1.0
     &        -22.808, -22.713, -22.633, -22.565, -22.507/               1.0
!
      data (crossch(it, 11), it = 1, 15) /
     &        -27.475, -26.000, -25.043, -24.382, -23.905,               1.1
     &        -23.548, -23.275, -23.062, -22.891, -22.753,               1.1
     &        -22.640, -22.546, -22.467, -22.400, -22.343/               1.1
!
      data (crossch(it, 12), it = 1, 15) /
     &        -27.221, -25.783, -24.844, -24.193, -23.723,               1.2
     &        -23.372, -23.102, -22.892, -22.724, -22.588,               1.2
     &        -22.476, -22.384, -22.306, -22.240, -22.184/               1.2
!
      data (crossch(it, 13), it = 1, 15) /
     &        -26.863, -25.506, -24.607, -23.979, -23.523,               1.3
     &        -23.182, -22.919, -22.714, -22.550, -22.417,               1.3
     &        -22.309, -22.218, -22.142, -22.078, -22.023/               1.3
!
      data (crossch(it, 14), it = 1, 15) /
     &        -26.685, -25.347, -24.457, -23.835, -23.382,               1.4
     &        -23.044, -22.784, -22.580, -22.418, -22.286,               1.4
     &        -22.178, -22.089, -22.014, -21.950, -21.896/               1.4
!
      data (crossch(it, 15), it = 1, 15) /
     &        -26.085, -24.903, -24.105, -23.538, -23.120,               1.5
     &        -22.805, -22.561, -22.370, -22.217, -22.093,               1.5
     &        -21.991, -21.906, -21.835, -21.775, -21.723/               1.5
!
      data (crossch(it, 16), it = 1, 15) /
     &        -25.902, -24.727, -23.936, -23.376, -22.964,               1.6
     &        -22.654, -22.415, -22.227, -22.076, -21.955,               1.6
     &        -21.855, -21.772, -21.702, -21.644, -21.593/               1.6
!
      data (crossch(it, 17), it = 1, 15) /
     &        -25.215, -24.196, -23.510, -23.019, -22.655,               1.7
     &        -22.378, -22.163, -21.992, -21.855, -21.744,               1.7
     &        -21.653, -21.577, -21.513, -21.459, -21.412/               1.7
!
      data (crossch(it, 18), it = 1, 15) /
     &        -24.914, -23.937, -23.284, -22.820, -22.475,               1.8
     &        -22.212, -22.007, -21.845, -21.715, -21.609,               1.8
     &        -21.522, -21.449, -21.388, -21.336, -21.292/               1.8
!
      data (crossch(it, 19), it = 1, 15) /
     &        -24.519, -23.637, -23.039, -22.606, -22.281,               1.9
     &        -22.030, -21.834, -21.678, -21.552, -21.450,               1.9
     &        -21.365, -21.295, -21.236, -21.185, -21.142/               1.9
!
      data (crossch(it, 20), it = 1, 15) /
     &        -24.086, -23.222, -22.650, -22.246, -21.948,               2.0
     &        -21.722, -21.546, -21.407, -21.296, -21.205,               2.0
     &        -21.131, -21.070, -21.018, -20.974, -20.937/               2.0
!
      data (crossch(it, 21), it = 1, 15) /
     &        -23.850, -23.018, -22.472, -22.088, -21.805,               2.1
     &        -21.590, -21.422, -21.289, -21.182, -21.095,               2.1
     &        -21.024, -20.964, -20.914, -20.872, -20.835/               2.1
!
      data (crossch(it, 22), it = 1, 15) /
     &        -23.136, -22.445, -21.994, -21.676, -21.440,               2.2
     &        -21.259, -21.117, -21.004, -20.912, -20.837,               2.2
     &        -20.775, -20.723, -20.679, -20.642, -20.611/               2.2
!
      data (crossch(it, 23), it = 1, 15) /
     &        -23.199, -22.433, -21.927, -21.573, -21.314,               2.3
     &        -21.119, -20.969, -20.851, -20.758, -20.682,               2.3
     &        -20.621, -20.571, -20.529, -20.493, -20.463/               2.3
!
      data (crossch(it, 24), it = 1, 15) /
     &        -22.696, -22.020, -21.585, -21.286, -21.071,               2.4
     &        -20.912, -20.791, -20.697, -20.622, -20.563,               2.4
     &        -20.514, -20.475, -20.442, -20.414, -20.391/               2.4
!
      data (crossch(it, 25), it = 1, 15) /
     &        -22.119, -21.557, -21.194, -20.943, -20.761,               2.5
     &        -20.624, -20.518, -20.434, -20.367, -20.313,               2.5
     &        -20.268, -20.231, -20.201, -20.175, -20.153/               2.5
!
      data (crossch(it, 26), it = 1, 15) /
     &        -21.855, -21.300, -20.931, -20.673, -20.485,               2.6
     &        -20.344, -20.235, -20.151, -20.084, -20.031,               2.6
     &        -19.988, -19.953, -19.924, -19.900, -19.880/               2.6
!
      data (crossch(it, 27), it = 1, 15) /
     &        -21.126, -20.673, -20.382, -20.184, -20.044,               2.7
     &        -19.943, -19.868, -19.811, -19.769, -19.736,               2.7
     &        -19.710, -19.690, -19.674, -19.662, -19.652/               2.7
!
      data (crossch(it, 28), it = 1, 15) /
     &        -20.502, -20.150, -19.922, -19.766, -19.657,               2.8
     &        -19.578, -19.520, -19.478, -19.446, -19.422,               2.8
     &        -19.404, -19.390, -19.379, -19.371, -19.365/               2.8
!
      data (crossch(it, 29), it = 1, 15) /
     &        -20.030, -19.724, -19.530, -19.399, -19.309,               2.9
     &        -19.245, -19.199, -19.166, -19.142, -19.125,               2.9
     &        -19.112, -19.103, -19.096, -19.091, -19.088/               2.9
!
      data (crossch(it, 30), it = 1, 15) /
     &        -19.640, -19.364, -19.189, -19.074, -18.996,               3.0
     &        -18.943, -18.906, -18.881, -18.863, -18.852,               3.0
     &        -18.844, -18.839, -18.837, -18.836, -18.836/               3.0
!
      data (crossch(it, 31), it = 1, 15) /
     &        -19.333, -19.092, -18.939, -18.838, -18.770,               3.1
     &        -18.725, -18.695, -18.675, -18.662, -18.655,               3.1
     &        -18.651, -18.649, -18.649, -18.651, -18.653/               3.1
!
      data (crossch(it, 32), it = 1, 15) /
     &        -19.070, -18.880, -18.756, -18.674, -18.621,               3.2
     &        -18.585, -18.562, -18.548, -18.540, -18.536,               3.2
     &        -18.536, -18.537, -18.539, -18.542, -18.546/               3.2
!
      data (crossch(it, 33), it = 1, 15) /
     &        -18.851, -18.708, -18.617, -18.558, -18.521,               3.3
     &        -18.498, -18.484, -18.477, -18.475, -18.476,               3.3
     &        -18.478, -18.482, -18.487, -18.493, -18.498/               3.3
!
      data (crossch(it, 34), it = 1, 15) /
     &        -18.709, -18.599, -18.533, -18.494, -18.471,               3.4
     &        -18.459, -18.454, -18.454, -18.457, -18.462,               3.4
     &        -18.469, -18.476, -18.483, -18.490, -18.498/               3.4
!
      data (crossch(it, 35), it = 1, 15) /
     &        -18.656, -18.572, -18.524, -18.497, -18.485,               3.5
     &        -18.480, -18.482, -18.486, -18.493, -18.501,               3.5
     &        -18.510, -18.519, -18.527, -18.536, -18.544/               3.5
!
      data (crossch(it, 36), it = 1, 15) /
     &        -18.670, -18.613, -18.582, -18.566, -18.561,               3.6
     &        -18.562, -18.568, -18.575, -18.583, -18.592,               3.6
     &        -18.601, -18.610, -18.619, -18.627, -18.635/               3.6
!
      data (crossch(it, 37), it = 1, 15) /
     &        -18.728, -18.700, -18.687, -18.683, -18.685,               3.7
     &        -18.691, -18.698, -18.706, -18.715, -18.723,               3.7
     &        -18.731, -18.739, -18.745, -18.752, -18.758/               3.7
!
      data (crossch(it, 38), it = 1, 15) /
     &        -18.839, -18.835, -18.836, -18.842, -18.849,               3.8
     &        -18.857, -18.865, -18.872, -18.878, -18.883,               3.8
     &        -18.888, -18.892, -18.895, -18.898, -18.900/               3.8
!
      data (crossch(it, 39), it = 1, 15) /
     &        -19.034, -19.041, -19.049, -19.057, -19.064,               3.9
     &        -19.069, -19.071, -19.071, -19.070, -19.068,               3.9
     &        -19.065, -19.061, -19.058, -19.054, -19.051/               3.9
!
      data (crossch(it, 40), it = 1, 15) /
     &        -19.372, -19.378, -19.382, -19.380, -19.372,               4.0
     &        -19.359, -19.341, -19.321, -19.300, -19.280,               4.0
     &        -19.261, -19.243, -19.227, -19.212, -19.199/               4.0
!
      data (crossch(it, 41), it = 1, 15) /
     &        -19.780, -19.777, -19.763, -19.732, -19.686,               4.1
     &        -19.631, -19.573, -19.517, -19.465, -19.419,               4.1
     &        -19.379, -19.344, -19.314, -19.288, -19.265/               4.1
!
      data (crossch(it, 42), it = 1, 15) /
     &        -20.151, -20.133, -20.087, -20.009, -19.911,               4.2
     &        -19.810, -19.715, -19.631, -19.559, -19.497,               4.2
     &        -19.446, -19.402, -19.365, -19.333, -19.306/               4.2
!
      data (crossch(it, 43), it = 1, 15) /
     &        -20.525, -20.454, -20.312, -20.138, -19.970,               4.3
     &        -19.825, -19.705, -19.607, -19.528, -19.464,               4.3
     &        -19.411, -19.367, -19.330, -19.300, -19.274/               4.3
!
      data (crossch(it, 44), it = 1, 15) /
     &        -20.869, -20.655, -20.366, -20.104, -19.894,               4.4
     &        -19.731, -19.604, -19.505, -19.426, -19.363,               4.4
     &        -19.312, -19.271, -19.236, -19.208, -19.184/               4.4
!
      data (crossch(it, 45), it = 1, 15) /
     &        -21.179, -20.768, -20.380, -20.081, -19.856,               4.5
     &        -19.686, -19.556, -19.454, -19.375, -19.311,               4.5
     &        -19.260, -19.218, -19.184, -19.155, -19.131/               4.5
!
      data (crossch(it, 46), it = 1, 15) /
     &        -21.167, -20.601, -20.206, -19.925, -19.719,               4.6
     &        -19.565, -19.447, -19.355, -19.283, -19.226,               4.6
     &        -19.180, -19.143, -19.112, -19.087, -19.066/               4.6
!
      data (crossch(it, 47), it = 1, 15) /
     &        -20.918, -20.348, -19.976, -19.720, -19.536,               4.7
     &        -19.401, -19.299, -19.220, -19.159, -19.112,               4.7
     &        -19.073, -19.043, -19.018, -18.998, -18.981/               4.7
!
      data (crossch(it, 48), it = 1, 15) /
     &        -20.753, -20.204, -19.847, -19.602, -19.427,               4.8
     &        -19.299, -19.203, -19.129, -19.072, -19.028,               4.8
     &        -18.993, -18.965, -18.942, -18.924, -18.909/               4.8
!
      data (crossch(it, 49), it = 1, 15) /
     &        -20.456, -19.987, -19.677, -19.460, -19.302,               4.9
     &        -19.186, -19.098, -19.030, -18.978, -18.937,               4.9
     &        -18.904, -18.878, -18.857, -18.841, -18.827/               4.9
!
      data (crossch(it, 50), it = 1, 15) /
     &        -20.154, -19.734, -19.461, -19.272, -19.136,               5.0
     &        -19.035, -18.960, -18.902, -18.858, -18.824,               5.0
     &        -18.797, -18.775, -18.759, -18.745, -18.735/               5.0
!
      data (crossch(it, 51), it = 1, 15) /
     &        -19.941, -19.544, -19.288, -19.114, -18.992,               5.1
     &        -18.903, -18.837, -18.788, -18.751, -18.723,               5.1
     &        -18.701, -18.684, -18.671, -18.661, -18.654/               5.1
!
      data (crossch(it, 52), it = 1, 15) /
     &        -19.657, -19.321, -19.104, -18.956, -18.853,               5.2
     &        -18.779, -18.724, -18.684, -18.655, -18.632,               5.2
     &        -18.615, -18.602, -18.592, -18.585, -18.579/               5.2
!
      data (crossch(it, 53), it = 1, 15) /
     &        -19.388, -19.109, -18.930, -18.810, -18.725,               5.3
     &        -18.664, -18.620, -18.586, -18.562, -18.543,               5.3
     &        -18.529, -18.518, -18.510, -18.503, -18.498/               5.3
!
      data (crossch(it, 54), it = 1, 15) /
     &        -19.201, -18.953, -18.794, -18.686, -18.611,               5.4
     &        -18.556, -18.515, -18.485, -18.462, -18.446,               5.4
     &        -18.433, -18.423, -18.416, -18.410, -18.406/               5.4
!
      data (crossch(it, 55), it = 1, 15) /
     &        -18.923, -18.719, -18.588, -18.500, -18.439,               5.5
     &        -18.396, -18.365, -18.344, -18.328, -18.318,               5.5
     &        -18.311, -18.307, -18.304, -18.303, -18.302/               5.5
!
      data (crossch(it, 56), it = 1, 15) /
     &        -18.614, -18.458, -18.361, -18.298, -18.258,               5.6
     &        -18.232, -18.216, -18.206, -18.202, -18.201,               5.6
     &        -18.202, -18.205, -18.208, -18.213, -18.218/               5.6
!
      data (crossch(it, 57), it = 1, 15) /
     &        -18.419, -18.295, -18.222, -18.178, -18.153,               5.7
     &        -18.139, -18.132, -18.131, -18.133, -18.138,               5.7
     &        -18.143, -18.150, -18.157, -18.164, -18.172/               5.7
!
      data (crossch(it, 58), it = 1, 15) /
     &        -18.296, -18.201, -18.148, -18.118, -18.101,               5.8
     &        -18.094, -18.091, -18.093, -18.096, -18.101,               5.8
     &        -18.107, -18.113, -18.120, -18.126, -18.132/               5.8
!
      data (crossch(it, 59), it = 1, 15) /
     &        -18.021, -17.992, -17.977, -17.970, -17.967,               5.9
     &        -17.968, -17.970, -17.974, -17.978, -17.983,               5.9
     &        -17.989, -17.994, -18.000, -18.005, -18.011/               5.9
!
      data (crossch(it, 60), it = 1, 15) /
     &        -17.694, -17.686, -17.686, -17.691, -17.698,               6.0
     &        -17.708, -17.718, -17.729, -17.740, -17.750,               6.0
     &        -17.761, -17.771, -17.781, -17.790, -17.798/               6.0
!
      data (crossch(it, 61), it = 1, 15) /
     &        -17.374, -17.384, -17.400, -17.420, -17.440,               6.1
     &        -17.462, -17.483, -17.503, -17.523, -17.541,               6.1
     &        -17.558, -17.575, -17.590, -17.603, -17.616/               6.1
!
      data (crossch(it, 62), it = 1, 15) /
     &        -17.169, -17.199, -17.230, -17.262, -17.293,               6.2
     &        -17.323, -17.351, -17.378, -17.404, -17.427,               6.2
     &        -17.449, -17.469, -17.488, -17.505, -17.520/               6.2
!
      data (crossch(it, 63), it = 1, 15) /
     &        -17.151, -17.184, -17.217, -17.250, -17.282,               6.3
     &        -17.313, -17.342, -17.369, -17.395, -17.418,               6.3
     &        -17.440, -17.461, -17.480, -17.497, -17.513/               6.3
!
      data (crossch(it, 64), it = 1, 15) /
     &        -17.230, -17.260, -17.290, -17.320, -17.348,               6.4
     &        -17.375, -17.401, -17.425, -17.448, -17.469,               6.4
     &        -17.489, -17.508, -17.525, -17.541, -17.556/               6.4
!
      data (crossch(it, 65), it = 1, 15) /
     &        -17.379, -17.403, -17.425, -17.446, -17.467,               6.5
     &        -17.486, -17.505, -17.524, -17.541, -17.558,               6.5
     &        -17.574, -17.588, -17.602, -17.615, -17.627/               6.5
!
      data (crossch(it, 66), it = 1, 15) /
     &        -17.596, -17.604, -17.609, -17.612, -17.616,               6.6
     &        -17.622, -17.628, -17.636, -17.644, -17.652,               6.6
     &        -17.661, -17.670, -17.679, -17.687, -17.695/               6.6
!
      data (crossch(it, 67), it = 1, 15) /
     &        -17.846, -17.823, -17.795, -17.770, -17.750,               6.7
     &        -17.735, -17.725, -17.719, -17.716, -17.715,               6.7
     &        -17.716, -17.719, -17.722, -17.726, -17.730/               6.7
!
      data (crossch(it, 68), it = 1, 15) /
     &        -18.089, -18.015, -17.942, -17.882, -17.836,               6.8
     &        -17.802, -17.777, -17.760, -17.748, -17.740,               6.8
     &        -17.736, -17.734, -17.733, -17.734, -17.736/               6.8
!
      data (crossch(it, 69), it = 1, 15) /
     &        -18.299, -18.156, -18.038, -17.947, -17.881,               6.9
     &        -17.833, -17.798, -17.774, -17.757, -17.745,               6.9
     &        -17.738, -17.733, -17.730, -17.729, -17.729/               6.9
!
      data (crossch(it, 70), it = 1, 15) /
     &        -18.441, -18.243, -18.096, -17.991, -17.915,               7.0
     &        -17.860, -17.821, -17.792, -17.772, -17.757,               7.0
     &        -17.746, -17.738, -17.733, -17.730, -17.728/               7.0
!
      data (crossch(it, 71), it = 1, 15) /
     &        -18.474, -18.262, -18.111, -18.004, -17.926,               7.1
     &        -17.869, -17.826, -17.795, -17.771, -17.753,               7.1
     &        -17.740, -17.730, -17.722, -17.717, -17.713/               7.1
!
      data (crossch(it, 72), it = 1, 15) /
     &        -18.387, -18.191, -18.053, -17.952, -17.878,               7.2
     &        -17.823, -17.782, -17.752, -17.729, -17.711,               7.2
     &        -17.698, -17.689, -17.681, -17.676, -17.672/               7.2
!
      data (crossch(it, 73), it = 1, 15) /
     &        -18.161, -17.990, -17.874, -17.793, -17.736,               7.3
     &        -17.696, -17.668, -17.648, -17.634, -17.625,               7.3
     &        -17.619, -17.616, -17.614, -17.614, -17.615/               7.3
!
      data (crossch(it, 74), it = 1, 15) /
     &        -17.908, -17.774, -17.690, -17.637, -17.604,               7.4
     &        -17.583, -17.572, -17.567, -17.566, -17.568,               7.4
     &        -17.571, -17.576, -17.581, -17.587, -17.593/               7.4
!
      data (crossch(it, 75), it = 1, 15) /
     &        -17.681, -17.589, -17.540, -17.515, -17.506,               7.5
     &        -17.505, -17.511, -17.520, -17.530, -17.542,               7.5
     &        -17.554, -17.566, -17.578, -17.589, -17.600/               7.5
!
      data (crossch(it, 76), it = 1, 15) /
     &        -17.647, -17.606, -17.584, -17.575, -17.573,               7.6
     &        -17.576, -17.582, -17.589, -17.597, -17.605,               7.6
     &        -17.614, -17.623, -17.631, -17.639, -17.646/               7.6
!
      data (crossch(it, 77), it = 1, 15) /
     &        -17.300, -17.291, -17.291, -17.297, -17.307,               7.7
     &        -17.319, -17.333, -17.347, -17.361, -17.375,               7.7
     &        -17.389, -17.402, -17.415, -17.427, -17.438/               7.7
!
      data (crossch(it, 78), it = 1, 15) /
     &        -16.786, -16.802, -16.825, -16.853, -16.883,               7.8
     &        -16.914, -16.944, -16.974, -17.003, -17.030,               7.8
     &        -17.055, -17.079, -17.101, -17.122, -17.141/               7.8
!
      data (crossch(it, 79), it = 1, 15) /
     &        -16.489, -16.533, -16.579, -16.625, -16.670,               7.9
     &        -16.713, -16.754, -16.793, -16.830, -16.864,               7.9
     &        -16.896, -16.925, -16.952, -16.977, -17.000/               7.9
!
      data (crossch(it, 80), it = 1, 15) /
     &        -16.694, -16.724, -16.756, -16.789, -16.823,               8.0
     &        -16.856, -16.888, -16.919, -16.949, -16.976,               8.0
     &        -17.002, -17.026, -17.048, -17.069, -17.088/               8.0
!
      data (crossch(it, 81), it = 1, 15) /
     &        -16.935, -16.951, -16.971, -16.993, -17.016,               8.1
     &        -17.040, -17.064, -17.088, -17.111, -17.132,               8.1
     &        -17.153, -17.172, -17.190, -17.206, -17.222/               8.1
!
      data (crossch(it, 82), it = 1, 15) /
     &        -17.200, -17.208, -17.220, -17.235, -17.251,               8.2
     &        -17.269, -17.286, -17.304, -17.322, -17.338,               8.2
     &        -17.354, -17.369, -17.384, -17.397, -17.409/               8.2
!
      data (crossch(it, 83), it = 1, 15) /
     &        -17.597, -17.591, -17.589, -17.590, -17.594,               8.3
     &        -17.600, -17.608, -17.617, -17.626, -17.635,               8.3
     &        -17.645, -17.654, -17.662, -17.671, -17.679/               8.3
!
      data (crossch(it, 84), it = 1, 15) /
     &        -18.166, -18.134, -18.107, -18.085, -18.068,               8.4
     &        -18.056, -18.047, -18.041, -18.038, -18.036,               8.4
     &        -18.035, -18.035, -18.036, -18.038, -18.039/               8.4
!
      data (crossch(it, 85), it = 1, 15) /
     &        -19.000, -18.917, -18.838, -18.770, -18.714,               8.5
     &        -18.669, -18.632, -18.603, -18.579, -18.560,               8.5
     &        -18.545, -18.532, -18.522, -18.514, -18.507/               8.5
!
      data (crossch(it, 86), it = 1, 15) /
     &        -20.313, -19.982, -19.754, -19.592, -19.472,               8.6
     &        -19.380, -19.309, -19.253, -19.208, -19.172,               8.6
     &        -19.143, -19.119, -19.099, -19.083, -19.069/               8.6
!
      data (crossch(it, 87), it = 1, 15) /
     &        -19.751, -19.611, -19.520, -19.461, -19.423,               8.7
     &        -19.398, -19.382, -19.372, -19.366, -19.364,               8.7
     &        -19.363, -19.364, -19.366, -19.368, -19.371/               8.7
!
      data (crossch(it, 88), it = 1, 15) /
     &        -19.581, -19.431, -19.337, -19.277, -19.240,               8.8
     &        -19.218, -19.207, -19.202, -19.203, -19.207,               8.8
     &        -19.212, -19.220, -19.228, -19.236, -19.245/               8.8
!
      data (crossch(it, 89), it = 1, 15) /
     &        -19.685, -19.506, -19.389, -19.311, -19.258,               8.9
     &        -19.222, -19.199, -19.184, -19.175, -19.170,               8.9
     &        -19.168, -19.169, -19.171, -19.174, -19.177/               8.9
!
      data (crossch(it, 90), it = 1, 15) /
     &        -19.977, -19.756, -19.606, -19.501, -19.425,               9.0
     &        -19.370, -19.330, -19.300, -19.278, -19.262,               9.0
     &        -19.250, -19.241, -19.235, -19.230, -19.227/               9.0
!
      data (crossch(it, 91), it = 1, 15) /
     &        -20.445, -20.158, -19.958, -19.815, -19.711,               9.1
     &        -19.633, -19.574, -19.528, -19.493, -19.465,               9.1
     &        -19.442, -19.425, -19.410, -19.398, -19.389/               9.1
!
      data (crossch(it, 92), it = 1, 15) /
     &        -20.980, -20.625, -20.391, -20.229, -20.110,               9.2
     &        -20.020, -19.949, -19.892, -19.846, -19.807,               9.2
     &        -19.775, -19.748, -19.724, -19.704, -19.687/               9.2
!
      data (crossch(it, 93), it = 1, 15) /
     &        -21.404, -21.023, -20.771, -20.594, -20.461,               9.3
     &        -20.358, -20.274, -20.205, -20.148, -20.099,               9.3
     &        -20.058, -20.022, -19.991, -19.965, -19.942/               9.3
!
      data (crossch(it, 94), it = 1, 15) /
     &        -21.309, -20.970, -20.753, -20.603, -20.495,               9.4
     &        -20.412, -20.348, -20.295, -20.252, -20.215,               9.4
     &        -20.185, -20.158, -20.135, -20.115, -20.098/               9.4
!
      data (crossch(it, 95), it = 1, 15) /
     &        -21.221, -20.906, -20.707, -20.574, -20.480,               9.5
     &        -20.412, -20.361, -20.322, -20.292, -20.268,               9.5
     &        -20.249, -20.233, -20.221, -20.210, -20.201/               9.5
!
      data (crossch(it, 96), it = 1, 15) /
     &        -21.441, -21.097, -20.878, -20.728, -20.623,               9.6
     &        -20.546, -20.489, -20.446, -20.413, -20.387,               9.6
     &        -20.368, -20.352, -20.340, -20.330, -20.322/               9.6
!
      data (crossch(it, 97), it = 1, 15) /
     &        -21.668, -21.305, -21.071, -20.911, -20.797,               9.7
     &        -20.713, -20.650, -20.602, -20.565, -20.536,               9.7
     &        -20.514, -20.496, -20.481, -20.470, -20.460/               9.7
!
      data (crossch(it, 98), it = 1, 15) /
     &        -21.926, -21.556, -21.316, -21.150, -21.031,               9.8
     &        -20.942, -20.874, -20.822, -20.782, -20.750,               9.8
     &        -20.724, -20.704, -20.687, -20.674, -20.663/               9.8
!
      data (crossch(it, 99), it = 1, 15) /
     &        -22.319, -21.937, -21.686, -21.510, -21.380,               9.9
     &        -21.282, -21.206, -21.147, -21.099, -21.061,               9.9
     &        -21.031, -21.006, -20.985, -20.968, -20.954/               9.9
!
      data (crossch(it, 100), it = 1, 15) /
     &        -22.969, -22.561, -22.288, -22.092, -21.945,              10.0
     &        -21.832, -21.743, -21.672, -21.616, -21.570,              10.0
     &        -21.533, -21.503, -21.477, -21.457, -21.439/              10.0
!
      data (crossch(it, 101), it = 1, 15) /
     &        -24.001, -23.527, -23.199, -22.957, -22.772,              10.1
     &        -22.629, -22.516, -22.427, -22.355, -22.297,              10.1
     &        -22.250, -22.212, -22.180, -22.153, -22.131/              10.1
!
      data (crossch(it, 102), it = 1, 15) /
     &        -24.233, -23.774, -23.477, -23.273, -23.128,              10.2
     &        -23.022, -22.943, -22.883, -22.837, -22.802,              10.2
     &        -22.774, -22.752, -22.735, -22.721, -22.710/              10.2
!
      data (crossch(it, 103), it = 1, 15) /
     &        -24.550, -23.913, -23.521, -23.266, -23.094,              10.3
     &        -22.976, -22.893, -22.836, -22.796, -22.768,              10.3
     &        -22.750, -22.737, -22.730, -22.726, -22.725/              10.3
!
      data (crossch(it, 104), it = 1, 15) /
     &        -24.301, -23.665, -23.274, -23.019, -22.848,              10.4
     &        -22.730, -22.648, -22.591, -22.552, -22.525,              10.4
     &        -22.507, -22.495, -22.489, -22.485, -22.485/              10.4
!
      data (crossch(it, 105), it = 1, 15) /
     &        -24.519, -23.883, -23.491, -23.237, -23.065,              10.5
     &        -22.948, -22.866, -22.809, -22.770, -22.743,              10.5
     &        -22.724, -22.713, -22.706, -22.703, -22.702/              10.5
!
      data partch/
     &     203.741,  249.643,  299.341,  353.477,  412.607,  477.237,
     &     547.817,  624.786,  708.543,  799.463,  897.912, 1004.227,
     &    1118.738, 1241.761, 1373.588, 1514.481, 1664.677, 1824.394,
     &    1993.801, 2173.050, 2362.234, 2561.424, 2770.674, 2989.930,
     &    3219.204, 3458.378, 3707.355, 3966.005, 4234.155, 4511.604,
     &    4798.135, 5093.554, 5397.593, 5709.948, 6030.401, 6358.646,
     &    6694.379, 7037.313, 7387.147, 7743.579, 8106.313/
!
      data freq1 / 0.0 /
!
!------------------------------ EXECUTION ------------------------------
!
      chop = 0.0d0
!
      if (freq .ne. freq1) then
         freq1 = freq
         evolt = waveno / 8065.479d0
         n = evolt * 10.0d0
         en = dble(n) * 0.1
!
         if (n .ge. 20 .and. n .lt. 105) then
!
            do it = 1, 15
               crosscht(it) = (crossch(it, n) + 
     &                        (crossch(it, n+1) - crossch(it,n)) *
     &                        (evolt - en) / 0.1)
            end do
!
         end if
!
      end if
!
      if ((t(j) .lt. 9000.0) .and. (n .ge. 20) .and. (n .lt. 105)) then
         it = (t(j) - 1000.0d0) / 200.0d0 + 1.0
         it = max0(it, 1)
         tn = dble(it) * 200.0d0 + 800.0d0
         part = partch(it) + (partch(it+1) - partch(it)) * 
     &          (t(j) - tn) / 200.0d0
         it = (t(j) - 2000.0d0) / 500.0d0 + 1.0
         it = max0(it, 1)
         tn = dble(it) * 500.0d0 + 1500.0d0
         chop = exp((crosscht(it) + (crosscht(it+1) - crosscht(it)) * 
     &          (t(j) - tn) / 500.0d0) * tenlog) * part
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   C H O P !!!!!!!!!!!!!!!!
!
      double precision function fe1op (j)
      implicit none
!
!.... CROSS-SECTION TIMES PARTITION FUNCTION
!
!----------------------------- DUMMY VARIABLE --------------------------
!
      integer j
!
!------------------------------- COMMONS -------------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
      include 'common.stateb'
!
!----------------------------- LOCAL VARIABLES -------------------------
!
      double precision  bolt(48, maxd), e(48), freq1, g(48), wno(48), 
     &                  xsect(48)
      integer i, itemp1, k
      save bolt, e, freq1, g, itemp1, wno, xsect
!
!----------------------------- INITIALIZATION --------------------------
!
      data itemp1 / 0 /, freq1/ 0.0d0 /
!
      data e /     500.0,  7500.0, 12500.0, 17500.0, 19000.0, 19500.0, 
     &           19500.0, 21000.0, 22000.0, 23000.0, 23000.0, 24000.0, 
     &           24000.0, 24500.0, 24500.0, 26000.0, 26500.0, 26500.0, 
     &           27000.0, 27500.0, 28500.0, 29000.0, 29500.0, 29500.0, 
     &           29500.0, 30000.0, 31500.0, 31500.0, 33500.0, 33500.0, 
     &           34000.0, 34500.0, 34500.0, 35000.0, 35500.0, 37000.0, 
     &           37000.0, 37000.0, 38500.0, 40000.0, 40000.0, 41000.0, 
     &           41000.0, 43000.0, 43000.0, 43000.0, 43000.0, 44000.0/
!
      data g /      25.0,    35.0,    21.0,    15.0,     9.0,    35.0, 
     &              33.0,    21.0,    27.0,    49.0,     9.0,    21.0, 
     &              27.0,     9.0,     9.0,    25.0,    33.0,    15.0, 
     &              35.0,     3.0,     5.0,    11.0,    15.0,    13.0, 
     &              15.0,     9.0,    21.0,    15.0,    21.0,    25.0, 
     &              35.0,     9.0,     5.0,    45.0,    27.0,    21.0, 
     &              15.0,    21.0,    15.0,    25.0,    21.0,    35.0,  
     &               5.0,    15.0,    45.0,    35.0,    55.0,    25.0 /
!
      data wno / 63500.0, 58500.0, 53500.0, 59500.0, 45000.0, 44500.0, 
     &           44500.0, 43000.0, 58000.0, 41000.0, 54000.0, 40000.0, 
     &           40000.0, 57500.0, 55500.0, 38000.0, 57500.0, 57500.0, 
     &           37000.0, 54500.0, 53500.0, 55000.0, 34500.0, 34500.0, 
     &           34500.0, 34000.0, 32500.0, 32500.0, 32500.0, 32500.0, 
     &           32000.0, 29500.0, 29500.0, 31000.0, 30500.0, 29000.0, 
     &           27000.0, 54000.0, 27500.0, 24000.0, 47000.0, 23000.0, 
     &           44000.0, 42000.0, 42000.0, 21000.0, 42000.0, 42000.0/
!
!------------------------------- EXECUTION -- --------------------------
!
      if (itemp .ne. itemp1) then
         freq1 = 0.0d0
         itemp1=itemp
!
         do k = 1, nrhox
!
            do i = 1, 48
               bolt(i, k) = g(i) * exp(-e(i) * 2.99792458d10 * hkt(k))
            end do
!
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1 = freq
!
         if (waveno .ge. 21000.0d0) then
!
            do i = 1, 48
               xsect(i) = 0.0d0
               if (wno(i) .lt. waveno) xsect(i) =  3.0d-18 / 
     &            (1.0d0 + ((wno(i) + 3000.0d0 - waveno) / 
     &            wno(i) / 0.1d0) ** 4)
            end do
!
         end if
!
      end if
!
      fe1op = 0.0d0
!
      if (waveno .ge. 21000.0d0) then
!
         do i = 1, 48
            fe1op = fe1op + xsect(i) * bolt(i,j)
         end do
!
      end if
!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   F E 1 O P !!!!!!!!!!!!!!
!
      double precision function mg1op (j)
      implicit none
!
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
!
!.... 1995 OCT - TEST n .le. 7 BEFORE freq .le. freqmg TO AVOID AN ARRAY
!                OUT OF BOUNDS
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision   d, d1, dt(maxd), flog(9), freq1, freqmg(8), 
     &                   peach(7, 15), tlg(7), x(7)
      integer  it, itemp1, k, n, nt(maxd)
      save     dt, flog, freq1, freqmg, itemp1, nt, peach, tlg, x
!
!------------------------------ INITIALIZATION ------------------------
!
      data freq1, itemp1 / 0.0d0,  0 /
!
!.... IN THE ORDER 4000 K, 5000 K, 6000 K, 7000 K, 8000 K, 9000 K, 10000 K
!
      data flog / 35.23123, 35.19844, 35.15334, 34.71490, 34.31318, 
     &            33.75728, 33.65788, 33.64994, 33.43947 /
!
      data freqmg / 1.9341452d15, 1.8488510d15, 1.1925797d15,
     &              7.9804046d14, 4.5772110d14, 4.1440977d14,
     &              4.1113514d14 , 0.0/
!
      data peach /
!          4000     5000     6000     7000     8000     9000    10000   WAVE(A)
     &   -42.474, -42.350, -42.109, -41.795, -41.467, -41.159, -40.883,  1500
     &   -41.808, -41.735, -41.582, -41.363, -41.115, -40.866, -40.631,  1550
     &   -41.273, -41.223, -41.114, -40.951, -40.755, -40.549, -40.347,  1621
     &   -45.583, -44.008, -42.957, -42.205, -41.639, -41.198, -40.841,  1622
     &   -44.324, -42.747, -41.694, -40.939, -40.370, -39.925, -39.566,  2513
     &   -50.969, -48.388, -46.630, -45.344, -44.355, -43.568, -42.924,  2514
     &   -50.633, -48.026, -46.220, -44.859, -43.803, -42.957, -42.264,  3756
     &   -53.028, -49.643, -47.367, -45.729, -44.491, -43.520, -42.736,  3757
     &   -51.785, -48.352, -46.050, -44.393, -43.140, -42.157, -41.363,  6549
     &   -52.285, -48.797, -46.453, -44.765, -43.486, -42.480, -41.668,  6550
     &   -52.028, -48.540, -46.196, -44.507, -43.227, -42.222, -41.408,  7234
     &   -52.384, -48.876, -46.513, -44.806, -43.509, -42.488, -41.660,  7235
     &   -52.363, -48.856, -46.493, -44.786, -43.489, -42.467, -41.639,  7291
     &   -54.704, -50.772, -48.107, -46.176, -44.707, -43.549, -42.611,  7292
     &   -54.359, -50.349, -47.643, -45.685, -44.198, -43.027, -42.418 / 9000
!
      data tlg / 8.29405, 8.51719, 8.69951, 8.85367, 8.98720, 9.10498,
     &           9.21034 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         freq1  = 0.0d0
         itemp1 = itemp
!
         do k = 1, nrhox
            n     = max( min(6, int(t(k) / 1000.0d0) -3) , 1)
            nt(k) = n
            dt(k) = (tlog(k) - tlg(n)) / (tlg(n + 1) - tlg(n))
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1   = freq
         n       = 1
!
!.... THIS LOOP RETURNS N <= 8
!
         do while ((n .le. 7) .and. (freq .le. freqmg(n)))
            n = n + 1
         end do
!
         d = (freqln - flog(n)) / (flog(n + 1) - flog(n))
         if (n .gt. 2) n = 2 * n - 2
         d1 = 1.0d0 - d
!
         do it = 1, 7
            x(it) = peach(it, n + 1) * d + peach(it, n) * d1
         end do
!
      end if
!
      n     = nt(j)
      mg1op = exp(x(n) * (1.0d0 - dt(j)) + x(n + 1) * dt(j))
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   M G 1 O P !!!!!!!!!!!!!!
!
      double precision function ohop (j)
      implicit none
!
!.... CROSS-SECTION TIMES PARTITION FUNCTION
!
!----------------------------- DUMMY VARIABLE --------------------------
!
      integer, intent(in) ::  j
!
!------------------------------- COMMONS -------------------------------
!
      include 'common.constb'
      include 'common.tempbl'
      include 'common.stateb'
!
!----------------------------- LOCAL VARIABLES -------------------------
!
      double precision  crossoh(15, 130), crossoht(15), en, evolt, 
     &                  freq1, part, partoh(41), tn
      integer it, n
      save crossoh, crossoht, freq1, partoh, n 
!
!----------------------------- INITIALIZATION --------------------------
!
      data (crossoh(it, 1), it = 1, 15) /
     &        -30.855, -29.121, -27.976, -27.166, -26.566,               2.1
     &        -26.106, -25.742, -25.448, -25.207, -25.006,               2.1
     &        -24.836, -24.691, -24.566, -24.457, -24.363/               2.1
!
      data (crossoh(it, 2), it = 1, 15) /
     &        -30.494, -28.760, -27.615, -26.806, -26.206,               2.2
     &        -25.745, -25.381, -25.088, -24.846, -24.645,               2.2
     &        -24.475, -24.330, -24.205, -24.097, -24.002/               2.2
!
      data (crossoh(it, 3), it = 1, 15) /
     &        -30.157, -28.425, -27.280, -26.472, -25.872,               2.3
     &        -25.411, -25.048, -24.754, -24.513, -24.312,               2.3
     &        -24.142, -23.997, -23.872, -23.764, -23.669/               2.3
!
      data (crossoh(it, 4), it = 1, 15) /
     &        -29.848, -28.117, -26.974, -26.165, -25.566,               2.4
     &        -25.105, -24.742, -24.448, -24.207, -24.006,               2.4
     &        -23.836, -23.692, -23.567, -23.458, -23.364/               2.4
!
      data (crossoh(it, 5), it = 1, 15) /
     &        -29.567, -27.837, -26.693, -25.885, -25.286,               2.5
     &        -24.826, -24.462, -24.169, -23.928, -23.727,               2.5
     &        -23.557, -23.412, -23.287, -23.179, -23.084/               2.5
!
      data (crossoh(it, 6), it = 1, 15) /
     &        -29.307, -27.578, -26.436, -25.628, -25.029,               2.6
     &        -24.569, -24.205, -23.912, -23.671, -23.470,               2.6
     &        -23.300, -23.155, -23.031, -22.922, -22.828/               2.6
!
      data (crossoh(it, 7), it = 1, 15) /
     &        -29.068, -27.341, -26.199, -25.391, -24.792,               2.7
     &        -24.332, -23.969, -23.676, -23.435, -23.234,               2.7
     &        -23.064, -22.920, -22.795, -22.687, -22.592/               2.7
!
      data (crossoh(it, 8), it = 1, 15) /
     &        -28.820, -27.115, -25.978, -25.172, -24.574,               2.8
     &        -24.115, -23.752, -23.459, -23.218, -23.017,               2.8
     &        -22.848, -22.703, -22.579, -22.470, -22.376/               2.8
!
      data (crossoh(it, 9), it = 1, 15) /
     &        -28.540, -26.891, -25.768, -24.968, -24.372,               2.9
     &        -23.914, -23.552, -23.259, -23.019, -22.818,               2.9
     &        -22.649, -22.504, -22.380, -22.272, -22.177/               2.9
!
      data (crossoh(it, 10), it = 1, 15) /
     &        -28.275, -26.681, -25.574, -24.779, -24.186,               3.0
     &        -23.729, -23.368, -23.076, -22.836, -22.636,               3.0
     &        -22.467, -22.322, -22.198, -22.090, -21.996/               3.0
!
      data (crossoh(it, 11), it = 1, 15) /
     &        -27.993, -26.470, -25.388, -24.602, -24.014,               3.1
     &        -23.560, -23.200, -22.909, -22.669, -22.470,               3.1
     &        -22.301, -22.157, -22.033, -21.925, -21.831/               3.1
!
      data (crossoh(it, 12), it = 1, 15) /
     &        -27.698, -26.252, -25.204, -24.433, -23.851,               3.2
     &        -23.401, -23.043, -22.754, -22.515, -22.316,               3.2
     &        -22.148, -22.005, -21.881, -21.773, -21.679/               3.2
!
      data (crossoh(it, 13), it = 1, 15) /
     &        -27.398, -26.026, -25.019, -24.267, -23.696,               3.3
     &        -23.251, -22.896, -22.609, -22.372, -22.174,               3.3
     &        -22.007, -21.864, -21.741, -21.634, -21.540/               3.3
!
      data (crossoh(it, 14), it = 1, 15) /
     &        -27.100, -25.791, -24.828, -24.102, -23.543,               3.4
     &        -23.106, -22.756, -22.472, -22.238, -22.041,               3.4
     &        -21.875, -21.733, -21.611, -21.504, -21.411/               3.4
!
      data (crossoh(it, 15), it = 1, 15) /
     &        -26.807, -25.549, -24.631, -23.933, -23.391,               3.5
     &        -22.964, -22.621, -22.341, -22.109, -21.915,               3.5
     &        -21.751, -21.610, -21.488, -21.383, -21.290/               3.5
!
      data (crossoh(it, 16), it = 1, 15) /
     &        -26.531, -25.310, -24.431, -23.761, -23.238,               3.6
     &        -22.823, -22.488, -22.214, -21.986, -21.795,               3.6
     &        -21.633, -21.494, -21.374, -21.269, -21.178/               3.6
!
      data (crossoh(it, 17), it = 1, 15) /
     &        -26.239, -25.066, -24.225, -23.585, -23.082,               3.7
     &        -22.681, -22.356, -22.089, -21.866, -21.679,               3.7
     &        -21.520, -21.383, -21.265, -21.162, -21.072/               3.7
!
      data (crossoh(it, 18), it = 1, 15) /
     &        -25.945, -24.824, -24.017, -23.405, -22.923,               3.8
     &        -22.538, -22.223, -21.964, -21.748, -21.565,               3.8
     &        -21.410, -21.276, -21.160, -21.059, -20.970/               3.8
!
      data (crossoh(it, 19), it = 1, 15) /
     &        -25.663, -24.587, -23.810, -23.222, -22.761,               3.9
     &        -22.391, -22.088, -21.838, -21.629, -21.452,               3.9
     &        -21.300, -21.170, -21.057, -20.958, -20.872/               3.9
!
      data (crossoh(it, 20), it = 1, 15) /
     &        -25.372, -24.350, -23.603, -23.038, -22.596,               4.0
     &        -22.241, -21.950, -21.710, -21.508, -21.337,               4.0
     &        -21.190, -21.064, -20.954, -20.858, -20.774/               4.0
!
      data (crossoh(it, 21), it = 1, 15) /
     &        -25.076, -24.111, -23.396, -22.853, -22.429,               4.1
     &        -22.088, -21.809, -21.578, -21.384, -21.220,               4.1
     &        -21.078, -20.957, -20.851, -20.758, -20.676/               4.1
!
      data (crossoh(it, 22), it = 1, 15) /
     &        -24.779, -23.870, -23.189, -22.669, -22.261,               4.2
     &        -21.934, -21.667, -21.445, -21.259, -21.101,               4.2
     &        -20.965, -20.848, -20.746, -20.656, -20.578/               4.2
!
      data (crossoh(it, 23), it = 1, 15) /
     &        -24.486, -23.629, -22.983, -22.486, -22.095,               4.3
     &        -21.781, -21.524, -21.311, -21.132, -20.980,               4.3
     &        -20.850, -20.737, -20.639, -20.553, -20.478/               4.3
!
      data (crossoh(it, 24), it = 1, 15) /
     &        -24.183, -23.382, -22.774, -22.302, -21.928,               4.4
     &        -21.627, -21.381, -21.177, -21.005, -20.859,               4.4
     &        -20.734, -20.625, -20.531, -20.449, -20.376/               4.4
!
      data (crossoh(it, 25), it = 1, 15) /
     &        -23.867, -23.127, -22.561, -22.116, -21.761,               4.5
     &        -21.474, -21.238, -21.043, -20.878, -20.738,               4.5
     &        -20.617, -20.513, -20.423, -20.344, -20.274/               4.5
!
      data (crossoh(it, 26), it = 1, 15) /
     &        -23.538, -22.862, -22.340, -21.926, -21.592,               4.6
     &        -21.320, -21.096, -20.909, -20.751, -20.617,               4.6
     &        -20.502, -20.402, -20.315, -20.239, -20.172/               4.6
!
      data (crossoh(it, 27), it = 1, 15) /
     &        -23.234, -22.604, -22.120, -21.734, -21.422,               4.7
     &        -21.166, -20.953, -20.776, -20.625, -20.497,               4.7
     &        -20.387, -20.291, -20.208, -20.135, -20.071/               4.7
!
      data (crossoh(it, 28), it = 1, 15) /
     &        -22.934, -22.347, -21.898, -21.541, -21.250,               4.8
     &        -21.010, -20.811, -20.643, -20.500, -20.378,               4.8
     &        -20.273, -20.182, -20.102, -20.033, -19.971/               4.8
!
      data (crossoh(it, 29), it = 1, 15) /
     &        -22.637, -22.092, -21.676, -21.345, -21.075,               4.9
     &        -20.853, -20.666, -20.508, -20.374, -20.259,               4.9
     &        -20.159, -20.073, -19.997, -19.931, -19.872/               4.9
!
      data (crossoh(it, 30), it = 1, 15) /
     &        -22.337, -21.835, -21.452, -21.147, -20.899,               5.0
     &        -20.693, -20.520, -20.373, -20.247, -20.139,               5.0
     &        -20.046, -19.964, -19.892, -19.830, -19.774/               5.0
!
      data (crossoh(it, 31), it = 1, 15) /
     &        -22.049, -21.584, -21.230, -20.950, -20.721,               5.1
     &        -20.531, -20.372, -20.236, -20.119, -20.019,               5.1
     &        -19.931, -19.855, -19.788, -19.729, -19.676/               5.1
!
      data (crossoh(it, 32), it = 1, 15) /
     &        -21.768, -21.337, -21.011, -20.754, -20.544,               5.2
     &        -20.370, -20.223, -20.098, -19.991, -19.898,               5.2
     &        -19.817, -19.746, -19.683, -19.628, -19.579/               5.2
!
      data (crossoh(it, 33), it = 1, 15) /
     &        -21.494, -21.096, -20.796, -20.559, -20.367,               5.3
     &        -20.208, -20.074, -19.960, -19.861, -19.776,               5.3
     &        -19.701, -19.636, -19.578, -19.527, -19.482/               5.3
!
      data (crossoh(it, 34), it = 1, 15) /
     &        -21.233, -20.861, -20.585, -20.368, -20.193,               5.4
     &        -20.048, -19.926, -19.821, -19.732, -19.654,               5.4
     &        -19.586, -19.526, -19.473, -19.426, -19.384/               5.4
!
      data (crossoh(it, 35), it = 1, 15) /
     &        -20.983, -20.635, -20.380, -20.181, -20.021,               5.5
     &        -19.889, -19.778, -19.683, -19.602, -19.531,               5.5
     &        -19.469, -19.415, -19.367, -19.324, -19.286/               5.5
!
      data (crossoh(it, 36), it = 1, 15) /
     &        -20.743, -20.418, -20.182, -19.999, -19.853,               5.6
     &        -19.733, -19.633, -19.547, -19.474, -19.410,               5.6
     &        -19.354, -19.305, -19.261, -19.223, -19.189/               5.6
!
      data (crossoh(it, 37), it = 1, 15) /
     &        -20.515, -20.210, -19.991, -19.824, -19.690,               5.7
     &        -19.581, -19.490, -19.413, -19.347, -19.290,               5.7
     &        -19.240, -19.196, -19.157, -19.122, -19.092/               5.7
!
      data (crossoh(it, 38), it = 1, 15) /
     &        -20.297, -20.011, -19.808, -19.654, -19.532,               5.8
     &        -19.434, -19.352, -19.282, -19.223, -19.172,               5.8
     &        -19.127, -19.088, -19.054, -19.023, -18.996/               5.8
!
      data (crossoh(it, 39), it = 1, 15) /
     &        -20.090, -19.822, -19.633, -19.491, -19.381,               5.9
     &        -19.291, -19.218, -19.156, -19.103, -19.057,               5.9
     &        -19.018, -18.983, -18.952, -18.925, -18.901/               5.9
!
      data (crossoh(it, 40), it = 1, 15) /
     &        -19.893, -19.642, -19.467, -19.337, -19.236,               6.0
     &        -19.155, -19.089, -19.034, -18.987, -18.946,               6.0
     &        -18.912, -18.881, -18.854, -18.831, -18.810/               6.0
!
      data (crossoh(it, 41), it = 1, 15) /
     &        -19.705, -19.472, -19.309, -19.190, -19.098,               6.1
     &        -19.025, -18.966, -18.917, -18.876, -18.840,               6.1
     &        -18.810, -18.783, -18.760, -18.739, -18.721/               6.1
!
      data (crossoh(it, 42), it = 1, 15) /
     &        -19.527, -19.310, -19.161, -19.051, -18.968,               6.2
     &        -18.903, -18.851, -18.807, -18.771, -18.740,               6.2
     &        -18.713, -18.690, -18.670, -18.653, -18.637/               6.2
!
      data (crossoh(it, 43), it = 1, 15) /
     &        -19.357, -19.159, -19.022, -18.922, -18.847,               6.3
     &        -18.789, -18.743, -18.704, -18.673, -18.646,               6.3
     &        -18.623, -18.603, -18.586, -18.571, -18.558/               6.3
!
      data (crossoh(it, 44), it = 1, 15) /
     &        -19.195, -19.016, -18.892, -18.803, -18.736,               6.4
     &        -18.684, -18.643, -18.610, -18.583, -18.560,               6.4
     &        -18.540, -18.523, -18.509, -18.496, -18.485/               6.4
!
      data (crossoh(it, 45), it = 1, 15) /
     &        -19.042, -18.883, -18.772, -18.693, -18.634,               6.5
     &        -18.589, -18.553, -18.525, -18.501, -18.481,               6.5
     &        -18.465, -18.451, -18.438, -18.428, -18.419/               6.5
!
      data (crossoh(it, 46), it = 1, 15) /
     &        -18.894, -18.758, -18.662, -18.593, -18.542,               6.6
     &        -18.503, -18.473, -18.448, -18.428, -18.412,               6.6
     &        -18.398, -18.386, -18.376, -18.367, -18.359/               6.6
!
      data (crossoh(it, 47), it = 1, 15) /
     &        -18.752, -18.639, -18.559, -18.501, -18.458,               6.7
     &        -18.426, -18.400, -18.380, -18.363, -18.350,               6.7
     &        -18.338, -18.328, -18.320, -18.313, -18.306/               6.7
!
      data (crossoh(it, 48), it = 1, 15) /
     &        -18.611, -18.523, -18.460, -18.415, -18.381,               6.8
     &        -18.355, -18.334, -18.318, -18.304, -18.293,               6.8
     &        -18.284, -18.276, -18.269, -18.263, -18.258/               6.8
!
      data (crossoh(it, 49), it = 1, 15) /
     &        -18.471, -18.408, -18.362, -18.329, -18.304,               6.9
     &        -18.285, -18.269, -18.257, -18.247, -18.238,               6.9
     &        -18.231, -18.224, -18.219, -18.214, -18.210/               6.9
!
      data (crossoh(it, 50), it = 1, 15) /
     &        -18.330, -18.290, -18.261, -18.239, -18.223,               7.0
     &        -18.211, -18.201, -18.192, -18.185, -18.179,               7.0
     &        -18.174, -18.169, -18.165, -18.162, -18.159/               7.0
!
      data (crossoh(it, 51), it = 1, 15) /
     &        -18.190, -18.168, -18.154, -18.143, -18.135,               7.1
     &        -18.129, -18.124, -18.120, -18.116, -18.112,               7.1
     &        -18.109, -18.106, -18.104, -18.102, -18.100/               7.1
!
      data (crossoh(it, 52), it = 1, 15) /
     &        -18.055, -18.047, -18.043, -18.042, -18.040,               7.2
     &        -18.039, -18.039, -18.038, -18.037, -18.036,               7.2
     &        -18.035, -18.034, -18.033, -18.033, -18.032/               7.2
!
      data (crossoh(it, 53), it = 1, 15) /
     &        -17.929, -17.931, -17.935, -17.939, -17.943,               7.3
     &        -17.946, -17.948, -17.950, -17.952, -17.953,               7.3
     &        -17.955, -17.956, -17.957, -17.958, -17.959/               7.3
!
      data (crossoh(it, 54), it = 1, 15) /
     &        -17.818, -17.826, -17.834, -17.842, -17.849,               7.4
     &        -17.855, -17.860, -17.865, -17.869, -17.872,               7.4
     &        -17.875, -17.878, -17.881, -17.883, -17.886/               7.4
!
      data (crossoh(it, 55), it = 1, 15) /
     &        -17.724, -17.736, -17.747, -17.758, -17.767,               7.5
     &        -17.775, -17.782, -17.788, -17.793, -17.798,               7.5
     &        -17.803, -17.807, -17.811, -17.815, -17.819/               7.5
!
      data (crossoh(it, 56), it = 1, 15) /
     &        -17.651, -17.665, -17.678, -17.690, -17.701,               7.6
     &        -17.710, -17.718, -17.725, -17.732, -17.738,               7.6
     &        -17.744, -17.749, -17.755, -17.760, -17.765/               7.6
!
      data (crossoh(it, 57), it = 1, 15) /
     &        -17.601, -17.615, -17.629, -17.642, -17.653,               7.7
     &        -17.663, -17.672, -17.680, -17.688, -17.695,               7.7
     &        -17.701, -17.708, -17.714, -17.720, -17.726/               7.7
!
      data (crossoh(it, 58), it = 1, 15) /
     &        -17.572, -17.587, -17.602, -17.614, -17.626,               7.8
     &        -17.636, -17.645, -17.654, -17.662, -17.670,               7.8
     &        -17.677, -17.684, -17.691, -17.698, -17.704/               7.8
!
      data (crossoh(it, 59), it = 1, 15) /
     &        -17.565, -17.581, -17.595, -17.607, -17.619,               7.9
     &        -17.629, -17.638, -17.647, -17.656, -17.664,               7.9
     &        -17.671, -17.679, -17.686, -17.693, -17.700/               7.9
!
      data (crossoh(it, 60), it = 1, 15) /
     &        -17.580, -17.594, -17.608, -17.620, -17.630,               8.0
     &        -17.640, -17.650, -17.658, -17.667, -17.675,               8.0
     &        -17.682, -17.690, -17.697, -17.704, -17.711/               8.0
!
      data (crossoh(it, 61), it = 1, 15) /
     &        -17.613, -17.626, -17.639, -17.649, -17.659,               8.1
     &        -17.669, -17.677, -17.686, -17.694, -17.701,               8.1
     &        -17.709, -17.716, -17.723, -17.730, -17.737/               8.1
!
      data (crossoh(it, 62), it = 1, 15) /
     &        -17.663, -17.675, -17.685, -17.695, -17.703,               8.2
     &        -17.711, -17.719, -17.727, -17.734, -17.741,               8.2
     &        -17.748, -17.755, -17.761, -17.768, -17.774/               8.2
!
      data (crossoh(it, 63), it = 1, 15) /
     &        -17.728, -17.737, -17.745, -17.752, -17.759,               8.3
     &        -17.766, -17.772, -17.778, -17.785, -17.791,               8.3
     &        -17.797, -17.803, -17.808, -17.814, -17.820/               8.3
!
      data (crossoh(it, 64), it = 1, 15) /
     &        -17.803, -17.809, -17.814, -17.818, -17.823,               8.4
     &        -17.828, -17.832, -17.837, -17.842, -17.847,               8.4
     &        -17.852, -17.856, -17.861, -17.866, -17.871/               8.4
!
      data (crossoh(it, 65), it = 1, 15) /
     &        -17.884, -17.886, -17.888, -17.889, -17.891,               8.5
     &        -17.893, -17.896, -17.899, -17.902, -17.905,               8.5
     &        -17.908, -17.912, -17.915, -17.919, -17.922/               8.5
!
      data (crossoh(it, 66), it = 1, 15) /
     &        -17.966, -17.964, -17.961, -17.959, -17.958,               8.6
     &        -17.958, -17.958, -17.959, -17.960, -17.961,               8.6
     &        -17.963, -17.964, -17.966, -17.968, -17.970/               8.6
!
      data (crossoh(it, 67), it = 1, 15) /
     &        -18.040, -18.034, -18.028, -18.023, -18.019,               8.7
     &        -18.016, -18.013, -18.012, -18.010, -18.010,               8.7
     &        -18.009, -18.009, -18.009, -18.009, -18.010/               8.7
!
      data (crossoh(it, 68), it = 1, 15) /
     &        -18.096, -18.087, -18.078, -18.071, -18.065,               8.8
     &        -18.059, -18.055, -18.051, -18.047, -18.045,               8.8
     &        -18.042, -18.040, -18.039, -18.037, -18.036/               8.8
!
      data (crossoh(it, 69), it = 1, 15) /
     &        -18.125, -18.115, -18.105, -18.097, -18.089,               8.9
     &        -18.082, -18.076, -18.070, -18.065, -18.061,               8.9
     &        -18.057, -18.053, -18.051, -18.048, -18.046/               8.9
!
      data (crossoh(it, 70), it = 1, 15) /
     &        -18.120, -18.112, -18.103, -18.095, -18.087,               9.0
     &        -18.079, -18.072, -18.066, -18.060, -18.055,               9.0
     &        -18.050, -18.046, -18.042, -18.039, -18.036/               9.0
!
      data (crossoh(it, 71), it = 1, 15) /
     &        -18.083, -18.078, -18.071, -18.064, -18.057,               9.1
     &        -18.050, -18.044, -18.037, -18.032, -18.026,               9.1
     &        -18.022, -18.017, -18.014, -18.010, -18.007/               9.1
!
      data (crossoh(it, 72), it = 1, 15) /
     &        -18.025, -18.022, -18.017, -18.012, -18.006,               9.2
     &        -18.000, -17.994, -17.989, -17.984, -17.979,               9.2
     &        -17.975, -17.971, -17.968, -17.965, -17.963/               9.2
!
      data (crossoh(it, 73), it = 1, 15) /
     &        -17.957, -17.955, -17.952, -17.948, -17.943,               9.3
     &        -17.938, -17.934, -17.929, -17.925, -17.922,               9.3
     &        -17.918, -17.916, -17.913, -17.911, -17.910/               9.3
!
      data (crossoh(it, 74), it = 1, 15) /
     &        -17.890, -17.889, -17.886, -17.882, -17.879,               9.4
     &        -17.875, -17.871, -17.867, -17.864, -17.862,               9.4
     &        -17.860, -17.858, -17.857, -17.856, -17.855/               9.4
!
      data (crossoh(it, 75), it = 1, 15) /
     &        -17.831, -17.829, -17.826, -17.822, -17.819,               9.5
     &        -17.815, -17.812, -17.810, -17.807, -17.806,               9.5
     &        -17.804, -17.803, -17.803, -17.803, -17.803/               9.5
!
      data (crossoh(it, 76), it = 1, 15) /
     &        -17.786, -17.782, -17.777, -17.773, -17.769,               9.6
     &        -17.766, -17.763, -17.761, -17.759, -17.758,               9.6
     &        -17.757, -17.757, -17.757, -17.758, -17.759/               9.6
!
      data (crossoh(it, 77), it = 1, 15) /
     &        -17.753, -17.747, -17.741, -17.735, -17.731,               9.7
     &        -17.727, -17.724, -17.722, -17.721, -17.720,               9.7
     &        -17.720, -17.720, -17.721, -17.722, -17.724/               9.7
!
      data (crossoh(it, 78), it = 1, 15) /
     &        -17.733, -17.724, -17.716, -17.709, -17.703,               9.8
     &        -17.699, -17.696, -17.694, -17.693, -17.692,               9.8
     &        -17.692, -17.693, -17.694, -17.695, -17.697/               9.8
!
      data (crossoh(it, 79), it = 1, 15) /
     &        -17.723, -17.711, -17.700, -17.691, -17.685,               9.9
     &        -17.680, -17.676, -17.674, -17.673, -17.672,               9.9
     &        -17.673, -17.673, -17.675, -17.676, -17.678/               9.9
!
      data (crossoh(it, 80), it = 1, 15) /
     &        -17.718, -17.702, -17.689, -17.679, -17.672,              10.0
     &        -17.667, -17.663, -17.660, -17.659, -17.659,              10.0
     &        -17.659, -17.660, -17.661, -17.663, -17.665/              10.0
!
      data (crossoh(it, 81), it = 1, 15) /
     &        -17.713, -17.695, -17.681, -17.670, -17.662,              10.1
     &        -17.656, -17.653, -17.650, -17.649, -17.649,              10.1
     &        -17.649, -17.650, -17.651, -17.653, -17.655/              10.1
!
      data (crossoh(it, 82), it = 1, 15) /
     &        -17.705, -17.686, -17.671, -17.660, -17.652,              10.2
     &        -17.647, -17.643, -17.641, -17.640, -17.640,              10.2
     &        -17.640, -17.641, -17.643, -17.645, -17.647/              10.2
!
      data (crossoh(it, 83), it = 1, 15) /
     &        -17.690, -17.671, -17.657, -17.647, -17.640,              10.3
     &        -17.635, -17.632, -17.630, -17.630, -17.630,              10.3
     &        -17.631, -17.632, -17.634, -17.636, -17.639/              10.3
!
      data (crossoh(it, 84), it = 1, 15) /
     &        -17.667, -17.649, -17.637, -17.629, -17.623,              10.4
     &        -17.619, -17.618, -17.617, -17.617, -17.618,              10.4
     &        -17.619, -17.621, -17.623, -17.626, -17.628/              10.4
!
      data (crossoh(it, 85), it = 1, 15) /
     &        -17.635, -17.621, -17.611, -17.605, -17.601,              10.5
     &        -17.600, -17.599, -17.599, -17.601, -17.602,              10.5
     &        -17.604, -17.607, -17.609, -17.612, -17.615/              10.5
!
      data (crossoh(it, 86), it = 1, 15) /
     &        -17.596, -17.585, -17.579, -17.576, -17.575,              10.6
     &        -17.575, -17.576, -17.578, -17.580, -17.582,              10.6
     &        -17.585, -17.588, -17.591, -17.595, -17.598/              10.6
!
      data (crossoh(it, 87), it = 1, 15) /
     &        -17.550, -17.544, -17.542, -17.542, -17.544,              10.7
     &        -17.546, -17.548, -17.552, -17.555, -17.558,              10.7
     &        -17.562, -17.566, -17.570, -17.573, -17.577/              10.7
!
      data (crossoh(it, 88), it = 1, 15) /
     &        -17.501, -17.500, -17.501, -17.504, -17.508,              10.8
     &        -17.513, -17.517, -17.521, -17.526, -17.530,              10.8
     &        -17.535, -17.539, -17.544, -17.548, -17.553/              10.8
!
      data (crossoh(it, 89), it = 1, 15) /
     &        -17.449, -17.452, -17.457, -17.463, -17.470,              10.9
     &        -17.476, -17.482, -17.488, -17.493, -17.499,              10.9
     &        -17.504, -17.509, -17.514, -17.519, -17.524/              10.9
!
      data (crossoh(it, 90), it = 1, 15) /
     &        -17.396, -17.403, -17.412, -17.420, -17.429,              11.0
     &        -17.437, -17.444, -17.451, -17.458, -17.464,              11.0
     &        -17.470, -17.476, -17.481, -17.487, -17.492/              11.0
!
      data (crossoh(it, 91), it = 1, 15) /
     &        -17.344, -17.355, -17.366, -17.377, -17.387,              11.1
     &        -17.396, -17.405, -17.413, -17.420, -17.427,              11.1
     &        -17.434, -17.440, -17.446, -17.452, -17.458/              11.1
!
      data (crossoh(it, 92), it = 1, 15) /
     &        -17.295, -17.307, -17.321, -17.333, -17.345,              11.2
     &        -17.355, -17.365, -17.373, -17.382, -17.389,              11.2
     &        -17.397, -17.404, -17.410, -17.417, -17.423/              11.2
!
      data (crossoh(it, 93), it = 1, 15) /
     &        -17.249, -17.264, -17.278, -17.292, -17.304,              11.3
     &        -17.316, -17.326, -17.335, -17.344, -17.352,              11.3
     &        -17.360, -17.368, -17.375, -17.382, -17.389/              11.3
!
      data (crossoh(it, 94), it = 1, 15) /
     &        -17.209, -17.225, -17.241, -17.255, -17.268,              11.4
     &        -17.280, -17.291, -17.301, -17.310, -17.319,              11.4
     &        -17.327, -17.335, -17.343, -17.350, -17.357/              11.4
!
      data (crossoh(it, 95), it = 1, 15) /
     &        -17.177, -17.194, -17.210, -17.225, -17.239,              11.5
     &        -17.251, -17.262, -17.272, -17.282, -17.291,              11.5
     &        -17.300, -17.308, -17.316, -17.324, -17.331/              11.5
!
      data (crossoh(it, 96), it = 1, 15) /
     &        -17.154, -17.172, -17.189, -17.204, -17.218,              11.6
     &        -17.230, -17.242, -17.252, -17.262, -17.272,              11.6
     &        -17.280, -17.289, -17.298, -17.306, -17.314/              11.6
!
      data (crossoh(it, 97), it = 1, 15) /
     &        -17.144, -17.162, -17.179, -17.194, -17.208,              11.7
     &        -17.220, -17.232, -17.242, -17.253, -17.262,              11.7
     &        -17.271, -17.280, -17.289, -17.297, -17.306/              11.7
!
      data (crossoh(it, 98), it = 1, 15) /
     &        -17.146, -17.164, -17.181, -17.196, -17.210,              11.8
     &        -17.222, -17.234, -17.245, -17.255, -17.265,              11.8
     &        -17.274, -17.283, -17.292, -17.301, -17.309/              11.8
!
      data (crossoh(it, 99), it = 1, 15) /
     &        -17.163, -17.180, -17.197, -17.212, -17.225,              11.9
     &        -17.237, -17.249, -17.260, -17.270, -17.280,              11.9
     &        -17.289, -17.298, -17.307, -17.316, -17.325/              11.9
!
      data (crossoh(it, 100), it = 1, 15) /
     &        -17.193, -17.211, -17.227, -17.241, -17.254,              12.0
     &        -17.266, -17.277, -17.288, -17.298, -17.308,              12.0
     &        -17.317, -17.327, -17.336, -17.345, -17.353/              12.0
!
      data (crossoh(it, 101), it = 1, 15) /
     &        -17.239, -17.256, -17.271, -17.284, -17.297,              12.1
     &        -17.309, -17.320, -17.330, -17.340, -17.350,              12.1
     &        -17.359, -17.369, -17.378, -17.387, -17.395/              12.1
!
      data (crossoh(it, 102), it = 1, 15) /
     &        -17.299, -17.315, -17.329, -17.342, -17.354,              12.2
     &        -17.365, -17.376, -17.386, -17.396, -17.405,              12.2
     &        -17.415, -17.424, -17.433, -17.442, -17.451/              12.2
!
      data (crossoh(it, 103), it = 1, 15) /
     &        -17.373, -17.388, -17.402, -17.414, -17.425,              12.3
     &        -17.436, -17.446, -17.456, -17.466, -17.475,              12.3
     &        -17.484, -17.493, -17.502, -17.511, -17.520/              12.3
!
      data (crossoh(it, 104), it = 1, 15) /
     &        -17.462, -17.476, -17.489, -17.500, -17.511,              12.4
     &        -17.521, -17.531, -17.541, -17.550, -17.559,              12.4
     &        -17.569, -17.578, -17.587, -17.595, -17.604/              12.4
!
      data (crossoh(it, 105), it = 1, 15) /
     &        -17.567, -17.581, -17.592, -17.603, -17.613,              12.5
     &        -17.623, -17.632, -17.641, -17.651, -17.660,              12.5
     &        -17.669, -17.678, -17.686, -17.695, -17.704/              12.5
!
      data (crossoh(it, 106), it = 1, 15) /
     &        -17.689, -17.701, -17.712, -17.722, -17.732,              12.6
     &        -17.741, -17.750, -17.759, -17.768, -17.777,              12.6
     &        -17.786, -17.795, -17.803, -17.812, -17.821/              12.6
!
      data (crossoh(it, 107), it = 1, 15) /
     &        -17.829, -17.840, -17.851, -17.860, -17.869,              12.7
     &        -17.878, -17.887, -17.896, -17.904, -17.913,              12.7
     &        -17.922, -17.930, -17.939, -17.948, -17.956/              12.7
!
      data (crossoh(it, 108), it = 1, 15) /
     &        -17.988, -18.000, -18.010, -18.019, -18.028,              12.8
     &        -18.036, -18.045, -18.053, -18.062, -18.070,              12.8
     &        -18.079, -18.087, -18.096, -18.104, -18.112/              12.8
!
      data (crossoh(it, 109), it = 1, 15) /
     &        -18.171, -18.183, -18.192, -18.201, -18.210,              12.9
     &        -18.218, -18.227, -18.235, -18.243, -18.252,              12.9
     &        -18.260, -18.268, -18.277, -18.285, -18.293/              12.9
!
      data (crossoh(it, 110), it = 1, 15) /
     &        -18.381, -18.393, -18.403, -18.413, -18.422,              13.0
     &        -18.430, -18.438, -18.447, -18.455, -18.463,              13.0
     &        -18.471, -18.479, -18.487, -18.495, -18.503/              13.0
!
      data (crossoh(it, 111), it = 1, 15) /
     &        -18.625, -18.638, -18.650, -18.660, -18.669,              13.1
     &        -18.678, -18.687, -18.695, -18.703, -18.711,              13.1
     &        -18.719, -18.726, -18.734, -18.742, -18.750/              13.1
!
      data (crossoh(it, 112), it = 1, 15) /
     &        -18.912, -18.929, -18.943, -18.955, -18.966,              13.2
     &        -18.975, -18.984, -18.993, -19.001, -19.008,              13.2
     &        -19.016, -19.023, -19.031, -19.038, -19.045/              13.2
!
      data (crossoh(it, 113), it = 1, 15) /
     &        -19.260, -19.283, -19.303, -19.320, -19.333,              13.3
     &        -19.345, -19.355, -19.364, -19.372, -19.380,              13.3
     &        -19.387, -19.394, -19.400, -19.407, -19.413/              13.3
!
      data (crossoh(it, 114), it = 1, 15) /
     &        -19.704, -19.740, -19.771, -19.796, -19.816,              13.4
     &        -19.832, -19.845, -19.855, -19.863, -19.870,              13.4
     &        -19.876, -19.882, -19.887, -19.892, -19.897/              13.4
!
      data (crossoh(it, 115), it = 1, 15) /
     &        -20.339, -20.386, -20.424, -20.454, -20.476,              13.5
     &        -20.492, -20.502, -20.509, -20.513, -20.516,              13.5
     &        -20.518, -20.520, -20.521, -20.523, -20.524/              13.5
!
      data (crossoh(it, 116), it = 1, 15) /
     &        -21.052, -21.075, -21.093, -21.105, -21.114,              13.6
     &        -21.120, -21.123, -21.125, -21.126, -21.127,              13.6
     &        -21.128, -21.130, -21.131, -21.133, -21.135/              13.6
!
      data (crossoh(it, 117), it = 1, 15) /
     &        -21.174, -21.203, -21.230, -21.255, -21.278,              13.7
     &        -21.299, -21.320, -21.339, -21.357, -21.375,              13.7
     &        -21.392, -21.408, -21.424, -21.439, -21.454/              13.7
!
      data (crossoh(it, 118), it = 1, 15) /
     &        -21.285, -21.317, -21.346, -21.372, -21.395,              13.8
     &        -21.416, -21.435, -21.452, -21.468, -21.483,              13.8
     &        -21.497, -21.511, -21.524, -21.536, -21.548/              13.8
!
      data (crossoh(it, 119), it = 1, 15) /
     &        -21.396, -21.429, -21.459, -21.486, -21.511,              13.9
     &        -21.532, -21.551, -21.569, -21.585, -21.600,              13.9
     &        -21.614, -21.627, -21.640, -21.652, -21.663/              13.9
!
      data (crossoh(it, 120), it = 1, 15) /
     &        -21.516, -21.549, -21.580, -21.609, -21.635,              14.0
     &        -21.658, -21.678, -21.696, -21.713, -21.728,              14.0
     &        -21.742, -21.755, -21.767, -21.779, -21.790/              14.0
!
      data (crossoh(it, 121), it = 1, 15) /
     &        -21.651, -21.681, -21.711, -21.738, -21.763,              14.1
     &        -21.785, -21.804, -21.821, -21.837, -21.851,              14.1
     &        -21.864, -21.876, -21.887, -21.898, -21.908/              14.1
!
      data (crossoh(it, 122), it = 1, 15) /
     &        -21.810, -21.831, -21.853, -21.874, -21.893,              14.2
     &        -21.910, -21.925, -21.938, -21.950, -21.961,              14.2
     &        -21.971, -21.980, -21.989, -21.998, -22.006/              14.2
!
      data (crossoh(it, 123), it = 1, 15) /
     &        -22.009, -22.016, -22.026, -22.037, -22.048,              14.3
     &        -22.058, -22.066, -22.074, -22.081, -22.088,              14.3
     &        -22.094, -22.099, -22.105, -22.111, -22.117/              14.3
!
      data (crossoh(it, 124), it = 1, 15) /
     &        -22.353, -22.317, -22.296, -22.284, -22.276,              14.4
     &        -22.270, -22.266, -22.262, -22.260, -22.258,              14.4
     &        -22.257, -22.257, -22.257, -22.258, -22.259/              14.4
!
      data (crossoh(it, 125), it = 1, 15) /
     &        -22.705, -22.609, -22.552, -22.515, -22.488,              14.5
     &        -22.468, -22.451, -22.438, -22.427, -22.418,              14.5
     &        -22.410, -22.405, -22.400, -22.397, -22.395/              14.5
!
      data (crossoh(it, 126), it = 1, 15) /
     &        -22.889, -22.791, -22.731, -22.690, -22.659,              14.6
     &        -22.634, -22.612, -22.594, -22.579, -22.566,              14.6
     &        -22.555, -22.546, -22.539, -22.533, -22.528/              14.6
!
      data (crossoh(it, 127), it = 1, 15) /
     &        -23.211, -23.109, -23.041, -22.989, -22.945,              14.7
     &        -22.906, -22.872, -22.842, -22.816, -22.793,              14.7
     &        -22.774, -22.757, -22.743, -22.732, -22.722/              14.7
!
      data (crossoh(it, 128), it = 1, 15) /
     &        -25.312, -24.669, -24.250, -23.959, -23.746,              14.8
     &        -23.587, -23.463, -23.366, -23.288, -23.225,              14.8
     &        -23.173, -23.131, -23.095, -23.066, -23.041/              14.8
!
      data (crossoh(it, 129), it = 1, 15) /
     &        -25.394, -24.752, -24.333, -24.041, -23.829,              14.9
     &        -23.669, -23.546, -23.449, -23.371, -23.308,              14.9
     &        -23.256, -23.214, -23.178, -23.149, -23.124/              14.9
!
      data (crossoh(it, 130), it = 1, 15) /
     &        -25.430, -24.787, -24.369, -24.077, -23.865,              15.0
     &        -23.705, -23.582, -23.484, -23.407, -23.344,              15.0
     &        -23.292, -23.249, -23.214, -23.185, -23.160/              15.0
!
      data partoh /
     &     145.979,  178.033,  211.618,  247.053,  284.584,  324.398,
     &     366.639,  411.425,  458.854,  509.012,  561.976,  617.823,
     &     676.626,  738.448,  803.363,  871.437,  942.735, 1017.330,
     &    1095.284, 1176.654, 1261.510, 1349.898, 1441.875, 1537.483,
     &    1636.753, 1739.733, 1846.434, 1956.883, 2071.080, 2189.029,
     &    2310.724, 2436.155, 2565.283, 2698.103, 2834.571, 2974.627,
     &    3118.242, 3265.366, 3415.912, 3569.837, 3727.077/
!
      data freq1 / 0.0d0 /
!
!--------------------------------- EXECUTION ---------------------------
!
      ohop = 0.0d0  
!
      if (freq .ne. freq1) then
         freq1  = freq
         evolt  = waveno / 8065.479d0
         n      = int(evolt * 10.0d0 - 20.0d0,kind(n) )
         en     = dble(n) * 0.1 + 2.0

         if (n .gt. 0 .and. n .lt. 130) then
!
            do it = 1, 15
               crossoht(it) = (crossoh(it,n) + 
     &                        (crossoh(it,n+1) - crossoh(it,n)) *
     &                        (evolt - en) / 0.1)
            end do
!
         end if
!
      end if
!
      if ((t(j) .lt. 9000.0)) then 

        if( (n .gt. 0) .and. (n .lt. 130)) then

         it = (t(j) - 1000.0d0) / 200.0d0 + 1.0
         it = max0(it, 1)
         tn = dble(it) * 200.0d0 + 800.0d0
         part = partoh(it) + (partoh(it+1) - partoh(it)) * 
     &          (t(j) - tn) / 200.0d0
         it = (t(j) - 2000.0d0) / 500.0d0 + 1.0
         it = max0(it, 1)
         tn = dble(it) * 500.0d0 + 1500.0d0
         ohop = exp((crossoht(it) + (crossoht(it+1) - crossoht(it)) *
     &          (t(j) - tn) / 500.0d0) * tenlog) * part
        end if 
      end if

      return

!
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   O H O P !!!!!!!!!!!!!!!!
!
      double precision function seaton (freq0, xsect, power, a)
      implicit none
!
!--------------------------- DUMMY VARIABLES ---------------------------
!
      double precision a, freq0, power, xsect
!
!
!------------------------------- EXECUTION -----------------------------
!
      seaton = xsect * (a + (1.0d0 - a) * (freq0 / freq)) *
     &         sqrt((freq0 / freq) ** (int(2.0d0 * power + 0.01d0)))

      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   S E A T O N !!!!!!!!!!!!
!
      double precision function si1op(j)
      implicit none
!
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
!
!.... 1995 OCT - TEST n .le. 9 BEFORE freq .le. freqsi TO AVOID AN ARRAY
!                OUT OF BOUNDS
!.... 1996 JAN - SAVED ALL RELEVANT LOCAL VARIABLES
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  d, d1, dt(maxd), flog(11), freq1, freqsi(10), 
     &                  peach(9, 19), tlg(9), x(9)
      integer  it, itemp1, k, n, nt(maxd)
      save     dt, flog, freq1, freqsi, peach, tlg, x, itemp1, nt

!
!------------------------------ INITIALIZATION ------------------------
!
      data freq1, itemp1 / 0.0d0, 0 /
!
      data flog / 35.45438, 35.30022, 35.21799, 35.11986, 34.95438,
     &            33.95402, 33.90947, 33.80244, 33.78835, 33.76626,
     &            33.70518 /
!
!.... 3P, 1D, 1S, 1D, 3D, 3F, 1D, 3P
!
      data freqsi / 2.1413750d15, 1.9723165d15, 1.7879689d15,
     &              1.5152920d15, 5.5723927d14, 5.3295914d14,
     &              4.7886458d14, 4.7216422d14, 4.6185133d14, 0.0 /
!
      data peach /
!....     4000   5000   6000   7000   8000   9000   10000  11000  12000  WAVE(A)
     &   38.136,38.138,38.140,38.141,38.143,38.144,38.144,38.145,38.145,  1200
     &   37.834,37.839,37.843,37.847,37.850,37.853,37.855,37.857,37.858,  1400
     &   37.898,37.898,37.897,37.897,37.897,37.896,37.895,37.895,37.894,  1519
     &   40.737,40.319,40.047,39.855,39.714,39.604,39.517,39.445,39.385,  1520
     &   40.581,40.164,39.893,39.702,39.561,39.452,39.366,39.295,39.235,  1676
     &   45.521,44.456,43.753,43.254,42.878,42.580,42.332,42.119,41.930,  1677
     &   45.520,44.455,43.752,43.251,42.871,42.569,42.315,42.094,42.896,  1978
     &   55.068,51.783,49.553,47.942,46.723,45.768,44.997,44.360,43.823,  1979
     &   53.868,50.369,48.031,46.355,45.092,44.104,43.308,42.652,42.100,  5379
     &   54.133,50.597,48.233,46.539,45.261,44.262,43.456,42.790,42.230,  5380
     &   54.051,50.514,48.150,46.454,45.176,44.175,43.368,42.702,42.141,  5624
     &   54.442,50.854,48.455,46.733,45.433,44.415,43.592,42.912,42.340,  5625
     &   54.320,50.722,48.313,46.583,45.277,44.251,43.423,42.738,42.160,  6260
     &   55.691,51.965,49.444,47.615,46.221,45.119,44.223,43.478,42.848,  6261
     &   55.661,51.933,49.412,47.582,46.188,45.085,44.189,43.445,42.813,  6349
     &   55.973,52.193,49.630,47.769,46.349,45.226,44.314,43.555,42.913,  6350
     &   55.922,52.141,49.577,47.715,46.295,45.172,44.259,43.500,42.858,  6491
     &   56.828,52.821,50.110,48.146,46.654,45.477,44.522,43.730,43.061,  6492
     &   56.657,52.653,49.944,47.983,46.491,45.315,44.360,43.569,42.901/  6900
!
      data tlg / 8.29405, 8.51719, 8.69951, 8.85367, 8.98720, 9.10498,
     &           9.21034, 9.30565, 9.39266 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         freq1  = 0.0d0
         itemp1 = itemp
!
         do k = 1, nrhox
            n     = max(min(8, int(t(k) / 1000.0d0) -3) , 1)
            nt(k) = n
            dt(k) = (tlog(k) - tlg(n)) / (tlg(n + 1) - tlg(n))
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1   = freq
         n       = 1
!
!.... THIS LOOP RETURNS N <= 10
!
         do while (n .le. 9  .and. freq .le. freqsi(n))
           n = n + 1
         end do
!
         d  = (freqln - flog(n)) / (flog(n + 1) - flog(n))
         d1 = 1.0 - d
         if (n .gt. 2) n = 2 * n - 2
!
         do it = 1, 9
            x(it) = peach(it, n + 1) * d + peach(it, n) * d1
         end do
!
      end if
!
      n     = nt(j)
      si1op = exp(-(x(n) * (1.0 - dt(j)) + x(n + 1) * dt(j))) * 9.0
      end function 
!
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   S I 1 O P !!!!!!!!!!!!!!
!
!--------------- Castelli version of si1op ----------------------------:


      subroutine si1opcastelli 

      implicit none 
!      include 'common.sizebl' 
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
      include 'common.opsblk'
      include 'common.constb'

!------ external 
!      external xkarsas, si1seaton3p3d3p,si1seaton3p4p1d,si1seaton3p4p3s
      
!---- local variables ------------------------------------------------!
      integer itemp1, l, i 
      double precision elev(33), glev(33), bolt(33,maxd), x(33) 
      double precision z, freq3, rydloc, zeff2, eps, reson1, degen, h1  
      double precision elim
!---- initialisation -------------------------------------------------!

      data elev/ 59962.284,    59100.,  59077.112,  58893.40, 58801.529,
     &              58777., 57488.974,  56503.346, 54225.621,  53387.34,
     &            53362.24, 51612.012,  50533.424, 50189.389, 49965.894,
     &           49399.670, 49128.131,  48161.459, 47351.554, 47284.061,
     &           45303.310, 40991.884,  39859.920, 15394.370,  6298.850,
     &             223.157,    77.115,      0.000,    94000.,   79664.0,
     &              72000., 56698.738,  33326.053/

      data glev/  9.,  56.,  15.,  7.,  3.,  28.,  21.,  5.,  15.,
     &            3.,   7.,   1.,  9.,  5.,  21.,   3.,  9.,  15.,
     &            5.,   3.,  15.,  3.,  9.,   1.,   5.,  5.,   3.,
     &            1.,   3.,   3.,  5., 12.,   5./

      data itemp1/ 0/

      if( itemp .ne. itemp1) then 
        itemp1=itemp
        do   l=1,nrhox
          do  i=1,33
            bolt(i,l)=glev(i)*exp(-elev(i)*hckt(l))
          end do  
        end do 

      end if

      z=1.0d0 
      freq3=2.815d29/freq/freq/freq*z**4
      rydloc = ryd/c_cm

      do i=1,33
        x(i)=0.0d0 
      end do 
!---------------------------
! elim1 =  65939.18  
! elim2 =  65747.55 
! elim3 =  65747.55+287.45 
! elim4 =  65747.5+42824.35 
! elim5 =  65747.55 
!---------------------------------------------------------------------!
!      si ii 3s2 3p 2p average
      elim=65939.18
!
!  3s2 3p4d 3p
      if(waveno.lt.elim-elev(1))go to 30
      zeff2=16./rydloc*(elim-elev(1))
      x(1)=xkarsas(zeff2,4,2)
!
!  3s2 3p4f (2p3/2)4f
      if(waveno.lt.elim-elev(2))go to 30
      zeff2=16./rydloc*(elim-elev(2))
      x(2)=xkarsas(zeff2,4,3)
!
!  3s2 3p4d 3d
      if(waveno.lt.elim-elev(3))go to 30
      zeff2=16./rydloc*(elim-elev(3))
      x(3)=xkarsas(zeff2,4,2)
!
!  3s2 3p4d 1f
      if(waveno.lt.elim-elev(4))go to 30
      zeff2=16./rydloc*(elim-elev(4))
      x(4)=xkarsas(zeff2,4,2)
!
!  3s2 3p4d 1p
!  elev=58801.529
      if(waveno.lt.elim-elev(5))go to 30
!     glev=3.
      zeff2=16./rydloc*(elim-elev(5))
      x(5)=xkarsas(zeff2,4,2)
!
!  3s2 3p4f (2p1/2)4f
!     elev=58777.
      if(waveno.lt.elim-elev(6))go to 30
!     glev=28.
      zeff2=16./rydloc*(elim-elev(6))
      x(6)=xkarsas(zeff2,4,3)
!
!  3s2 3p4d 3f
!  elev=57488.974
      if(waveno.lt.elim-elev(7))go to 30
!     glev=21.
      zeff2=16./rydloc*(elim-elev(7))
      x(7)=xkarsas(zeff2,4,2)
!-----------
!  3s2 3p4d 1d
!  elev=56503.346
      if(waveno .lt. elim-elev(8))go to 30
!     glev=5.
      zeff2=16./rydloc*(elim-elev(8))
      x(8)=xkarsas(zeff2,4,2)
!------------
!  3s2 3p3d 3d
!  elev=54225.621
      if(waveno .lt. elim-elev(9))go to 30
!     glev=15.
      zeff2=9./rydloc*(elim-elev(9))
      x(9)=xkarsas(zeff2,3,2)
!-------------
!  3s2 3p3d 1p
!  elev=53387.34
      if(waveno .lt. elim-elev(10))go to 30
!     glev=3.
      zeff2=9./rydloc*(elim-elev(10))
      x(10)=xkarsas(zeff2,3,2)
!------------
!  3s2 3p3d 1f
!  elev=53362.24
      if(waveno.lt.elim-elev(11))go to 30
!     glev=7.
      zeff2=9./rydloc*(elim-elev(11))
      x(11)=xkarsas(zeff2,3,2)
!-------------
!  3s2 3p4p 1s
!     elev=51612.012
      if(waveno.lt.elim-elev(12))go to 30
!     glev=1.
      zeff2=16./rydloc*(elim-elev(12))
      x(12)=xkarsas(zeff2,4,1)
!-------------
!  3s2 3p3d 3p
!  elev=50533.424
      if(waveno.lt.elim-elev(13))go to 30
!     glev=9.
      zeff2=9./rydloc*(elim-elev(13))
      x(13)=si1seaton3p3d3p()

!------------
!  3s2 3p4p 1d
!     elev=50189.389
      if(waveno.lt.elim-elev(14))go to 30
!     glev=5.
      zeff2=16./rydloc*(elim-elev(14))
      x(14)=si1seaton3p4p1d()
!-----------
!  3s2 3p3d 3f
!  elev=49965.894
      if(waveno.lt.elim-elev(15))go to 30
!     glev=21.
      zeff2=9./rydloc*(elim-elev(15))
      x(15)=42.68d-18*(15973.206/waveno)**2.4
!----------
!  3s2 3p4p 3s
!  elev=49399.670
      if(waveno.lt.elim-elev(16))go to 30
!     glev=3.
      zeff2=16./rydloc*(elim-elev(16))
      x(16)=si1seaton3p4p3s()
!-----------
!  3s2 3p4p 3p
!  elev=49128.131
      if(waveno.lt.elim-elev(17))go to 30
!     glev=9.
      zeff2=16./rydloc*(elim-elev(17))
      if(waveno.ge.16810.969) x(17)=15.67d-18*(16810.969/waveno)**0.7
      if(waveno.gt.23803.) x(17)=11.90d-18*(23803/waveno)**3.7

!--------------
!  3s2 3p4p 3d
!     elev=48161.459
      if(waveno.lt.elim-elev(18))go to 30
!     glev=15.
      zeff2=16./rydloc*(elim-elev(18))
      x(18)=18.7d-18*(17777.641/waveno)**2.5
!------------
!  3s2 3p3d 1d
!  elev=47351.554
      if(waveno.lt.elim-elev(19))go to 30
!     glev=5.
      zeff2=9./rydloc*(elim-elev(19))
      x(19)=37d-18*(18587.564/waveno)**2
!------------
!  2s2 3p4p 1p
!  elev=47284.061
      if(waveno.lt.elim-elev(20))go to 30
!     glev=3.
      zeff2=16./rydloc*(elim-elev(20))
      x(20)=12.41d-18*(18655.039/waveno)**2.3
!------------
!  2s2p3 3d
!  elev=45303.310
      if(waveno.lt.elim-elev(21))go to 60
!     glev=15.
      zeff2=9./rydloc*(elim-elev(21))
      x(21)=16.64d-18*(20636.37/waveno)**1.4
!-----------
!  3s2 3p4s 1p
!  elev=40991.884
      if(waveno.lt.elim-elev(22))go to 30
!   glev=3.
      zeff2=16./rydloc*(elim-elev(22))
      x(22)=4.09d-18*(24947.216/waveno)**2
!-----------
!  3s2 3p4s 3p
!     elev=39859.920
      if(waveno.lt.elim-elev(23))go to 30
!     glev=9.
      zeff2=16./rydloc*(elim-elev(23))
      x(23)=1.25d-18*(26079.180/waveno)**2
!--------------
!   si ii 3s2 3p 2p1/2
30    elim=65747.55
!
!-------------------------------------------------------------
!   3s2 3p2 1s
!   elev=15394.370
      if(waveno.lt.elim-elev(24))go to 40
!     glev=1.
      eps=(waveno-70000.)*2./6500.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      reson1=(97.d-18*eps+94.d-18)/(eps**2+1.)
      x(24)=(37.d-18*(50353.180/waveno)**2.40+reson1)/3.
!---------
!  3s2 3p2 1d
!  elev=6298.850
      if(waveno.lt.elim-elev(25))go to 40
!     glev=5.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      eps=(waveno-78600.)*2./13000.
      reson1=(-10.e-18*eps+77.d-18)/(eps**2+1.)
      x(25)=(24.5d-18*(59448.700/waveno)**1.85+reson1)/3.
!--------
!  3s2 3p2 3p2
!  elev=223.157
      if(waveno.lt.elim-elev(26))go to 40
!     glev=5.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      if(waveno.le.74000.) x(26)=72.d-18*(65524.393/waveno)**1.90/3.
      if(waveno.gt.74000.) x(26)=93.d-18*(65524.393/waveno)**4.00/3.
!---------
!  3s2 3p2 3p1
!  elev=77.115
      if(waveno.lt.elim-elev(27))go to 40
!     glev=3.

!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      if(waveno.le.74000.)x(27)=72.d-18*(65524.393/waveno)**1.90*2./3.
      if(waveno.gt.74000.)x(27)=93.d-18*(65524.393/waveno)**4.00*2./3.
!----------
!  3s2 3p2 3p0
!  elev=0.00
      if(waveno.lt.elim-elev(28))go to 40
!     glev=1.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      if(waveno.le.74000.)x(28)=72.d-18*(65524.393/waveno)**1.90/3.
      if(waveno.gt.74000.)x(28)=93.d-18*(65524.393/waveno)**4.00/3.
!
!      si ii 3s2 3p 2p3/2
40    elim=65747.55+287.45
!--------------------------------------------------------------------
!
!----------------
!  3s2 3p2 1s
!  elev=15394.370
      if(waveno.lt.elim-elev(24))go to 50
!     glev=1.
      eps=(waveno-70000.)*2./6500.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      reson1=(97.e-18*eps+94.d-18)/(eps**2+1.)
      x(24)=x(24)+(37.d-18*(50640.630/waveno)**2.40+reson1)*2./3.
!----------------
!  3s2 3p2 1d
!  elev=6298.850
      if(waveno.lt.elim-elev(25))go to 50
!     glev=5.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      eps=(waveno-78600.)*2./13000.
      reson1=(-10.d-18*eps+77.d-18)/(eps**2+1.)
      x(25)=x(25)+(24.5d-18*(59448.700/waveno)**1.85+reson1)*2./3.
!---------------
!  3s2 3p2 3p2
!     elev=223.157
      if(waveno.lt.elim-elev(26))go to 50
!     glev=5.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      if(waveno.le.74000.)x(26)=
     & x(26)+72.d-18*(65524.393/waveno)**1.90*2./3.
      if(waveno.gt.74000.)x(26)=
     & x(26)+93.d-18*(65524.393/waveno)**4.00*2./3.
!-------------------
!
!  3s2 3p2 3p1
!  elev=77.115
      if(waveno.lt.elim-elev(27))go to 50
!     glev=3.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      if(waveno.le.74000.)x(27)=
     & x(27)+72.d-18*(65524.393/waveno)**1.90*2./3.
      if(waveno.gt.74000.)x(27)=
     & x(27)+93.d-18*(65524.393/waveno)**4.00*2./3.
!
!----------------
!  3s2 3p2 3p0
!  elev=0.00
      if(waveno.lt.elim-elev(28))go to 50
!     glev=1.
!     fits to nahar, s.n. and pradhan, a.k. j.phys.b 26, 1109-1127, 1993.
      if(waveno.le.74000.)x(28)=
     & x(28)+72.d-18*(65524.393/waveno)**1.90*2./3.
      if(waveno.gt.74000.)x(28)=
     & x(28)+93.d-18*(65524.393/waveno)**4.00*2./3.
!----------------
!       si II 3s 3p2 4p1/2
50    elim=65747.5+42824.35
!----------------------------------------------------------------------
!
!  3s3p3 1p
!  elev=94000.
      if(waveno.lt.elim-elev(29))go to 60
!     glev=3.
      degen=3.
      zeff2=9./rydloc*(elim-elev(29))
      x(29)=xkarsas(zeff2,3,1)*degen
!------------
!  3s3p3 3s
!     guess
!     elev=79664.0
      if(waveno.lt.elim-elev(30))go to 60
!     glev=3.
      degen=3.
      zeff2=9./rydloc*(elim-elev(30))
      x(30)=xkarsas(zeff2,3,1)*degen
!----------
!  3s3p3 1d
!     guess
!  elev=72000.
      if(waveno.lt.elim-elev(31))go to 60
!     glev=5.
      zeff2=9./rydloc*(elim-elev(31))
      x(31)=xkarsas(zeff2,3,1)*degen
!-------------
!  3s3p3 3p
!  elev=56698.738
      if(waveno.lt.elim-elev(32))go to 60
!     glev=12.
      zeff2=9./rydloc*(elim-elev(32))
      x(32)=xkarsas(zeff2,3,1)*degen
!------------
!  2s2p3 5s
!     elev=33326.053
      if(waveno.lt.elim-elev(33))go to 60
!     glev=5.
      zeff2=9./rydloc*(elim-elev(33))
      x(33)=52.4d-18*(75423.767/waveno)**3.5
!
60    continue
     
      do  l=1,nrhox
         elim=65747.55
!---------------
!  n=5 to infinity
!    gfactor = 6.0d0 
         h1=freq3*6.0d0*2./2./(rydloc*z**2*hckt(l))*
     &     (exp(-max(elim-rydloc*z**2/5.**2,elim-waveno)*hckt(l))
     &     -exp(-elim*hckt(l)))

        do i=1,33
          h1=h1+x(i)*bolt(i,l)
        end do 
        asi1(l)=h1*stim(l)
      end do
      end subroutine 
!--------------- end of main si1op castelli routine -------------!


      double precision function si1seaton3p4p3s()

        implicit none 
!        include 'common.sizebl'
!        include 'common.freqbl'

        integer i    
        double precision x, flam3p4p3s(29), x3p4p3s(29),freqs(29)

        data flam3p4p3s/
     &  2069.28, 2081.70, 2104.23, 2217.88, 2258.00, 2284.87, 2314.16,
     &  2324.52, 2335.20, 2346.21, 2357.57, 2369.29, 2393.88, 2406.80,
     &  2420.15, 2448.24, 2494.23, 2768.43, 2959.12, 3170.60, 3525.02,
     &  4071.66, 4471.89, 4814.83, 5245.71, 5505.27, 5802.79, 6147.16,
     &   6591.82/

        data x3p4p3s/
     &   0.943,   1.129,   1.348,   2.062,   2.416,   2.793,   3.400,
     &   3.429,   2.702,   1.107,   0.265,   0.293,   0.793,   1.004,
     &   1.179,   1.453,   1.764,   3.019,   3.955,   5.120,   7.154,
     &   9.578,  10.390,  10.490,   9.964,   9.350,   8.413,   7.093,
     &   5.228/


        si1seaton3p4p3s=0.0d0 

        do  i=1,29
          freqs(i)=2.9979245d18/flam3p4p3s(i)
        end do 

        if(freq .le. freqs(1)) then 
          i = 2
          do while ((i .lt. 29) .and. (freq .le. freqs(i))) 
            i = i + 1 
          end do 


          x=(freq-freqs(i))/(freqs(i-1)-freqs(i))*
     &       (x3p4p3s(i-1)-x3p4p3s(i))+x3p4p3s(i)
 
            si1seaton3p4p3s=x*1.d-18
        end if
        end function 

!------------------------------------------------------------------!
      double precision  function si1seaton3p4p1d()

        implicit none
!        include 'common.sizebl'
!        include 'common.freqbl'
        integer i  
        double precision x, flam3p4p1d(22), x3p4p1d(22),freqs(22)

        data flam3p4p1d/
     &   1886.54, 1888.90, 1891.36, 1893.92, 1897.28, 1907.68, 1915.27,
     &   1923.68, 1937.57, 1985.57, 2188.16, 2356.14, 2503.61, 2641.37,
     &   2776.86, 3024.33, 3841.03, 4712.79, 5322.58, 5874.98, 6611.74,
     &   7128.99/

        data x3p4p1d/
     &   0.446,   0.430,   0.440,   0.476,   0.553,   0.875,   1.059,
     &   1.181,   1.274,   1.367,   1.688,   2.014,   2.407,   2.897,
     &   3.524,   4.925,   7.984,  10.400,  12.180,  14.080,  17.100,
     &   19.490/

        si1seaton3p4p1d=0.0d0 
    
        do  i=1,22
          freqs(i)=2.99792458d18/flam3p4p1d(i)
        end do 
        
        if( freq .le. freqs(1)) then 
            i = 2 
          do while ((i .lt. 22) .and. (freq .le. freqs(i))) 
            i = i + 1
          end do
 

          x=(freq-freqs(i))/(freqs(i-1)-freqs(i))*
     &       (x3p4p1d(i-1)-x3p4p1d(i))+x3p4p1d(i)

          si1seaton3p4p1d=x*1.d-18
        end if
        end function 

!---------------------------------------------------------------------!

      double precision function si1seaton3p3d3p() 
        implicit  none 
        integer i 
        double precision  x, flam3p3d3p(41), x3p3d3p(41),freqs(41)

        data flam3p3d3p/
     &   2007.57, 2010.35, 2028.40, 2031.61, 2034.88, 2038.22, 2041.61,
     &   2045.06, 2052.21, 2055.88, 2059.63, 2063.45, 2067.35, 2071.33,
     &   2075.39, 2079.54, 2083.78, 2235.67, 2278.14, 2306.66, 2327.12,
     &   2337.82, 2348.85, 2360.23, 2371.98, 2384.11, 2396.63, 2409.58,
     &   2422.96, 2436.80, 2451.12, 2465.94, 2481.30, 2513.73, 2719.51,
     &   2860.71, 3041.18, 3279.06, 4201.09, 5258.93, 6612.72/

        data x3p3d3p/
     &   0.649,   0.656,   0.568,   0.567,   0.625,   0.791,   0.894,
     &   0.830,   0.590,   0.485,   0.461,   0.590,   0.800,   0.962,
     &   1.052,   1.096,   1.117,   1.262,   1.243,   1.154,   1.031,
     &   1.031,   1.166,   1.247,   1.181,   1.027,   0.836,   0.712,
     &   0.823,   1.183,   1.578,   1.871,   2.068,   2.304,   3.167,
     &   3.874,   4.942,   6.521,  11.900,  14.870,  16.440/


        si1seaton3p3d3p=0.0d0 

        do  i=1,41
         freqs(i)=2.99792458d18/flam3p3d3p(i)
        end do 

        if(freq .le. freqs(1)) then 
          i = 2
          do while ((i .lt. 41) .and. (freq .le. freqs(i))) 
            i = i +1
          end do
   
          x=(freq-freqs(i))/(freqs(i-1)-freqs(i))*
     &      (x3p3d3p(i-1)-x3p3d3p(i))+x3p3d3p(i)
          si1seaton3p3d3p=x*1.d-18

        end if
        end function 

!-----------------------------------------------------------------------------------------------
!           End castelli version of SiOP opactities ....                                       !
!-----------------------------------------------------------------------------------------------


      subroutine lukeop
      implicit none
!
!.... SI2, MG2, CA2, N1, O1
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.xnfpbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      integer j
!
!------------------------------- EXECUTION ----------------------------
!


      do j = 1, nrhox
         aluke(j) = stim(j) *
     &      (n1op((j))  * xnfpn(j,1)  + o1op((j))  * xnfpo(j,1)  +
     &       mg2op((j)) * xnfpmg(j,2) + si2op((j)) * xnfpsi(j,2) +
     &       ca2op((j)) * xnfpca(j,2)) / rho(j)
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L U K E O P !!!!!!!!!!!!
!
      double precision function ca2op (j)
      implicit none
!
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  c1218(maxd), c1420(maxd), freq1,  x1044,
     &                  x1218, x1420
      integer itemp1, k
      save    c1218, c1420, x1044, x1218, x1420, freq1, itemp1
!
!------------------------------ INITIALIZATION ------------------------
!
      data freq1, itemp1 / 0.0d0, 0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do k = 1, nrhox
            c1218(k) = 10.0d0 * exp(-1.697 / tkev(k))
            c1420(k) =  6.0d0 * exp(-3.142 / tkev(k))
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1 = freq
         x1044 = 0.0d0
         x1218 = 0.0d0
         x1420 = 0.0d0
!
         if (freq .ge. 2.870454d15) x1044 = 5.4d-20 * 
     &                                      (2.870454d15 / freq) ** 3
         if (freq .ge. 2.460127d15) x1218 = 1.64d-17 * 
     &                                      sqrt(2.460127d15 / freq)
         if (freq .ge. 2.110779d15) x1420 = seaton (2.110779d15, 
     &                                      4.13d-18, 3.0d0, 0.69d0)
      end if
!
      ca2op = x1044 * 2.0d0 + x1218 * c1218(j) + x1420 * c1420(j)
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   C A 2 O P !!!!!!!!!!!!!!
!
      double precision function mg2op (j)
      implicit none
!
!.... CROSS-SECTION TIMES PARTITION FUNCTION
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  c1169(maxd), freq1,  x824, x1169
      integer  itemp1, k
      save     c1169, freq1, itemp1, x824, x1169
!
!------------------------------ INITIALIZATION ------------------------
!
      data freq1, itemp1 / 0.0d0, 0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do k = 1, nrhox
            c1169(k) = 6.0 * exp(-4.43 / tkev(k) )
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1 = freq
         x824  = 0.0d0
         x1169 = 0.0d0
!
         if (freq .ge. 3.635492d15) x824  = seaton (3.635492d15, 
     &                                      1.40d-19, 4.0d0, 6.7d0)
         if (freq .ge. 2.564306d15) x1169 = 5.11d-19 * 
     &                                      (2.564306d15 / freq) ** 3
      end if
!
      mg2op = x824 * 2.0d0 + x1169 * c1169(j)
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   M G 2 O P !!!!!!!!!!!!!!
!
      double precision function n1op (j)
      implicit none
!
!.... CROSS-SECTION TIMES PARTITION FUNCTION
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  c1130(maxd), c1020(maxd), freq1,  
     &                  x853, x1020, x1130
      integer itemp1, k
      save    c1020, c1130, freq1, itemp1, x853, x1020, x1130
!
!------------------------------ INITIALIZATION ------------------------
!
      data freq1, itemp1 / 0.0d0, 0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do k = 1, nrhox
            c1020(k) = 10.0d0 * exp(-2.384 / tkev(k) )
            c1130(k) =  6.0 * exp(-3.575 / tkev(k) )
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1 = freq
         x853  = 0.0d0
         x1020 = 0.0d0
         x1130 = 0.0d0
!
         if (freq .ge. 3.517915d15) x853  = seaton (3.517915d15, 
     &                                      1.142d-17, 2.0d0, 4.29d0)
         if (freq .ge. 2.941534d15) x1020 = seaton (2.941534d15, 
     &                                      4.41d-18, 1.5d0, 3.85d0)
         if (freq .ge. 2.653317d15) x1130 = seaton (2.653317d15, 
     &                                      4.2d-18, 1.5d0, 4.34d0)
      end if
!
      n1op = x853 * 4.0 + x1020 * c1020(j) + x1130 * c1130(j)
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   N 1 O P !!!!!!!!!!!!!!!!
!
      double precision function o1op ( j )
      implicit none
!
!.... FROM DEANE PETERSON AFTER PEACH
!.... CROSS-SECTION TIMES PARTITION FUNCTION
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  freq1,  x911
      save     freq1, x911
!
!
!------------------------------ INITIALIZATION ------------------------
!
      data freq1 / 0.0d0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (freq .ne. freq1) then
         freq1 = freq
         x911  = 0.0d0
         if (freq .ge. 3.28805d15) x911 = seaton (3.28805d15, 2.94d-18, 
     &                                            1.0d0, 2.66d0)
      end if
!
      o1op = x911 * 9.0
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   O 1 O P !!!!!!!!!!!!!!!!
!
      double precision function si2op (j)
      implicit none
!
!.... CROSS-SECTION TIMES THE PARTITION FUNCTION
!
!.... 1995 OCT - TEST n .le. 7 BEFORE freq .le. freqsi TO AVOID AN ARRAY
!                OUT OF BOUNDS
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer j
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  d, d1, dt(maxd), flog(9),  freq1, freqsi(8), 
     &                  peach(6, 14), tlg(6), x(6)
      integer  it, itemp1, k, n, nt(maxd)
      save     dt, flog, freq1, freqsi, itemp1, nt, peach, tlg, x
!
!------------------------------ INITIALIZATION ------------------------
!
      data itemp1 / 0 /
!
!....                2P,       2D,       2P,       2D,       2P
      data flog / 36.32984, 36.14752, 35.91165, 34.99216, 34.95561,
     &            34.45951, 34.36234, 34.27572, 34.20161 /
!
      data freqsi /    4.9965417d15, 3.9466738d15, 1.5736321d15,
     &   1.5171539d15, 9.2378947d14, 8.3825004d14, 7.6869872d14 , 0.0/
!
      data peach /
!....       10000     12000     14000     16000     18000     20000     WAVE(A)
     &    -43.8941, -43.8941, -43.8941, -43.8941, -43.8941, -43.8941,     500
     &    -42.2444, -42.2444, -42.2444, -42.2444, -42.2444, -42.2444,     600
     &    -40.6054, -40.6054, -40.6054, -40.6054, -40.6054, -40.6054,     759
     &    -54.2389, -52.2906, -50.8799, -49.8033, -48.9485, -48.2490,     760
     &    -50.4108, -48.4892, -47.1090, -46.0672, -45.2510, -44.5933,    1905
     &    -52.0936, -50.0741, -48.5999, -47.4676, -46.5649, -45.8246,    1906
     &    -51.9548, -49.9371, -48.4647, -47.3340, -46.4333, -45.6947,    1975
     &    -54.2407, -51.7319, -49.9178, -48.5395, -47.4529, -46.5709,    1976
     &    -52.7355, -50.2218, -48.4059, -47.0267, -45.9402, -45.0592,    3245
     &    -53.5387, -50.9189, -49.0200, -47.5750, -46.4341, -45.5082,    3246
     &    -53.2417, -50.6234, -48.7252, -47.2810, -46.1410, -45.2153,    3576
     &    -53.5097, -50.8535, -48.9263, -47.4586, -46.2994, -45.3581,    3577
     &    -54.0561, -51.2365, -49.1980, -47.6497, -46.4302, -45.4414,    3900
     &    -53.8469, -51.0256, -48.9860, -47.4368, -46.2162, -45.2266 /   4200
!
      data tlg / 9.21034, 9.39266, 9.54681, 9.68034, 9.79813, 9.90349 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         freq1  = 0.0d0
         itemp1 = itemp
!
         do k = 1, nrhox
            n     = max( min(5, int(t(k) / 2000.0d0) -4), 1)
            nt(k) = n
            dt(k) = (tlog(k) - tlg(n)) / (tlg(n + 1) - tlg(n))
         end do
!
      end if
!
      if (freq .ne. freq1) then
         freq1   = freq
         n       = 1
!
!.... THIS LOOP RETURNS N <= 8
!
         do while (n .le. 7 .and. freq .le. freqsi(n))
           n=n+1
         end do
!
         d  = (freqln - flog(n)) / (flog(n + 1) - flog(n))
         d1 = 1.0 - d
         if (n .gt. 2)  n = 2 * n - 2
         if (n .eq. 14) n = 13
!
         do it = 1, 6
            x(it) = peach(it, n + 1) * d + peach(it, n) * d1
         end do
!
      end if
!
      n     = nt(j)
      si2op = exp(x(n) * (1.0 - dt(j)) + x(n + 1) * dt(j) ) * 6.0
      end function 
!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   S I 2 O P !!!!!!!!!!!!!!
!
      subroutine hotop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
!      include 'common.freqbl'
      include 'common.opsblk'
      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xnfblk'
      include 'common.xnfpbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision  a(7, 60), ac2op(maxd), ac3op(maxd), ac4op(maxd), 
     &                  an2op(maxd), an3op(maxd), an4op(maxd), 
     &                  an5op(maxd), ao2op(maxd), ao3op(maxd), 
     &                  ao4op(maxd), ao5op(maxd), ao6op(maxd),  
     &                  free, ratio, xnfp(maxd, 21), xsect, xx
      integer  i, id, itemp1, j, nuend, num
      logical  start
      save     a, itemp1, start, xnfp
!
!-------------------------------- EXTERNAL -----------------------------
!
!      external c2op, c3op, c4op,  n2op, n3op, n4op, n5op, 
!     &         o2op, o3op, o4op, o5op, o6op
!
!------------------------------ INITIALIZATION ------------------------
!
      data start / .true. /
!
      data itemp1, num / 0, 60 /
!
      data ((a(i,j), i = 1, 7), j = 1, 15) /
!           FREQUENCY     ALPHA       A      P      G      E      ID       CODE
     &     4.149945d15,  6.90d-18,  1.000,  6.0,   6.0,  13.71,   2.0,     6.01
     &     4.574341d15,  2.50d-18,  1.000,  4.0,   2.0,  11.96,   2.0,     6.01
     &     5.220770d15,  1.08d-17,  1.000,  4.0,  10.0,   9.28,   2.0,     6.01
     &     5.222307d15,  5.35d-18,  3.769,  2.0,   1.0,   0.00,  16.0,    10.00
     &     5.892577d15,  4.60d-18,  1.950,  6.0,   6.0,   0.00,   2.0,     6.01
     &     6.177022d15,  3.50d-18,  1.000,  4.0,  12.0,   5.33,   2.0,     6.01
     &     6.181062d15,  6.75d-18,  3.101,  5.0,   1.0,   4.05,   6.0,     7.01
     &     6.701879d15,  6.65d-18,  2.789,  5.0,   5.0,   1.90,   6.0,     7.01
     &     7.158382d15,  6.65d-18,  2.860,  6.0,   9.0,   0.00,   6.0,     7.01
     &     7.284488d15,  3.43d-18,  4.174,  5.0,   6.0,   5.02,  11.0,     8.01
     &     7.693612d15,  3.53d-18,  3.808,  5.0,  10.0,   3.33,  11.0,     8.01
     &     7.885955d15,  2.32d-18,  3.110,  5.0,   6.0,   5.02,  11.0,     8.01
     &     8.295079d15,  3.97d-18,  3.033,  5.0,  10.0,   3.33,  11.0,     8.01
     &     8.497686d15,  7.32d-18,  3.837,  5.0,   4.0,   0.00,  11.0,     8.01
     &     8.509966d15,  2.00d-18,  1.750,  7.0,   3.0,  12.69,   3.0/     6.02
!
      data ((a(i,j), i = 1, 7), j = 16, 30) /
!           FREQUENCY     ALPHA       A      P      G      E      ID       CODE
     &     8.572854d15,  1.68d-18,  3.751,  5.0,   6.0,   5.02,  11.0,     8.01
     &     9.906370d15,  4.16d-18,  2.717,  3.0,   6.0,   0.00,  17.0,    10.01
     &     1.000693d16,  2.40d-18,  1.750,  7.0,   9.0,   6.50,   3.0,     6.02
     &     1.046078d16,  4.80d-18,  1.000,  4.0,  10.0,  12.53,   7.0,     7.02
     &     1.067157d16,  2.71d-18,  2.148,  3.0,   6.0,   0.00,  17.0,    10.01
     &     1.146734d16,  2.06d-18,  1.626,  6.0,   6.0,   0.00,   7.0,     7.02
     &     1.156813d16,  5.20d-19,  2.126,  3.0,   6.0,   0.00,  17.0,    10.01
     &     1.157840d16,  9.10d-19,  4.750,  4.0,   1.0,   0.00,   3.0,     6.02
     &     1.177220d16,  5.30d-18,  1.000,  4.0,  12.0,   7.10,   7.0,     7.02
     &     1.198813d16,  3.97d-18,  2.780,  6.0,   1.0,   5.35,  12.0,     8.02
     &     1.267503d16,  3.79d-18,  2.777,  6.0,   5.0,   2.51,  12.0,     8.02
     &     1.327649d16,  3.65d-18,  2.014,  6.0,   9.0,   0.00,  12.0,     8.02
     &     1.361466d16,  7.00d-18,  1.000,  2.0,   5.0,   7.48,  12.0,     8.02
     &     1.365932d16,  9.30d-19,  1.500,  7.0,   6.0,   8.00,   4.0,     6.03
     &     1.481487d16,  1.10d-18,  1.750,  7.0,   3.0,  16.20,   8.0 /    7.03
!
      data ((a(i,j), i = 1, 7), j = 31, 45) /
!           FREQUENCY     ALPHA       A      P      G      E      ID       CODE
     &     1.490032d16,  5.49d-18,  3.000,  5.0,   1.0,   6.91,  18.0,    10.02
     &     1.533389d16,  1.80d-18,  2.277,  4.0,   9.0,   0.00,  18.0,    10.02
     &     1.559452d16,  8.70d-19,  3.000,  6.0,   2.0,   0.00,   4.0,     6.03
     &     1.579688d16,  4.17d-18,  2.074,  4.0,   5.0,   3.20,  18.0,    10.02
     &     1.643205d16,  1.39d-18,  2.792,  5.0,   5.0,   3.20,  18.0,    10.02
     &     1.656208d16,  2.50d-18,  2.346,  5.0,   9.0,   0.00,  18.0,    10.02
     &     1.671401d16,  1.30d-18,  1.750,  7.0,   9.0,   8.35,   8.0,     7.03
     &     1.719725d16,  1.48d-18,  2.225,  5.0,   9.0,   0.00,  18.0,    10.02
     &     1.737839d16,  2.70d-18,  1.000,  4.0,  10.0,  15.74,  13.0,     8.03
     &     1.871079d16,  1.27d-18,  0.831,  6.0,   6.0,   0.00,  13.0,     8.03
     &     1.873298d16,  9.10d-19,  3.000,  4.0,   1.0,   0.00,   8.0,     7.03
     &     1.903597d16,  2.90d-18,  1.000,  4.0,  12.0,   8.88,  13.0,     8.03
     &     2.060738d16,  4.60d-18,  1.000,  3.0,  12.0,  22.84,  19.0,    10.03
     &     2.125492d16,  5.90d-19,  1.000,  6.0,   6.0,   9.99,   9.0,     7.04
     &     2.162610d16,  1.69d-18,  1.937,  5.0,   6.0,   7.71,  19.0 /   10.03
!
      data ((a(i,j), i = 1, 7), j = 46, 60) /
!           FREQUENCY     ALPHA       A      P      G      E      ID       CODE
     &     2.226127d16,  1.69d-18,  1.841,  5.0,  10.0,   5.08,  19.0,    10.03
     &     2.251163d16,  9.30d-19,  2.455,  6.0,   6.0,   7.71,  19.0,    10.03
     &     2.278001d16,  7.90d-19,  1.000,  6.0,   9.0,  10.20,  14.0,     8.04
     &     2.317678d16,  1.65d-18,  2.277,  6.0,  10.0,   5.08,  19.0,    10.03
     &     2.348946d16,  3.11d-18,  1.963,  6.0,   4.0,   0.00,  19.0,    10.03
     &     2.351911d16,  7.30d-19,  1.486,  5.0,   6.0,   7.71,  19.0,    10.03
     &     2.366973d16,  5.00d-19,  1.000,  4.0,   2.0,   0.00,   9.0,     7.04
     &     2.507544d16,  6.90d-19,  1.000,  6.0,   3.0,  19.69,  14.0,     8.04
     &     2.754065d16,  7.60d-19,  1.000,  2.0,   1.0,   0.00,  14.0,     8.04
     &     2.864850d16,  1.54d-18,  2.104,  6.0,   1.0,   7.92,  20.0,    10.04
     &     2.965598d16,  1.53d-18,  2.021,  6.0,   5.0,   3.76,  20.0,    10.04
     &     3.054151d16,  1.40d-18,  1.471,  6.0,   9.0,   0.00,  20.0,    10.04
     &     3.085141d16,  2.80d-18,  1.000,  4.0,   5.0,  11.01,  20.0,    10.04
     &     3.339687d16,  3.60d-19,  1.000,  6.0,   2.0,   0.00,  15.0,     8.05
     &     3.818757d16,  4.90d-19,  1.145,  6.0,   6.0,   0.00,  21.0 /   10.05
!
!------------------------------- EXECUTION ----------------------------
!
!.... COPY PARTS OF XNFPC/N/O INTO XNFP FOR USE WITH ID AND AN INDEX
!
      if (start) then
         start = .false.
!.... CARBON
!
         do i = 1, 4
!
            do j = 1, nrhox
               xnfp (j, i) = xnfpc (j, i)
            end do
!
         end do
!
!.... NITROGEN
!
         do i = 1, 5
!
            do j = 1, nrhox
               xnfp (j, i + 4) = xnfpn (j, i)
            end do
!
         end do
!
!.... OXYGEN
!
         do i = 1, 6
!
            do j = 1, nrhox
               xnfp (j, i + 9) = xnfpo (j, i)
            end do
!
         end do
!
!.... NEON
!
         do i = 1, 6
!
            do j = 1, nrhox
               xnfp (j, i + 15) = xnfpne (j, i)
            end do
!
         end do
!
      end if
!
!.... FREE-FREE
!
      do j = 1, nrhox
         free = coulff(j,1) * 1.0 ** 2 * (xnfc(j,2) + xnfn(j,2) + 
     &             xnfo(j,2) + xnfne(j,2) + xnfmg(j,2) + xnfsi(j,2) + 
     &             xnfs(j,2) + xnffe(j,2)) + 
     &          coulff(j,2) * 2.0 ** 2 * (xnfc(j,3) + xnfn(j,3) + 
     &             xnfo(j,3) + xnfne(j,3) + xnfmg(j,3) + xnfsi(j,3) + 
     &             xnfs(j,3) + xnffe(j,3)) + 
     &          coulff(j,3) * 3.0 ** 2 * (xnfc(j,4) + xnfn(j,4) + 
     &             xnfo(j,4) + xnfne(j,4) + xnfmg(j,4) + xnfsi(j,4) + 
     &             xnfs(j,4) + xnffe(j,4)) + 
     &          coulff(j,4) * 4.0 ** 2 * (xnfc(j,5) + xnfn(j,5) + 
     &             xnfo(j,5) + xnfne(j,5) + xnfmg(j,5) + xnfsi(j,5) + 
     &             xnfs(j,5) + xnffe(j,5)) + 
     &          coulff(j,5) * 5.0 ** 2 *(xnfc(j,6) + xnfn(j,6) + 
     &             xnfo(j,6) + xnfne(j,6) + xnfmg(j,6) + xnfsi(j,6) + 
     &             xnfs(j,6))
         ahot(j) = free * 3.6919d8 / freq * xne(j) /freq / freq / 
     &             sqrt(t(j))
      end do
!
      do j = 1, nrhox
         ac2op(j) = 0.0d0
         ac3op(j) = 0.0d0
         ac4op(j) = 0.0d0
         an2op(j) = 0.0d0
         an3op(j) = 0.0d0
         an4op(j) = 0.0d0
         an5op(j) = 0.0d0
         ao2op(j) = 0.0d0
         ao3op(j) = 0.0d0
         ao4op(j) = 0.0d0
         ao5op(j) = 0.0d0
         ao6op(j) = 0.0d0
      end do
!
      call c2op (ac2op)
      call c3op (ac3op)
      call c4op (ac4op)
      call n2op (an2op)
      call n3op (an3op)
      call n4op (an4op)
      call n5op (an5op)
      call o2op (ao2op)
      call o3op (ao3op)
      call o4op (ao4op)
      call o5op (ao5op)
      call o6op (ao6op)
!
      do j = 1, nrhox
         ahot(j) =  ahot(j) + ac2op(j) + ac3op(j) + ac4op(j) + 
     &             an2op(j) + an3op(j) + an4op(j) + an5op(j) +
     &             ao2op(j) + ao3op(j) + ao4op(j) + ao5op(j) + 
     &             ao6op(j)
      end do
!
      nuend = 1
!
      do while (freq .ge. a(1, nuend) .and. nuend .lt. 60)
         nuend = nuend + 1
      end do
!
      nuend = nuend - 1
!
      do i = 1, nuend
         ratio = a(1, i) / freq
         xsect = a(2, i) * (a(3, i) + ratio - a(3, i) * ratio) *
     &           sqrt(ratio ** int( a(4, i)))
         id = a (7, i)
!
         do j = 1, nrhox
            xx = xsect * xnfp(j, id) * a (5, i)
            if (xx .gt. 0.01 * ahot(j)) ahot(j) = ahot(j) + xx /
     &               exp(a(6, i) / tkev(j) )
         end do
!
      end do
!
      do j = 1, nrhox
         ahot(j) = ahot(j) * stim(j) / rho(j)
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H O T O P !!!!!!!!!!!!!!
!
      subroutine c2op(ac2op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ac2op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   C 2 O P !!!!!!!!!!!!!!!!
!
      subroutine c3op(ac3op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ac3op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   C 3 O P !!!!!!!!!!!!!!!!
!
      subroutine c4op(ac4op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ac4op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   C 4 O P !!!!!!!!!!!!!!!!
!
      subroutine n2op(an2op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision an2op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   N 2 O P !!!!!!!!!!!!!!!!
!
      subroutine n3op(an3op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision an3op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   N 3 O P !!!!!!!!!!!!!!!!
!
      subroutine n4op(an4op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision an4op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   N 4 O P !!!!!!!!!!!!!!!!
!
      subroutine n5op(an5op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision an5op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   N 5 O P !!!!!!!!!!!!!!!!
!
      subroutine o2op(ao2op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ao2op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   O 2 O P !!!!!!!!!!!!!!!!
!
      subroutine o3op(ao3op)
      implicit none
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ao3op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   O 3 O P !!!!!!!!!!!!!!!!
!
      subroutine o4op(ao4op)
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ao4op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   O 4 O P !!!!!!!!!!!!!!!!
!
      subroutine o5op(ao5op)
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ao5op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   O 5 O P !!!!!!!!!!!!!!!!
!
      subroutine o6op(ao6op)
!
!-------------------------------- DUMMY VARIABLES ----------------------
!
      double precision ao6op(1)
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   O 6 O P !!!!!!!!!!!!!!!!
!
      subroutine elecop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.opsblk'
      include 'common.stateb'
      logical cast_version 
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      integer j
!
!------------------------------- EXECUTION ----------------------------
!
      cast_version = .false. 
      if (cast_version) then 
        do j = 1, nrhox
         sigel(j) = 0.6653d-24 * xne(j) / rho(j)
        end do
      else

        do j = 1, nrhox
           sigel(j) = sige * xne(j) / rho(j)
        end do
      end if
!!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   E L E C O P !!!!!!!!!!!!
!
      subroutine h2raop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.depart'
      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!---------------------------- LOCAL VARIABLES  ------------------------
!
      double precision sig, wave, w4, ww, xnh2(maxd)
      integer itemp1, j
      save itemp1, xnh2
!
!------------------------------ INITIALIZATION ------------------------
!
      data itemp1 / 0 /
!
!------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
!
         do j = 1, nrhox
            xnh2(j) = (xnfph(j, 1) * 2.0 * bhyd(j, 1) ) ** 2 *
     &                exp(4.477 / tkev(j) - 4.6628d1 + 
!!!!     &                (1.8031d-3 + (-5.0739d-7 + (8.1424d-11 - 
     &                (1.8031d-3 + (-5.0239d-7 + (8.1424d-11 - 
     &                 5.0501d-15 * t(j) ) * t(j) ) * t(j) ) * t(j) -
     &                 1.5 * tlog(j) ) / rho(j)
         end do
!
      end if
!
      wave = 2.997925d18 / min (freq, 2.922d15)
      ww = 1.0 / (wave ** 2 )
      w4 = ww * ww
      sig = (8.14d-13 + 1.28d-6 * ww + 1.61 * w4 ) * w4
!
      do j = 1, nrhox
         sigh2(j) = sig * xnh2(j)
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H 2 R A O P !!!!!!!!!!!!
!
      subroutine hlinop
      implicit none
!
!.... REQUIRES STARK AND COULX
!
!--------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.depart'
!      include 'common.freqbl'
      include 'common.ionsbl'
      include 'common.opsblk'
      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!------------------------------ LOCAL VARIABLES ------------------------
!
      double precision  a, bhydjm, bolt(maxd, 4),  h, mlast(maxd), 
     &                  s
      integer itemp1, j, m, m1, m2, mfreq, n
      logical more39
      save itemp1, bolt, mlast
!
!
!------------------------------- INITIALIZATION ------------------------
!
      data itemp1 / 0 /
!
!-------------------------------- EXECUTION ----------------------------
!
      if (itemp .ne. itemp1) then
         itemp1=itemp
!
         do j = 1, nrhox
            mlast(j) = 1100.0d0 / xne(j) ** 0.133333333
!
            do n = 1, 4
               bolt(j,n) = exp(-(13.595 - 13.595 / dble(n*n)) / 
     &                     tkev(j)) * 2.0 * dble(n*n) * bhyd(j,n) * 
     &                     xnfph(j,1) / rho(j)
            end do
!
         end do
!
      end if
!
      n = sqrt(3.28805d15 / freq)
!
      if ( (n .eq. 1 .and. freq .ge. 2.00d15) .or. 
     &     (n .eq. 2 .and. freq .ge. 4.44d14) .or.
     &     (n .eq. 3)                         .or.
     &     (n .eq. 4)                        ) then
         mfreq = sqrt(3.28805d15 / (3.28805d15 / dble(n*n) - freq))
!
         do j = 1, nrhox
            m1 = mfreq
            m2 = m1 + 1
            m1 = max0(m1, n + 1)
            h = 0.0d0
            s = 0.0d0
            more39 = .false.
!
            if (m1 .le. 6) then
               more39 = .true.
!
            else if (m1 .gt. mlast(j)) then
               ahline(j) = coulx(n, 3.28806d15 / dble(n*n), 1.0d0) * 
     &                     (1.0 - ehvkt(j) / bhyd(j,n)) * bolt(j,n)
               shline(j) = bnu(j) * stim(j) / (bhyd(j,n) - ehvkt(j))
!
            else
               m1 = m1 - 1
               m2 = m2 + 3
!
               if (n .lt. 4 .or. m1 .gt. 8) then
                  more39 = .true.
!
               else
                  more39 = .true.
                  h = stark(3,4,j) * (1.0 - ehvkt(j) * bhyd(j,4) / 
     &                bhyd(j,3)) * bolt(j,3)
                  s = h * bnu(j) * stim(j) / 
     &                (bhyd(j,3) / bhyd(j,4) - ehvkt(j))
               end if
!
            end if
!
            if (more39) then
!
               do m = m1, m2
                  bhydjm = 1.0
                  if (m .le. 6) bhydjm = bhyd(j,m)
!
!.... ASSUMING FREQ APROXIMATELY FREQNM
!
                  a = stark(n,m,j) * (1.0 - ehvkt(j) * bhydjm / 
     &                bhyd(j,n)) * bolt(j,n)
                  h = h + a
                  s = s + a * bnu(j) * stim(j) / (bhyd(j,n) / 
     &                bhydjm - ehvkt(j))
               end do
!
               ahline(j) = h
               shline(j) = s / h
            end if
!
         end do
!
      end if
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   H L I N O P !!!!!!!!!!!!
!
      double precision function stark (n, m, j)
      implicit none
!
!------------------------------- DUMMY VARIABLES -----------------------
!
      integer j, m, n
!
!--------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!
!      include 'common.freqbl'
      !include 'common.rhoxbl'
      include 'common.stateb'
      include 'common.tempbl'
!
!------------------------------ LOCAL VARIABLES ------------------------
!
      double precision  beta, dbeta, del, dioi, exint, exy2, f0(maxd), 
     &                  fnm, freqnm, fstark(10, 4), impact, knm, 
     &                  knmtab(5, 4), nn, mm, prof, qstat, ratio, 
     &                  x, xm, xn, xx, y1, y2
      integer itemp1, l, mminn
      save itemp1, f0
!
!------------------------------- INITIALIZATION ------------------------
!
      data itemp1 / 0 /
!
      data knmtab / 0.000356, 0.000523, 0.00109, 0.00149, 0.00225,
     &              0.0125,   0.0177,   0.028,   0.0348,  0.0493,
     &              0.124,    0.171,    0.223,   0.261,   0.342,
     &              0.683,    0.866,    1.02,    1.19,    1.46   /
!
      data fstark / 0.1387,   0.07910,  0.02126,  0.01394,  0.006462, 
     &              0.004814, 0.002779, 0.002216, 0.001443, 0.001201, 
     &              0.3921,   0.1193,   0.03766,  0.02209,  0.01139,
     &              0.008036, 0.005007, 0.003850, 0.002658, 0.002151, 
     &              0.6103,   0.1506,   0.04931,  0.02768,  0.01485, 
     &              0.01023,  0.006588, 0.004996, 0.003542, 0.002838, 
     &              0.8163,   0.1788,   0.05985,  0.03189,  0.01762, 
     &              0.01196,  0.007825, 0.005882, 0.004233, 0.003375 /
!
!--------------------------- STATEMENT FUNCTION ------------------------
!
      exint(x) = -log(x) - 0.57516d0 + (0.97996d0 - (0.21654d0 - 
     &            (0.033572d0 - (0.0029222d0 - 10.05439d-4 * x) * x) * 
     &            x) * x) * x
!
!----------------------------- EXECUTION -------------------------------
!
      if (itemp .ne. itemp1) then
         itemp1=itemp
!
         do l = 1, nrhox
            f0(l) = 1.25d-9 * xne(l) ** 0.6666667
         end do
!
      end if
!
      xn = n
      xm = m
      x = xn / xm
      xx = x ** 2
      nn = n * n
      mm = m * m
      mminn = m - n
!
      if (mminn .gt. 5) then
         knm = 5.5d-5 * (nn * mm) ** 2 / (mm - nn)
!
      else
         knm = knmtab(mminn, n)
      end if
!
      if (mminn .gt. 10) then
         fnm = fstark(10, n) * ((20.0 * xn + 100.0) / (xn + 10.0) / 
     &         xm / (1.0 - xx)) ** 3
!
      else
         fnm = fstark(mminn, n)
      end if
!
      freqnm = ryd * (1.0 / nn - 1.0 / mm)
      del    = abs(freq - freqnm)
      dbeta  = 2.997925d18 / freqnm ** 2 / f0(j) / knm
      beta   = dbeta * del
      y1     = mm * del * hkt(j)* 0.5
      y2     = (3.14159 * 3.14159 / 2.0 / 0.0265384 / 2.997925d10) * 
     &          del ** 2 / xne(j)
      qstat  = 1.5 + 0.5 * (y1 ** 2 - 1.384) / (y1 ** 2 + 1.384)
      impact = 0.0d0
!
      if (y1 .le. 8.0 .and. y1 .ge. y2) then
         exy2 = 0.0d0
         if(y2 .le. 8.0) exy2 = exint(y2)
         impact = 1.438 * sqrt(y1 * (1.0 - xx)) * (0.4 * exp(-y1) + 
     &            exint(y1) - 0.5 * exy2)
      end if
!
      if (beta .gt. 20.0d0) then
         prof = 1.5 / beta / beta / sqrt(beta)
         dioi = 6.28 * 1.48d-25 * (2.0 * mm * ryd / del) * xne(j) * 
     &          (sqrt(2.0 * mm * ryd / del) * 
     &          (1.3 * qstat + 0.30 * impact) - 3.9 * ryd * hkt(j))
         ratio = qstat * min(1.0d0 + dioi, 1.25d0) + impact
!
      else
         prof = 8.0 / (80.0d0 + beta ** 3)
         ratio = qstat + impact
      end if
!
      stark = 0.0265384 * fnm * prof * dbeta * ratio
      end function 
!
!!!!!!!!!!!!! E N D   O F   F U N C T I O N   S T A R K !!!!!!!!!!!!!!

      subroutine interPT(n, stepwt)
      implicit none
!
!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE OPACITY FILE IS GIVEN
!.... ONLY FOR THAT VTURB
!.... ODF were already read into kapw which is a global variable from
!.... the ODF.nc file
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      integer, intent(in) :: n
      integer  n1
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.turbpr'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision, intent(out) :: stepwt

      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd),
     &                  plog, tl, wave, x, y


      integer    j, ip, ipj(maxd), it, itemp1, itj(maxd), iwl 

      save       co1, co2, co3, co4,  ipj, itj, iwl 
      
!
!----------------------------- INITIALIZATION --------------------------
!
      data itemp1 / 0 /
!
!---------------------------------- EXECUTION --------------------------
!
! - iwl is the loop index over all ibin s
!
! isubbin - nsteps
! wlend it the freq where the bin finishes, so inifreset(2:ibin+1)  
! global binwith variable, sbwt,  and globale boarders of subbins sbvalues

     

      if (itemp .ne. itemp1) then
         itemp1 = itemp
         iwl    = 1
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
!
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         end do
!
      end if
!
      wave = 2.997925d17 / freq
!
      do while (inifreset(iwl+1) .le. wave)
         iwl = iwl + 1
      end do


!    we start with the smallest sub-bin
      n1= isubbin+1 -n
      stepwt = sbwt(n1,iwl) 
 
! but now there are globally known!  (sbwt(n, iwl) = wt(n), nsteps = isubbin)
!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         a  = exp( co1(j) * dble(kapw(n1, iwl, ip-1, it-1)) +
     &             co2(j) * dble(kapw(n1, iwl, ip  , it-1)) +
     &             co3(j) * dble(kapw(n1, iwl, ip-1, it)) +
     &             co4(j) * dble(kapw(n1, iwl, ip  , it)) )
         alines(j) = a
      end do
!
      end subroutine 



!----------------------------------------------------------------------!
      subroutine linopb (n, nsteps, stepwt)
      implicit none
!
!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE OPACITY FILE IS GIVEN
!.... ONLY FOR THAT VTURB
!.... This one is for the xxxbig0 etc files. call linopl for the xxxlit0 
!.... files
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision  stepwt
      integer n, nsteps
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.turbpr'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd),
     &                  plog, tl, wave, 
     &                  wlend(328), wlb, wle, x, y
      integer(kind=2)   kap(25, 57, 12), kapwo(25, 57, 12, 328)
      integer    ip, ipj(maxd), it, itemp1, itj(maxd), iwl, 
     &           j, n12, istep, ios, nsize 
      logical    first
      save       co1, co2, co3, co4, first, ipj, itj, iwl,
     &           kapwo,  wlend
!
!----------------------------- INITIALIZATION --------------------------
!
      data first / .true. /
!
      data itemp1 / 0 /
!
!
!.... VALUES ARE 1/60, 2/60, 3/60, 6/60
!
!      data wt / 1.66666667d-2, 3.33333333d-2, 5.00000000d-2,
!     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1,
!     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1, 
!     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1 /



      if (first) then 

      first = .false.


         call get_bin_size(nsize, wlend, 328) ! 
!
       do iwl = 1, nsize 
             do it = 1,numt
              read(19)((kap(ip,it,istep),istep=1,isubbin),ip=1,numpres)
             end do
        
              do ip = 1, numpres
                do it = 1, numt 
                  do istep = 1, isubbin
                    kapwo(ip,it,istep,iwl) = kap(ip,it,istep)
                  end do
                end do
              end do
         end do
!
         print *, ' finished reading kapwo '

         iwl = iwl - 1
         wlend(iwl) = 1000000.0d0
         iwl = 1
      end if
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
         iwl    = 1
!
         do j = 1, nrhox
            tl = min(max(tlog(j) / tenlog,tabt(1)),tabt(numt))
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. (numt))
               it = it + 1
            end do
!
            plog = min(max(log10(p(j)), tabp(1)),tabp(numpres))
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. (numpres))
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
!
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         end do
!
      end if
!
      wave = 2.997925d17 / freq
!
      do while (wlend(iwl) .le. wave)
         iwl = iwl + 1
      end do

!
      nsteps = isubbin 
      n12    = isubbin+1 - n
      stepwt = sbwt(n,iwl)
!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         a  = exp( co1(j) * dble(kapwo(ip-1, it-1, n12, iwl)) +
     &             co2(j) * dble(kapwo(ip  , it-1, n12, iwl)) +
     &             co3(j) * dble(kapwo(ip-1, it  , n12, iwl)) +
     &             co4(j) * dble(kapwo(ip  , it  , n12, iwl)) )
         alines(j) = a
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L I N O P B !!!!!!!!!!!!
!
      subroutine linopl (n, nsteps, stepwt)

      implicit none
!
!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE OPACITY FILE IS GIVEN
!.... ONLY FOR THAT VTURB
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision  stepwt
      integer n, nsteps
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.turbpr'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd),
     &                  plog,  tl, wave,
     &                  wlend(1212), x, y
      integer(kind=2)  kap(25, 57, 12), kapwo(25, 57, 12, 1212)
      integer    ip, ipj(maxd), it, itemp1, itj(maxd), iwl,
     &           iwllast, j, n12, istep,  nsize 
      logical    first, cast_version 
      save       co1, co2, co3, co4, first, ipj, itj, iwl, iwllast,
     &           kapwo,  wlend
#ifdef DYNBIN 
      integer i 
      double precision binbound(1212, 12)   
      integer nsubbin(1212)
      save binbound, nsubbin 
#endif
!
!----------------------------- INITIALIZATION --------------------------
!
      data first / .true. /
!
      data itemp1 / 0 /
!
!
!---------------------------------- EXECUTION --------------------------
!
      cast_version = .true.

      if (first) then
         first = .false.
!
         call get_bin_size(nsize, wlend, 1212 ) ! 1 for lit and 0 for big! 

         do iwl = 1, nsize 
             do it = 1,numt 
              read(19) ((kap(ip,it,istep),istep=1,isubbin),ip=1,numpres)
             end do 
!            swap do loops ip and it
             do it = 1, numt 
                do ip = 1, numpres 
                   do istep = 1, isubbin 
                      kapwo(ip,it,istep,iwl) = kap(ip,it,istep)
                   end do
                end do
             end do
         end do

#ifdef DYNBIN 
      open (unit = 78, type = 'old', file='subbin_info.dat',  
     & form = 'formatted',    shared, readonly )

       do i = 1, nsize 
         read(78,*) nsubbin(i)
         read(78,*) (binbound(i,j), j=1, nsubbin(i))
       end do 

       close(unit=78) 

#endif


!
        print *, ' finished reading kapwo '

         iwl = iwl - 1
         wlend(iwl) = 1000000.0d0
         iwl = 1
      end if
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
         iwl    = 1
!
         do j = 1, nrhox
              tl = min(max(tlog(j) / tenlog, tabt(1)), tabt(numt))
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. numt)
               it = it + 1
            end do
              plog = min(max(log10(p(j)),tabp(1)), tabp(numpres))
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. numpres )
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
!
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         end do
!
      end if
!
      wave = 2.997925d17 / freq
!
      do while (wlend(iwl) .le. wave)
         iwl = iwl + 1
      end do

!
#ifdef DYNBIN
      nsteps = nsubbin(iwl)
      n12    = nsubbin(iwl)+1-n
      if (n12 .gt. 1) then
      stepwt = binbound(iwl,n12)-binbound(iwl,n12-1) 
      else
      stepwt = binbound(iwl,n12)
      endif

#else
      nsteps = isubbin
      n12    = isubbin+1 - n
      stepwt = sbwt(n, iwl)
#endif


!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         a  = exp( co1(j) * dble(kapwo(ip-1, it-1, n12, iwl)) +
     &             co2(j) * dble(kapwo(ip  , it-1, n12, iwl)) +
     &             co3(j) * dble(kapwo(ip-1, it  , n12, iwl)) +
     &             co4(j) * dble(kapwo(ip  , it  , n12, iwl)) )

         alines(j) = a
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L I N O P L !!!!!!!!!!!!
!
      subroutine linopm (n, nsteps, stepwt)
      implicit none
!
!.... ASSUMES THAT VTURB IS CONSTANT AND THAT THE OPACITY FILE IS GIVEN
!.... ONLY FOR THAT VTURB - THE ODF FOR ALL WAVELENGTHS IS LOADED INTO MEMORY
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision  stepwt
      integer n, nsteps
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.turbpr'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd), 
     &                  plog,  tl, wave, wlbeg, 
     &                  wlend(328), wt(12), x, y
      integer(kind=2)  kap(25, 57, 12), kapwol(25, 57, 12, 328)
      integer    ios, ip, ipj(maxd), istep, it, itemp1, itj(maxd), iwl, 
     &           j, n12
      logical    first
      save       co1, co2, co3, co4, first, ipj, itj, kapwol, wlend
!
!----------------------------- INITIALIZATION --------------------------
!
      data first / .true. /
!
      data itemp1 / 0 /
!
!.... VALUES ARE 1/60, 2/60, 3/60, 6/60
!
      data wt / 1.66666667d-2, 3.33333333d-2, 5.00000000d-2,
     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1,
     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1, 
     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1 /
!
!---------------------------------- EXECUTION --------------------------
!
      if(first) then
         first = .false.
         ios   = 0
         iwl   = 0
!
         do while (ios .eq. 0)
            iwl = iwl + 1
            read(19, rec = iwl, iostat = ios) wlbeg, wlend(iwl), kap
!
            if (ios .eq. 0) then
!
               do ip = 1, 25
!
                  do it = 1, 57
!
                     do istep = 1, 12
                        kapwol(ip, it, istep, iwl) = kap(ip, it, istep)
                     end do
!
                  end do
!
               end do
!
            end if
!
         end do
!
         wlend(328) = 1000000.0d0
         close (unit = 19)
      end if
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
         iwl    = 1
!
         do j = 1, nrhox
            tl = min(max(tlog(j) / tenlog, 3.32d0), 5.3d0)
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. 57)
               it = it + 1
            end do
!
            plog = min(max(log10(p(j)), -2.0d0), 8.0d0)
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. 25)
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
!
!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         end do
!
      end if
!
      wave = 2.997925d17 / freq
!
      do while (wlend(iwl) .le. wave)
         iwl = iwl + 1
      end do
!
      nsteps = 12
      n12    = 13 - n
      stepwt = wt(n)
!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         a  = exp( co1(j) * dble(kapwol(ip-1, it-1, n12, iwl)) +
     &             co2(j) * dble(kapwol(ip  , it-1, n12, iwl)) +
     &             co3(j) * dble(kapwol(ip-1, it  , n12, iwl)) +
     &             co4(j) * dble(kapwol(ip  , it  , n12, iwl)) )
         alines(j) = a
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L I N O P M !!!!!!!!!!!!
!
      subroutine linopv (n, nsteps, stepwt)
      implicit none
!
!.... ODF'S AT 5 MICROTURBULENCES ARE USED
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision  stepwt
      integer n, nsteps
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.constb'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.turbpr'
!
!--------------------------------- CONSTANTS ---------------------------
!
      double precision  ttenlg
      parameter ( ttenlg = 0.001d0 * tenlog )
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision  a, co1(maxd), co2(maxd), co3(maxd), co4(maxd), 
     &                  co5(maxd), co6(maxd), plog, 
     &                  tabv(5), tl, vtur, wave, wlb, wle, wlend(1221),
     &                  wt(12), x, y, z
      integer(kind=2) kap0(25, 57, 12),kap1(25,57,12),kap2(25,57,12),
     &           kap4(25, 57, 12), kap8(25, 57, 12)
      integer(kind=2) kap(25, 57, 12, 5), kapwo(25, 57, 12, 328, 5)
      integer    ios, ip, ipj(maxd), it, itemp1, itj(maxd), iv, 
     &           ivj(maxd), iwl, iwllast, j, n12, istep
      logical    first
      save       co1, co2, co3, co4, co5, co6, first, ipj, itj, ivj, 
     &           iwl, iwllast, kapwo, 
     &           wlb, wle, wlend, wt
!     kap, kap0, kap1, kap2, kap4, kap8,
!
!----------------------------- EQUIVALENCES ----------------------------
!
!.... TO MAP THE DIFFERENT TURBULENT VELOCITIES INTO A SINGLE VARIABLE
!
!     commented 02/12/14 replace by reading in all opacities to kapwo
!      equivalence (kap(1, 1, 1, 1), kap0(1, 1, 1))
!      equivalence (kap(1, 1, 1, 2), kap1(1, 1, 1))
!      equivalence (kap(1, 1, 1, 3), kap2(1, 1, 1))
!      equivalence (kap(1, 1, 1, 4), kap4(1, 1, 1))
!      equivalence (kap(1, 1, 1, 5), kap8(1, 1, 1))
!
!----------------------------- INITIALIZATION --------------------------
!
      data first / .true. /
!
      data itemp1 / 0 /
!
!
      data tabv / 0.0d0, 1.0d5, 2.0d5, 4.0d5, 8.0d5 /
!
!.... VALUES ARE 1/60, 2/60, 3/60, 6/60
!
      data wt / 1.66666667d-2, 3.33333333d-2, 5.00000000d-2,
     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1,
     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1, 
     &          1.00000000d-1, 1.00000000d-1, 1.00000000d-1 /
!
!---------------------------------- EXECUTION --------------------------
!

      if (first) then
         first   = .false.
         ios     = 0
         iwl     = 0
         iwllast = 0
!
!     added 02/12/14 to read all opacities at once.
         do iwl = 1, 328
             read(20,*,iostat=ios) wlb, wlend(iwl)
             read(20,*,iostat=ios) kap0
             read(21,*,iostat=ios) wlb, wlend(iwl)
             read(21,*,iostat=ios) kap1
             read(22,*,iostat=ios) wlb, wlend(iwl)
             read(22,*,iostat=ios) kap2
             read(24,*,iostat=ios) wlb, wlend(iwl)
             read(24,*,iostat=ios) kap4
             read(28,*,iostat=ios) wlb, wlend(iwl)
             read(28,*,iostat=ios) kap8
             do ip = 1, 25
             do it = 1, 57
               do istep = 1, 12
                      kapwo(ip,it,istep,iwl,1) = kap0(ip,it,istep)
                      kapwo(ip,it,istep,iwl,2) = kap1(ip,it,istep)
                      kapwo(ip,it,istep,iwl,3) = kap2(ip,it,istep)
                      kapwo(ip,it,istep,iwl,4) = kap4(ip,it,istep)
                      kapwo(ip,it,istep,iwl,5) = kap8(ip,it,istep)
                   end do
               end do
             end do
           end do
         print *, ' finished reading kapwo '

!         do while (ios .eq. 0)
!            iwl = iwl + 1
!            read(20, rec = iwl, iostat = ios) wlb, wlend(iwl)
!         end do
!
         iwl        = iwl - 1
         wlend(iwl) = 1000000.0d0
         iwl   = 1
      end if
!
      if (itemp .ne. itemp1) then
         itemp1 = itemp
         iwl    = 1
!
         do j = 1, nrhox
            tl = min(max(tlog(j) / tenlog, 3.32d0), 5.3d0)
            it = 2
!
            do while (tabt(it) .le. tl .and. it .lt. 57)
               it = it + 1
            end do
!
            plog = min(max(log10(p(j)), -2.0d0), 8.0d0)
            ip   = 2
!
            do while (tabp(ip) .le. plog .and. ip .lt. 25)
               ip = ip + 1
            end do
!
            vtur = min(max(vturb(j), 0.0d0), 8.0d5)
            iv   = 2
!
            do while (tabv(iv) .le. vtur .and. iv .lt. 5)
               iv = iv + 1
            end do
!
            ipj(j) = ip
            itj(j) = it     
            ivj(j) = iv
            x      = (tl   - tabt(it-1)) / (tabt(it) - tabt(it-1))
            y      = (plog - tabp(ip-1)) / (tabp(ip) - tabp(ip-1))
            z      = (vtur - tabv(iv-1)) / (tabv(iv) - tabv(iv-1))
            co1(j) = (1.0d0 - x) * (1.0d0 - y) * ttenlg
            co2(j) = (1.0d0 - x) * y * ttenlg
            co3(j) = x * (1.0d0 - y) * ttenlg
            co4(j) = x * y * ttenlg
            co5(j) = z
            co6(j) = 1.0d0 - z

!.... THE STEPS HAVE BEEN SCALED BY 1000
!
         end do
!
      end if
!
      wave = 2.997925d17 / freq
!
      do while (wlend(iwl) .le. wave)
         iwl = iwl + 1
      end do
!
!      if (iwl .ne. iwllast) then
!         iwllast = iwl
!         read(20, rec = iwl) wlb, wle, kap0
!         read(21, rec = iwl) wlb, wle, kap1
!         read(22, rec = iwl) wlb, wle, kap2
!         read(24, rec = iwl) wlb, wle, kap4
!         read(28, rec = iwl) wlb, wle, kap8
!      end if
!
      nsteps = 12
      n12    = 13 - n
      stepwt = wt(n)
!
      do j = 1, nrhox
         it = itj(j)
         ip = ipj(j)
         iv = ivj(j)
         a  = exp((co1(j) * dble(kapwo(ip-1, it-1, n12, iwl, iv  )) +
     &             co2(j) * dble(kapwo(ip  , it-1, n12, iwl, iv  )) +
     &             co3(j) * dble(kapwo(ip-1, it  , n12, iwl, iv  )) +
     &             co4(j) * dble(kapwo(ip  , it  , n12, iwl, iv  ))) * 
     &             co5(j) + 
     &            (co1(j) * dble(kapwo(ip-1, it-1, n12, iwl, iv-1)) +
     &             co2(j) * dble(kapwo(ip  , it-1, n12, iwl, iv-1)) +
     &             co3(j) * dble(kapwo(ip-1, it  , n12, iwl, iv-1)) +
     &             co4(j) * dble(kapwo(ip  , it  , n12, iwl, iv-1))) * 
     &             co6(j))
         alines(j) = a
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L I N O P V !!!!!!!!!!!!
!
      subroutine linsop (j, nsteps, stepwt)
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.opsblk'
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision stepwt
      integer j, nsteps
!
!.... DUMMY FOR LINE ABSORPTION DISTRIBUTION FUNCTIONS    S = J
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   L I N S O P !!!!!!!!!!!!
!
      double precision function xkarsas (zeff2, n, l)
      implicit none
!
!.... 1995 JUL - REMOVED freq FROM THE SUBROUTINE CALL
!              - ADDED common.freqbl TO SUPPLY THE VALUE OF freq
!              - CREATED A LOCAL VARIABLE, flog10 = log10(freq/zeff2)
!
!---------------------------- DUMMY VARIABLES --------------------------
!
      double precision  zeff2
      integer n, l
!
!--------------------------------- COMMON -----------------------------
!
!      include 'common.sizebl'
      include 'common.constb'
!      include 'common.freqbl'
!
!---------------------------- LOCAL VARIABLES --------------------------
!
      double precision  ekarsas(29), flog10, freqn(29, 15), freqn15(29),
     &                  x, xl(29, 6, 6), xn(29, 15)
      integer i
!
!------------------------------ INITIALIZATION -------------------------
!
       data ekarsas / 10000.0,   4444.0, 2500.0,   1111.0,   400.0, 
     &                  204.1,    100.0,   44.44,    25.0,    16.0, 
     &                   11.11,     6.25,   4.0,      2.778,   2.041,
     &                    1.562,    1.235,  1.0,      0.6944,  0.4444, 
     &                    0.25,     0.1111, 0.04,     0.02041, 0.01, 
     &                    0.004444, 0.0025, 0.001111, 0.0 /
!
!.... LEVEL N = 1
!
      data (freqn(i, 1), i = 1, 29) /
     &    19.516982,   19.164810,   18.915052,   18.563043,   18.120083,
     &    17.828904,   17.521260,   17.174377,   16.931912,   16.747387,
     &    16.600083,   16.377277,   16.215909,   16.094200,   15.999955,
     &    15.925518,   15.866216,   15.817969,   15.745954,   15.676626,
     &    15.613849,   15.562692,   15.533972,   15.525713,   15.521260,
     &    15.518864,   15.518023,   15.517421,   15.516939/
!
      data (xn(i, 1), i = 1, 29) /
     &   -30.274422,  -29.048572,  -28.181067,  -26.962272,  -25.437868,
     &   -24.444170,  -23.404269,  -22.248421,  -21.454163,  -20.858944,
     &   -20.390346,  -19.694283,  -19.200905,  -18.835387,  -18.556686,
     &   -18.339364,  -18.168213,  -18.030238,  -17.826632,  -17.633456,
     &   -17.461067,  -17.322353,  -17.245241,  -17.223162,  -17.211266,
     &   -17.204840,  -17.202587,  -17.200999,  -17.199715/
!
!.... L = S
!
      data (xl(i, 1, 1), i = 1, 29) /
     &   -30.274422,  -29.048572,  -28.181067,  -26.962272,  -25.437868,
     &   -24.444170,  -23.404269,  -22.248421,  -21.454163,  -20.858944,
     &   -20.390346,  -19.694283,  -19.200905,  -18.835387,  -18.556686,
     &   -18.339364,  -18.168213,  -18.030238,  -17.826632,  -17.633456,
     &   -17.461067,  -17.322353,  -17.245241,  -17.223162,  -17.211266,
     &   -17.204840,  -17.202587,  -17.200999,  -17.199715/
!
!.... LEVEL N = 2
!
      data (freqn(i, 2), i = 1, 29) /
     &    19.516949,   19.164737,   18.914922,   18.562750,   18.119270,
     &    17.827313,   17.518023,   17.167149,   16.919200,   16.727792,
     &    16.572317,   16.329852,   16.145327,   15.998094,   15.876964,
     &    15.775097,   15.688665,   15.613849,   15.492095,   15.358548,
     &    15.215909,   15.074566,   14.979337,   14.948961,   14.931912,
     &    14.922531,   14.919200,   14.916804,   14.914879/
!
      data (xn(i, 2), i = 1, 29) /
     &   -31.779474,  -30.553459,  -29.685827,  -28.466543,  -26.940432,
     &   -25.943993,  -24.898608,  -23.729491,  -22.917021,  -22.298979,
     &   -21.803393,  -21.042629,  -20.473370,  -20.025469,  -19.660029,
     &   -19.355246,  -19.098003,  -18.876442,  -18.517855,  -18.127425,
     &   -17.714170,  -17.308930,  -17.038908,  -16.953361,  -16.905447,
     &   -16.879127,  -16.869826,  -16.863085,  -16.857754/
!
!.... L = S
!
      data (xl(i, 1, 2), i = 1, 29) /
     &   -31.177414,  -29.951530,  -29.083850,  -27.864712,  -26.339031,
     &   -25.343652,  -24.299685,  -23.134693,  -22.327692,  -21.716473,
     &   -21.228927,  -20.487480,  -19.941059,  -19.517455,  -19.178033,
     &   -18.899376,  -18.668043,  -18.471683,  -18.160149,  -17.830286,
     &   -17.492277,  -17.172499,  -16.965517,  -16.901255,  -16.865263,
     &   -16.845632,  -16.838714,  -16.833696,  -16.829681/
!
!.... L = P
!
      data (xl(i, 2, 2), i = 1, 29) /
     &   -35.779538,  -34.184208,  -33.083933,  -31.512708,  -29.543604,
     &   -28.256123,  -26.903279,  -25.387738,  -24.333408,  -23.531477,
     &   -22.889415,  -21.907557,  -21.178842,  -20.610306,  -20.152156,
     &   -19.774043,  -19.458248,  -19.189136,  -18.759267,  -18.299831,
     &   -17.823327,  -17.365980,  -17.066362,  -16.972218,  -16.919695,
     &   -16.890892,  -16.880696,  -16.873357,  -16.867478/
!
!.... LEVEL N = 3
!
      data (freqn(i, 3), i = 1, 29) /
     &    19.516943,   19.164723,   18.914898,   18.562696,   18.119119,
     &    17.827018,   17.517421,   17.165797,   16.916804,   16.724064,
     &    16.566974,   16.320472,   16.130898,   15.977703,   15.849803,
     &    15.740463,   15.646019,   15.562696,   15.423010,   15.261631,
     &    15.074579,   14.863704,   14.696235,   14.635934,   14.600123,
     &    14.579728,   14.572359,   14.567017,   14.562696/
!
      data (xn(i, 3), i = 1, 29) /
     &   -32.659912,  -31.433874,  -30.566210,  -29.346836,  -27.820290,
     &   -26.823453,  -25.777089,  -24.605440,  -23.789519,  -23.167057,
     &   -22.666147,  -21.891933,  -21.306393,  -20.839041,  -20.451712,
     &   -20.122889,  -19.840361,  -19.591597,  -19.176587,  -18.699419,
     &   -18.149566,  -17.533628,  -17.049033,  -16.875774,  -16.773227,
     &   -16.714935,  -16.693926,  -16.678663,  -16.666369/
!
!.... L = S
!
      data (xl(i, 1, 3), i = 1, 29) /
     &   -31.705705,  -30.479739,  -29.612265,  -28.392746,  -26.866974,
     &   -25.871133,  -24.826672,  -23.659806,  -22.850344,  -22.235989,
     &   -21.744734,  -20.993964,  -20.435725,  -19.998364,  -19.643303,
     &   -19.347420,  -19.097776,  -18.881962,  -18.529746,  -18.137370,
     &   -17.701228,  -17.231454,  -16.873769,  -16.748412,  -16.674666,
     &   -16.633129,  -16.617776,  -16.606984,  -16.598091/
!
!.... L = P
!
      data (xl(i, 2, 3), i = 1, 29) /
     &   -36.234105,  -34.655854,  -33.538432,  -31.967064,  -29.997698,
     &   -28.709867,  -27.356451,  -25.839127,  -24.782259,  -23.977343,
     &   -23.331485,  -22.340276,  -21.599900,  -21.017917,  -20.544424,
     &   -20.149344,  -19.815760,  -19.527654,  -19.058410,  -18.538322,
     &   -17.967020,  -17.364676,  -16.918642,  -16.765111,  -16.675798,
     &   -16.625318,  -16.607492,  -16.594210,  -16.583614/
!
!.... L = D
!
      data (xl(i, 3, 3), i = 1, 29) /
     &   -41.364414,  -39.434006,  -38.066663,  -36.143204,  -33.730242,
     &   -32.150245,  -30.487089,  -28.617809,  -27.311427,  -26.313205,
     &   -25.509946,  -24.270587,  -23.339149,  -22.602299,  -21.924436,
     &   -21.493723,  -21.063954,  -20.691590,  -20.080654,  -19.397357,
     &   -18.637161,  -17.823176,  -17.209853,  -16.996234,  -16.871214,
     &   -16.800539,  -16.775144,  -16.756765,  -16.741919/
!
!.... LEVEL N = 4
!
      data (freqn(i, 4), i = 1, 29) /
     &    19.516941,   19.164719,   18.914889,   18.562677,   18.119066,
     &    17.826915,   17.517210,   17.165323,   16.915963,   16.722752,
     &    16.565089,   16.317140,   16.125732,   15.970333,   15.839881,
     &    15.727658,   15.630046,   15.543267,   15.395977,   15.221861,
     &    15.011789,   14.756488,   14.527662,   14.435545,   14.377277,
     &    14.342650,   14.329852,   14.320471,   14.312819/
!
      data (xn(i, 4), i = 1, 29) /
     &   -33.284599,  -32.058554,  -31.190879,  -29.971473,  -28.444826,
     &   -27.447836,  -26.401066,  -25.228582,  -24.411413,  -23.787317,
     &   -23.284581,  -22.505775,  -21.914353,  -21.439606,  -21.044235,
     &   -20.705972,  -20.413135,  -20.153596,  -19.714525,  -19.197426,
     &   -18.576241,  -17.824248,  -17.155428,  -16.887819,  -16.719154,
     &   -16.619216,  -16.582315,  -16.555295,  -16.533276/
!
!.... L = S
!
      data (xl(i, 1, 4), i = 1, 29) /
     &   -32.080641,  -30.854674,  -29.986801,  -28.767697,  -27.241693,
     &   -26.245685,  -25.200974,  -24.033538,  -23.223063,  -22.607845,
     &   -22.115266,  -21.360872,  -20.798453,  -20.355878,  -19.995174,
     &   -19.692644,  -19.435600,  -19.211713,  -18.841933,  -18.420428,
     &   -17.932110,  -17.363567,  -16.873130,  -16.680219,  -16.559751,
     &   -16.488746,  -16.462241,  -16.443053,  -16.427763/
!
!.... L = P
!
      data (xl(i, 2, 4), i = 1, 29) /
     &   -36.585694,  -35.007703,  -33.890016,  -32.318668,  -30.349350,
     &   -29.061334,  -27.707618,  -26.189677,  -25.132040,  -24.325956,
     &   -23.678826,  -22.684226,  -21.939671,  -21.352566,  -20.873369,
     &   -20.471723,  -20.130813,  -19.835172,  -19.348733,  -18.800381,
     &   -18.178384,  -17.480038,  -16.904760,  -16.685329,  -16.550262,
     &   -16.471169,  -16.442151,  -16.420831,  -16.403759/
!
!.... L = D
!
      data (xl(i, 3, 4), i = 1, 29) /
     &   -41.585694,  -39.655304,  -38.288039,  -36.364454,  -33.951410,
     &   -32.371226,  -30.707789,  -28.837992,  -27.530994,  -26.531796,
     &   -25.727043,  -24.484484,  -23.549206,  -22.807462,  -22.198909,
     &   -21.686891,  -21.250382,  -20.870478,  -20.243060,  -19.532238,
     &   -18.722925,  -17.815346,  -17.075994,  -16.798160,  -16.628568,
     &   -16.529588,  -16.493472,  -16.467238,  -16.445815/
!
!.... L = F
!
      data (xl(i, 4, 4), i = 1, 29) /
     &   -47.062815,  -44.780358,  -43.163100,  -40.887314,  -38.030685,
     &   -36.158301,  -34.185235,  -31.963719,  -30.407089,  -29.214529,
     &   -28.252197,  -26.761810,  -25.634821,  -24.737662,  -23.998757,
     &   -23.374580,  -22.839980,  -22.373323,  -21.598611,  -20.713453,
     &   -19.693804,  -18.530997,  -17.563112,  -17.193424,  -16.965517,
     &   -16.832288,  -16.783370,  -16.747717,  -16.718672/
!
!.... LEVEL N = 5
!
      data (freqn(i, 5), i = 1, 29) /
     &    19.516940,   19.164717,   18.914886,   18.562668,   18.119042,
     &    17.826867,   17.517112,   17.165103,   16.915573,   16.722143,
     &    16.564213,   16.315589,   16.123320,   15.966880,   15.835211,
     &    15.721601,   15.622449,   15.533972,   15.382871,   15.202143,
     &    14.979337,   14.696203,   14.420029,   14.298047,   14.215909,
     &    14.164752,   14.145327,   14.130897,   14.118999/
!
      data (xn(i, 5), i = 1, 29) /
     &   -33.769146,  -32.543097,  -31.675417,  -30.455996,  -28.929303,
     &   -27.932243,  -26.885239,  -25.712408,  -24.894628,  -24.269941,
     &   -23.766226,  -22.985245,  -22.390846,  -21.912586,  -21.513577,
     &   -21.170761,  -20.873304,  -20.608270,  -20.156957,  -19.619181,
     &   -18.958075,  -18.121143,  -17.308727,  -16.951892,  -16.712503,
     &   -16.563827,  -16.507488,  -16.465627,  -16.431184/
!
!.... L = S
!
      data (xl(i, 1, 5), i = 1, 29) /
     &   -32.371142,  -31.145245,  -30.277611,  -29.058332,  -27.532386,
     &   -26.536299,  -25.491539,  -24.323724,  -23.512880,  -22.897091,
     &   -22.403960,  -21.648140,  -21.083702,  -20.638728,  -20.275002,
     &   -19.969127,  -19.708598,  -19.480857,  -19.102318,  -18.665521,
     &   -18.148008,  -17.516456,  -16.921283,  -16.663742,  -16.492247,
     &   -16.386117,  -16.345903,  -16.316173,  -16.291778/
!
!.... L = P
!
      data (xl(i, 2, 5), i = 1, 29) /
     &   -36.866137,  -35.287883,  -34.170413,  -32.599199,  -30.629663,
     &   -29.341564,  -27.987755,  -26.469536,  -25.411517,  -24.604882,
     &   -23.957191,  -22.961135,  -22.214481,  -21.625034,  -21.142933,
     &   -20.738297,  -20.393941,  -20.094254,  -19.599261,  -19.036165,
     &   -18.385686,  -17.626125,  -16.948476,  -16.665818,  -16.480643,
     &   -16.367024,  -16.324502,  -16.292865,  -16.266917/
!
!.... L = D
!
      data (xl(i, 3, 5), i = 1, 29) /
     &   -41.816885,  -39.886598,  -38.519116,  -36.595706,  -34.182651,
     &   -32.602365,  -30.938792,  -29.068803,  -27.761491,  -26.761551,
     &   -25.956256,  -24.712472,  -23.775049,  -23.031086,  -22.420027,
     &   -21.905038,  -21.464940,  -21.081321,  -20.445565,  -19.720393,
     &   -18.883701,  -17.916497,  -17.077571,  -16.738117,  -16.519620,
     &   -16.387033,  -16.337715,  -16.301341,  -16.271391/
!
!.... L = F
!
      data (xl(i, 4, 5), i = 1, 29) /
     &   -47.128880,  -44.846322,  -43.229046,  -40.953347,  -38.096716,
     &   -36.224291,  -34.250943,  -32.029199,  -30.472360,  -29.279276,
     &   -28.316408,  -26.824527,  -25.695751,  -24.796176,  -24.054627,
     &   -23.427631,  -22.889877,  -22.419401,  -21.636478,  -20.737351,
     &   -19.690904,  -18.469715,  -17.404053,  -16.973748,  -16.697901,
     &   -16.531879,  -16.469784,  -16.423961,  -16.386588/
!
!.... L = G
!
      data (xl(i, 5, 5), i = 1, 29) /
     &   -52.894711,  -50.260082,  -48.392958,  -45.765034,  -42.464679,
     &   -40.300146,  -38.017153,  -35.443424,  -33.636754,  -32.250427,
     &   -31.129593,  -29.389103,  -28.068001,  -27.012118,  -26.138711,
     &   -25.398332,  -24.761042,  -24.202462,  -23.268415,  -22.188504,
     &   -20.919298,  -19.415147,  -18.073478,  -17.521544,  -17.163795,
     &   -16.946562,  -16.865194,  -16.805098,  -16.755865/
!
!.... LEVEL N = 6
!
      data (freqn(i, 6), i = 1, 29) /
     &    19.516940,   19.164715,   18.914883,   18.562663,   18.119029,
     &    17.826841,   17.517059,   17.164984,   16.915361,   16.721812,
     &    16.563737,   16.314744,   16.122004,   15.964992,   15.832652,
     &    15.718275,   15.618265,   15.528838,   15.375583,   15.191044,
     &    14.960636,   14.659571,   14.348026,   14.199875,   14.094175,
     &    14.025088,   13.998063,   13.977668,   13.960636/
!
      data (xn(i, 6), i = 1, 29) /
     &   -34.165051,  -32.939000,  -32.071317,  -30.851888,  -29.325169,
     &   -28.328071,  -27.280986,  -26.107892,  -25.289843,  -24.664705,
     &   -24.160564,  -23.378190,  -22.782394,  -22.302428,  -21.901012,
     &   -21.555896,  -21.255472,  -20.987585,  -20.529803,  -19.979782,
     &   -19.295022,  -18.402541,  -17.482757,  -17.047424,  -16.737838,
     &   -16.536084,  -16.457331,  -16.397931,  -16.348398/
!
!.... L = S
!
      data (xl(i, 1, 6), i = 1, 29) /
     &   -32.608820,  -31.382756,  -30.515126,  -29.295866,  -27.769793,
     &   -26.773814,  -25.728819,  -24.560932,  -23.750086,  -23.133811,
     &   -22.640288,  -21.883631,  -21.318035,  -20.871913,  -20.506426,
     &   -20.198858,  -19.936428,  -19.706400,  -19.322760,  -18.877373,
     &   -18.342274,  -17.669792,  -16.995256,  -16.680122,  -16.457336,
     &   -16.312694,  -16.256489,  -16.214113,  -16.178612/
!
!.... L = P
!
      data (xl(i, 2, 6), i = 1, 29) /
     &   -37.098169,  -35.519950,  -34.402525,  -32.831070,  -30.861699,
     &   -29.573885,  -28.219694,  -26.701459,  -25.643044,  -24.836230,
     &   -24.188105,  -23.191275,  -22.443490,  -21.852666,  -21.369042,
     &   -20.962634,  -20.616374,  -20.314553,  -19.814673,  -19.242970,
     &   -18.575541,  -17.775947,  -17.020568,  -16.681448,  -16.445735,
     &   -16.294606,  -16.235710,  -16.191866,  -16.154983/
!
!.... L = D
!
      data (xl(i, 3, 6), i = 1, 29) /
     &   -42.024362,  -40.094064,  -38.726686,  -36.803137,  -34.390124,
     &   -32.809866,  -31.146180,  -29.276029,  -27.968300,  -26.968324,
     &   -26.162701,  -24.918051,  -23.979662,  -23.234506,  -22.621799,
     &   -22.105162,  -21.663212,  -21.277514,  -20.637026,  -19.903484,
     &   -19.050185,  -18.044511,  -17.129904,  -16.735338,  -16.467566,
     &   -16.298269,  -16.232977,  -16.184230,  -16.143922/
!
!.... L = F
!
      data (xl(i, 4, 6), i = 1, 29) /
     &   -47.267412,  -44.984913,  -43.367636,  -41.091842,  -38.235239,
     &   -36.362731,  -34.389528,  -32.167518,  -30.610443,  -29.417223,
     &   -28.453971,  -26.961283,  -25.831491,  -24.930907,  -24.187725,
     &   -23.559075,  -23.019383,  -22.547066,  -21.759545,  -20.852145,
     &   -19.789541,  -18.530522,  -17.390884,  -16.906727,  -16.582667,
     &   -16.380139,  -16.302886,  -16.245236,  -16.197380/
!
!.... L = G
!
      data (xl(i, 5, 6), i = 1, 29) /
     &   -52.845039,  -50.210247,  -48.343069,  -45.715131,  -42.414728,
     &   -40.250164,  -37.967149,  -35.393156,  -33.586496,  -32.199833,
     &   -31.078643,  -29.337458,  -27.969702,  -26.958401,  -26.083595,
     &   -25.341555,  -24.702345,  -24.141808,  -23.203287,  -22.115356,
     &   -20.830007,  -19.288694,  -17.874057,  -17.268729,  -16.863465,
     &   -16.610369,  -16.513883,  -16.442010,  -16.382570/
!
!.... L = H
!
      data (xl(i, 6, 6), i = 1, 29) /
     &   -58.850334,  -55.863542,  -53.746437,  -50.766409,  -47.022317,
     &   -44.565391,  -41.972509,  -39.046704,  -36.990356,  -35.410261,
     &   -34.131188,  -32.140740,  -30.626018,  -29.411767,  -28.404701,
     &   -27.548439,  -26.808936,  -26.159088,  -25.067378,  -23.795088,
     &   -22.279431,  -20.436907,  -18.711058,  -17.957760,  -17.446882,
     &   -17.124901,  -17.001376,  -16.909196,  -16.832806/
!
!.... LEVEL N = 7
!
      data (freqn(i, 7), i = 1, 29) /
     &    19.516939,   19.164715,   18.914882,   18.562661,   18.119021,
     &    17.826825,   17.517027,   17.164912,   16.915233,   16.721612,
     &    16.563450,   16.314234,   16.121209,   15.963850,   15.831103,
     &    15.716257,   15.615723,   15.525712,   15.371128,   15.184212,
     &    14.948958,   14.635891,   14.298034,   14.127792,   13.999929,
     &    13.912303,   13.876929,   13.849764,   13.826742/
!
      data (xn(i, 7), i = 1, 29) /
     &   -34.499784,  -33.273731,  -32.406047,  -31.186614,  -29.659879,
     &   -28.662758,  -27.615624,  -26.442410,  -25.624138,  -24.998790,
     &   -24.494343,  -23.711394,  -23.114332,  -22.633333,  -22.230699,
     &   -21.884181,  -21.582185,  -21.312152,  -20.849982,  -20.292819,
     &   -19.593097,  -18.663739,  -17.663648,  -17.161477,  -16.785637,
     &   -16.528798,  -16.425342,  -16.345983,  -16.278790/
!
!.... LEVEL N = 8
!
      data (freqn(i, 8), i = 1, 29) /
     &    19.516939,   19.164714,   18.914881,   18.562659,   18.119016,
     &    17.826815,   17.517006,   17.164865,   16.915150,   16.721482,
     &    16.563263,   16.313903,   16.120692,   15.963107,   15.830094,
     &    15.714942,   15.614066,   15.523672,   15.368212,   15.179720,
     &    14.941207,   14.619801,   14.262209,   14.073663,   13.925602,
     &    13.819464,   13.775217,   13.740590,   13.710759/
!
      data (xn(i, 8), i = 1, 29) /
     &   -34.789743,  -33.563690,  -32.696004,  -31.476568,  -29.949823,
     &   -28.952576,  -27.905521,  -26.732230,  -25.913849,  -25.288312,
     &   -24.783697,  -24.000359,  -23.402741,  -22.921064,  -22.517235,
     &   -22.169801,  -21.866776,  -21.595595,  -21.130798,  -20.568503,
     &   -19.858590,  -18.903358,  -17.843146,  -17.285660,  -16.849210,
     &   -16.537235,  -16.407454,  -16.306014,  -16.218699/
!
!.... LEVEL N = 9
!
      data (freqn(i, 9), i = 1, 29) /
     &    19.516939,   19.164714,   18.914881,   18.562657,   18.119012,
     &    17.826808,   17.516992,   17.164833,   16.915093,   16.721394,
     &    16.563135,   16.313676,   16.120337,   15.962597,   15.829401,
     &    15.714039,   15.612925,   15.522267,   15.366202,   15.176613,
     &    14.935812,   14.608414,   14.235819,   14.032225,   13.866132,
     &    13.741981,   13.688539,   13.645876,   13.608454/
!
      data (xn(i, 9), i = 1, 29) /
     &   -35.045505,  -33.819451,  -32.951765,  -31.732326,  -30.205575,
     &   -29.208318,  -28.161241,  -26.987832,  -26.169441,  -25.543807,
     &   -25.039029,  -24.255440,  -23.657439,  -23.175297,  -22.770919,
     &   -22.422852,  -22.118723,  -21.846749,  -21.380133,  -20.814545,
     &   -20.097359,  -19.123314,  -18.017622,  -17.414518,  -16.923750,
     &   -16.558183,  -16.401026,  -16.275647,  -16.165911/
!
!.... LEVEL N = 10
!
      data (freqn(i, 10), i = 1, 29) /
     &    19.516939,   19.164714,   18.914880,   18.562657,   18.119009,
     &    17.826803,   17.516982,   17.164810,   16.915052,   16.721330,
     &    16.563043,   16.313513,   16.120083,   15.962231,   15.828904,
     &    15.713391,   15.612108,   15.521260,   15.364758,   15.174377,
     &    14.931912,   14.600083,   14.215909,   13.999955,   13.817969,
     &    13.676626,   13.613849,   13.562692,   13.516939/
!
      data (xn(i, 10), i = 1, 29) /
     &   -35.274293,  -34.048238,  -33.180551,  -31.961111,  -30.434355,
     &   -29.437090,  -28.389998,  -27.216550,  -26.398051,  -25.772354,
     &   -25.267495,  -24.483312,  -23.885464,  -23.402587,  -22.997820,
     &   -22.649302,  -22.344664,  -22.072514,  -21.604193,  -21.035827,
     &   -20.313639,  -19.326284,  -18.184568,  -17.544349,  -17.005732,
     &   -16.588554,  -16.403642,  -16.253350,  -16.118795/
!
!.... LEVEL N = 11
!
      data (freqn(i, 11), i = 1, 29) /
     &    19.516939,   19.164713,   18.914880,   18.562656,   18.119008,
     &    17.826799,   17.516974,   17.164793,   16.915022,   16.721283,
     &    16.562976,   16.313392,   16.119895,   15.961961,   15.828537,
     &    15.712911,   15.611502,   15.520513,   15.363687,   15.172715,
     &    14.929003,   14.593814,   14.200566,   13.974434,   13.778545,
     &    13.621032,   13.548931,   13.488931,   13.434153/
!
      data (xn(i, 11), i = 1, 29) /
     &   -35.481256,  -34.255201,  -33.387514,  -32.168073,  -30.641313,
     &   -29.644043,  -28.596939,  -27.423463,  -26.604924,  -25.979176,
     &   -25.474255,  -24.689915,  -24.091864,  -23.608739,  -23.203681,
     &   -22.854826,  -22.549810,  -22.276842,  -21.807547,  -21.237407,
     &   -20.511071,  -19.513620,  -18.342986,  -17.672186,  -17.092253,
     &   -16.625647,  -16.412652,  -16.237373,  -16.076228/    
!---- 23jul03 FIVE bad numbers in X level 11 -> see Castelli
!     &   -20.511071,  -19.513620,  -18.342150,  -17.667949,  -17.093121,
!     &   -16.627232,  -16.414294,  -16.237373,  -16.076228/
!
!.... LEVEL N = 12
!
      data (freqn(i, 12), i = 1, 29) /
     &    19.516939,   19.164713,   18.914880,   18.562655,   18.119006,
     &    17.826796,   17.516969,   17.164780,   16.914999,   16.721247,
     &    16.562924,   16.313301,   16.119752,   15.961755,   15.828257,
     &    15.712546,   15.611041,   15.519944,   15.362870,   15.171447,
     &    14.926778,   14.588984,   14.188523,   13.953966,   13.745966,
     &    13.573403,   13.492115,   13.423028,   13.358576/
!
      data (xn(i, 12), i = 1, 29) /
     &   -35.670198,  -34.444144,  -33.576456,  -32.357014,  -30.830251,
     &   -29.832977,  -28.785864,  -27.612367,  -26.793798,  -26.168012,
     &   -25.663043,  -24.878583,  -24.280378,  -23.797065,  -23.391784,
     &   -23.042673,  -22.737368,  -22.464078,  -21.994040,  -21.422148,
     &   -20.692935,  -19.687256,  -18.494545,  -17.795069,  -17.182159,
     &   -16.669643,  -16.429381,  -16.227310,  -16.037494/

!----- see Castelli version bad number corrections!
!
!     &   -20.692935,  -19.687256,  -18.494545,  -17.795069,  -17.183891,
!     &   -16.673156,  -16.431990,  -16.227310,  -16.037494/
!
!.... LEVEL N = 13
!
      data (freqn(i, 13), i = 1, 29) /
     &    19.516939,   19.164713,   18.914880,   18.562655,   18.119005,
     &    17.826794,   17.516964,   17.164770,   16.914981,   16.721219,
     &    16.562884,   16.313230,   16.119641,   15.961595,   15.828039,
     &    15.712262,   15.610681,   15.519501,   15.362233,   15.170457,
     &    14.925038,   14.585188,   14.178914,   13.937343,   13.718804,
     &    13.532347,   13.442104,   13.363780,   13.289052/
!
      data (xn(i, 13), i = 1, 29) /
     &   -35.844009,  -34.617954,  -33.750266,  -32.530823,  -31.004058,
     &   -30.006781,  -28.959661,  -27.786148,  -26.967555,  -26.341739,
     &   -25.836687,  -25.051753,  -24.453445,  -23.969994,  -23.564544,
     &   -23.215236,  -22.909707,  -22.636559,  -22.165546,  -21.592592,
     &   -20.861125,  -19.849269,  -18.640363,  -17.921966,  -17.273191,
     &   -16.719020,  -16.451969,  -16.222218,  -16.001878/

!---- see Castelli version; bad number corrections!
!
!     &   -20.861125,  -19.849269,  -18.639111,  -17.918166,  -17.276217,
!     &   -16.722786,  -16.454970,  -16.222218,  -16.001878/
!
!.... LEVEL N = 14
!
      data (freqn(i, 14), i = 1, 29) /
     &    19.516939,   19.164713,   18.914879,   18.562655,   18.119004,
     &    17.826792,   17.516961,   17.164762,   16.914967,   16.721197,
     &    16.562852,   16.313173,   16.119552,   15.961468,   15.827866,
     &    15.712036,   15.610396,   15.519149,   15.361728,   15.169670,
     &    14.923652,   14.582152,   14.171135,   13.923684,   13.695974,
     &    13.496762,   13.397869,   13.310243,   13.224682/
!
      data (xn(i, 14), i = 1, 29) /
     &   -36.004932,  -34.778877,  -33.911189,  -32.691746,  -31.164979,
     &   -30.167699,  -29.120574,  -27.947047,  -27.128436,  -26.502596,
     &   -25.997515,  -25.212506,  -24.614103,  -24.130536,  -23.724949,
     &   -23.375482,  -23.069774,  -22.796032,  -22.324557,  -21.750758,
     5   -21.017491,  -20.000677,  -18.777116,  -18.041065,  -17.364348,
     6   -16.772813,  -16.479089,  -16.221551,  -15.968930/
!
!--- see above corrections
!     &   -21.017491,  -20.000677,  -18.776282,  -18.037692,  -17.368650,
!     &   -16.776515,  -16.482568,  -16.221551,  -15.968930/
!
!.... LEVEL N = 15
!
      data (freqn(i, 15), i = 1, 29) /
     &    19.516939,   19.164713,   18.914879,   18.562654,   18.119003,
     &    17.826791,   17.516958,   17.164756,   16.914956,   16.721179,
     &    16.562826,   16.313127,   16.119481,   15.961365,   15.827726,
     &    15.711854,   15.610166,   15.518864,   15.361319,   15.169034,
     &    14.922532,   14.579688,   14.164756,   13.912343,   13.676639,
     &    13.465764,   13.358576,   13.261657,   13.164756/
!
      data (xn(i, 15), i = 1, 29) /
     &   -36.154748,  -34.928693,  -34.061005,  -32.841561,  -31.314793,
     &   -30.317511,  -29.270382,  -28.096844,  -27.278218,  -26.652358,
     &   -26.147254,  -25.362186,  -24.763705,  -24.280044,  -23.874346,
     &   -23.524751,  -23.218899,  -22.944996,  -22.473148,  -21.898667,
     &   -21.163944,  -20.143099,  -18.907170,  -18.155759,  -17.454858,
     &   -16.827663,  -16.509932,  -16.224591,  -15.938340/
!----------------!
!     &   -21.163944,  -20.143099,  -18.906962,  -18.152646,  -17.460462,
!     &   -16.832855,  -16.513888,  -16.224591,  -15.938340/
!
!----------------------------- EXECUTION -------------------------------
!
      flog10 = log10(freq / zeff2)
      xkarsas = 0.0d0
!
      if (l .ge. n .or. n .gt. 6) then
!
         if (n .gt. 15) then
            freqn15(29) = log10(109677.576 * 2.99792458d10 / n**2)
!
            if(flog10 .ge. freqn15(29)) then
               i = 2
               freqn15(i) = log10((ekarsas(i) + 1.0 / n**2) * 
     &                      109677.576 * 2.99792458d10)
!
               do while(flog10 .le. freqn15(i) .and. i .lt. 29)
                  i = i + 1
                  if (i .lt. 29) freqn15(i) = log10((ekarsas(i) + 1.0 / 
     &               n**2) * 109677.576 * 2.99792458d10)
               end do
!
               x = (flog10 - freqn15(i)) / (freqn15(i-1)-freqn15(i))*
     &             (xn(i-1,15) - xn(i,15)) + xn(i,15)
               xkarsas = exp(x * tenlog) / zeff2
            end if
!
         else if (flog10 .ge. freqn(29, n)) then
!
            i = 2
!
            do while(flog10 .le. freqn(i, n) .and. i .lt. 29)
               i = i + 1
            end do
!
            x = (flog10 - freqn(i,n)) / (freqn(i-1,n) - freqn(i,n)) *
     &          (xn(i-1,n) - xn(i,n)) + xn(i,n)
            xkarsas = exp(x * tenlog) / zeff2
         end if
!
      else if (flog10 .ge. freqn(29,n)) then
         i = 2
!
         do while (flog10 .le. freqn(i,n) .and. i .lt. 29)
            i = i + 1
         end do
!
         x = (flog10 - freqn(i,n)) / (freqn(i-1,n) - freqn(i,n)) *
     &       (xl(i-1,l+1,n) - xl(i,l+1,n)) + xl(i,l+1,n)
         xkarsas = exp(x * tenlog) / zeff2
      end if
!
      end function 

!
!!!!!!!!!!!!!!! E N D   O F   F U N C T I O N   X K A R S A S !!!!!!!!!!
!
#ifdef NHOP
      subroutine nhop
      implicit none
! -- 2Sig- X
! -- 2Pi - X
! -- APi - X
! -- cross-section taken from 
! -- Gang Shen, A. Kuri, J. M. Fontenla, P. C. Stancil, and Jianguo Wang 
! -- (unpublished) 

!---- COMMON
      include 'common.constb'
      include 'common.freqbl'
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.xnfpbl'
     
      integer j
!--- all 3 X-sections added, only until 520nm 
      if (freq .ge. 5.76525d14) then 
        do j = 1, nrhox
          anhop(j) = xnfpnh(j) * nhallx(j)*(1.0d-16) * stim(j) / rho(j) 
        end do 
      end if 

!--- three seperate cross-section calculations : 

!---- 2Sig -x and APi - X until 520nm
!      if (freq .ge. 5.76525d14) then
!
!        do j = 1, nrhox
!         anhop(j) = xnfpnh(j)*(nhsigx(j)+nhapix(j) )*(1.0d-16)*stim(j)
!     &                / rho(j)
!        end do 
!---- --- 2Pi -X until 346.7 nm 
!
!        if (freq .ge. 8.64705e+14 ) then 
!          do j = 1, nrhox
!           anhop(j)=anhop(j)+xnfpnh(j)*nhpix(j)*(1.0d-16)*stim(j)/rho(j)
!         end do 
!        end if 
!      end if 
 

      end subroutine 


      double precision function nhsigx(j) 
      implicit none

      include 'common.tempbl'
      include 'common.constb'
!
!----------------------------- DUMMY VARIABLE --------------------------
!
      integer, intent(in) ::  j
!---- Cross-section 2 ^3 SIGMA -- X^3 SIGMA - from Stancil-Kurucz-communication

      real(kind=8) wavearr1(4701)
      real(kind=8) crosspix(4701,15)
      real(kind=8) tarr(15), freq1 
      real(kind=8) waveiA(1), crosst(15) 

      integer itemp1, numbT, i, numbL
      integer nT,  l 
!-----external
      external map1

      integer map1, maxtemp

!--- initialisation 

      save wavearr1, crosspix, tarr , numbT, numbL 
      save itemp1 
      save freq1, crosst
      data itemp1 / 0 /


      if (itemp1 .ne. itemp) then
        itemp1 = itemp
        freq1 = 0.0d0 
        open (unit = 78, file = '../data/NHXsec/2Sig_X.dat',
     &   status = 'old', form = 'formatted')
        read(78,*) numbL
        read(78,*) numbT
        read(78,*) (tarr(i), i=1, numbT)
        do l= 1, numbL
             read(78,*) wavearr1(l), (crosspix(l,i), i = 1, numbT)
        end do
        close(unit =78)
        
      end if
 
      if (t(j) .le. 10000.0) then 
         if (freq1 .ne. freq ) then 
           freq1 = freq
           waveiA(1) = c_nm*10.0/freq
           do i = 1, numbT 
            call linter(wavearr1,crosspix(1,i),numbL,waveiA,crosst(i),1)
           end do  
         end if 

         maxtemp = map1(tarr, crosst, numbT, t(j), nhsigx, 1)
         if (nhsigx .lt. 0.0d0) nhsigx =0.0d0
      else 
         nhsigx= 0.0d0
      end if 
      return 

7000  continue
      print*, 'NH - crosssection files not found'
      stop 


      end function 


      double precision function nhpix(j)
       implicit none

       include 'common.tempbl' 
       include 'common.constb'
!
!----------------------------- DUMMY VARIABLE --------------------------
!
      integer, intent(in) ::  j
!---- Cross-section 2 ^3 PI -- X^3 SIGMA - from Stancil-Kurucz-communication

      real(kind=8) wavearr1(4701)
      real(kind=8) crosspix(4701,15)
      real(kind=8) tarr(15), freq1 
      real(kind=8) waveiA(1), crosst(15)

      integer itemp1, numbT, i, numbL
      integer nT, l
!-----external
      external map1

      integer map1, maxtemp

!--- initialisation
      save wavearr1, crosspix, tarr , numbT, numbL
      save itemp1
      save freq1, crosst
      data itemp1 / 0 /


      if (itemp1 .ne. itemp) then 
        itemp1 = itemp
        freq1 = 0.0d0
        open (unit = 78, file = '../data/NHXsec/2Pi_X.dat',
     &   status = 'old', form = 'formatted', err=7000)
        read(78,*) numbL
        read(78,*) numbT
        read(78,*) (tarr(i), i=1, numbT)
        do l= 1, numbL 
          read(78,*) wavearr1(l), (crosspix(l,i), i = 1, numbT)        
        end do 
        close(unit =78)

      end if 
     
      if(t(j) .le. 10000) then 
        if (freq1 .ne. freq) then 
          freq1 = freq
          waveiA(1) = c_nm*10.0/freq

          do i = 1, numbT
            call linter(wavearr1,crosspix(1,i),numbL,waveiA,crosst(i),1)
          end do
        end if

         maxtemp = map1(tarr, crosst, numbT, t(j), nhpix, 1)
         if (nhpix .lt. 0.0d0 ) nhpix =0.0d0

      else  
         nhpix = 0.0d0 
      end if 
      return
7000  continue
      print*, 'NH - crosssection files not found'
      stop
 

      end function 

      double precision function nhapix(j)
       implicit none

       include 'common.tempbl'
       include 'common.constb'
!
!----------------------------- DUMMY VARIABLE --------------------------
!
      integer, intent(in) ::  j
!---- Cross-section A ^3 Pi -- X^3 SIGMA - from Stancil-Kurucz-communication
      real(kind=8) wavearr1(4701)
      real(kind=8) crosspix(4701,15)
      real(kind=8) tarr(15), freq1 
      real(kind=8) waveiA(1), crosst(15)

      integer itemp1, numbT, i, numbL
      integer nT,  l
!-----external
      external map1

      integer map1, maxtemp

!--- initialisation

      save wavearr1, crosspix, tarr  , numbT, numbL
      save itemp1
      save freq1, crosst
      data itemp1 / 0 /


      if (itemp1 .ne. itemp) then
        itemp1 = itemp
        freq1 = 0.0d0 
        open (unit = 78, file = '../data/NHXsec/A_APi_X.dat',
     &   status = 'old', form = 'formatted', err = 7000)
        read(78,*) numbL
        read(78,*) numbT
        read(78,*) (tarr(i), i=1, numbT)
        do l= 1, numbL
           read(78,*) wavearr1(l), (crosspix(l,i), i = 1, numbT)
        end do
        close(unit =78)
      
      end if
      
      if (t(j) .le. 10000) then 
        if (freq1 .ne. freq) then 
          freq1 = freq
          waveiA(1) = c_nm*10.0/freq
          do i = 1, numbT
           call linter(wavearr1,crosspix(1,i),numbL,waveiA,crosst(i),1)
          end do
        end if 

         maxtemp = map1(tarr, crosst, numbT, t(j), nhapix, 1)
         if (nhapix .lt. 0.0d0) nhapix = 0.0d0 
 
      else 
        nhapix = 0.0d0 
      endif 
      return

7000  continue
      print*, 'NH - crosssection files not found'
      stop


      end function 

      double precision function nhallx(j)
       implicit none

       include 'common.tempbl'
       include 'common.constb'
!
!----------------------------- DUMMY VARIABLE --------------------------
!
      integer, intent(in) ::  j
!---- Cross-section A ^3 Pi -- X^3 SIGMA - from Stancil-Kurucz-communication
      real(kind=8) wavearr1(4701)
      real(kind=8) crosspix(4701,15), crossnh(4701,15)
      real(kind=8) tarr(15), freq1
      real(kind=8) waveiA(1), crosst(15)

      integer itemp1, numbT, i, numbL
      integer nT, m, l
!-----external
      external map1

      integer map1, maxtemp

!--- initialisation

      save wavearr1, crossnh, tarr  , numbT, numbL
      save itemp1
      save freq1, crosst
      data itemp1 / 0 /


      if (itemp1 .ne. 1) then
        itemp1 = 1 
        freq1 = 0.0d0

!--- first read smallest lambda range : 2Pi -- X^2
        open (unit = 78, file = '../data/NHXsec/2Pi_X.dat',
     &   status = 'old', form = 'formatted', err=7000)
        read(78,*) numbL
        read(78,*) numbT
        read(78,*) (tarr(i), i=1, numbT)
        do l= 1, numbL
          read(78,*) wavearr1(l), (crossnh(l,i), i = 1, numbT)
        end do
        close(unit =78)
! --- 2nd read 2 Sig ^3 --- X^2  

        open (unit = 78, file = '../data/NHXsec/2Sig_X.dat',
     &   status = 'old', form = 'formatted')
        read(78,*) numbL
        read(78,*) numbT
        read(78,*) (tarr(i), i=1, numbT)
        do l= 1, numbL
            read(78,*) wavearr1(l), (crosspix(l,i), i = 1, numbT)
            do m = 1, numbT
              crossnh(l, m) = crossnh(l,m) + crosspix(l,m)
            end do  

        end do
        close(unit =78)


!---- read A ^3 Pi -- X^2 
        open (unit = 78, file = '../data/NHXsec/A_APi_X.dat',
     &   status = 'old', form = 'formatted', err = 7000)
        read(78,*) numbL
        read(78,*) numbT
        read(78,*) (tarr(i), i=1, numbT)
        do l= 1, numbL
           read(78,*) wavearr1(l), (crosspix(l,i), i = 1, numbT)
           do m = 1, numbT
            crossnh(l,m) = crossnh(l,m) + crosspix(l,m)
           end do 
        end do
        close(unit =78)

      end if

      if (t(j) .le. 10000) then
        if (freq1 .ne. freq) then
          freq1 = freq
          waveiA(1) = c_nm*10.0/freq
          do i = 1, numbT
           call linter(wavearr1,crossnh(1,i),numbL,waveiA,crosst(i),1)
          end do
        end if

         maxtemp = map1(tarr, crosst, numbT, t(j), nhallx, 1)
         if (nhallx .lt. 0.0d0) nhallx = 0.0d0

      else
        nhallx = 0.0d0
      endif
      return

7000  continue
      print*, 'NH - crosssection files not found'
      stop


      end function
#endif



!
!!!!!!!!!!!!!!! E N D   O F   S U B    A N D   F U N C T I O N   NHOP  !!!!!!!!!!
!
      subroutine xlinop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
      include 'common.opsblk'
!
!.... DUMMY LINE OPACITY ROUTINE
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   X L I N O P !!!!!!!!!!!!
!
      subroutine xlisop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.opsblk'
!
!.... DUMMY LINE SCATTERING ROUTINE
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   X L I S O P !!!!!!!!!!!!
!
      subroutine xconop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
      include 'common.opsblk'
      include 'common.stateb'
      include 'common.tempbl'
      include 'common.turbpr'
!
!------------------------------- LOCAL VARIABLES -----------------------
!
      double precision  rosstab
      integer j
!
!------------------------------- EXTERNAL FUNCTION ---------------------
!
      external rosstab
!
!------------------------------- EXECUTION -----------------------------
!
      do j = 1, nrhox
         axcont(j) = rosstab(t(j), p(j), vturb(j))
         sxcont(j) = 5.667d-5 / 12.5664d0 * t(j) ** 4 * 4.0d0
      end do
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   X C O N O P !!!!!!!!!!!!
!
      subroutine xsop
      implicit none
!
!-------------------------------- COMMONS -----------------------------
!
!      include 'common.sizebl'
      include 'common.opsblk'
!
!.... DUMMY SCATTERING ROUTINE
!
      end subroutine 
!
!!!!!!!!!!! E N D   O F   S U B R O U T I N E   X S O P !!!!!!!!!!!!!!!!

      end module kappa_cal 
