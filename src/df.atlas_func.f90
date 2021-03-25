      module func_atlas
      implicit none 

      contains

      double precision function expi (n, x)
      implicit none
!
!.... EXPONENTIAL INTEGRAL FOR POSITIVE ARGUMENTS AFTER CODY AND
!... THACHER, MATH. OF COMP.,22,641(1968)
!
!--------------------------- DUMMY ARGUMENTS ---------------------------
!
      double precision  x
      integer  n
!
!--------------------------- LOCAL VARIABLES --------------------------

      real(kind=8)  a0, a1, a2, a3, a4, a5 
      real(kind=8)  b0, b1, b2, b3, b4
      real(kind=8)  c0, c1, c2, c3, c4, c5, c6
      real(kind=8)  d1, d2, d3, d4, d5, d6 
      real(kind=8)  e0, e1, e2, e3, e4, e5, e6, ex, ex1
      real(kind=8)  f1, f2, f3, f4, f5, f6 
      real(kind=8)  x1

      integer  i
      save     ex1, x1
!
!--------------------------- INITIALIZATION ---------------------------
!
       x1 = -1.0d20 
!
!      data a0, a1, a2, a3, a4, a5 /

       a0 =  -44178.5471728217d0
       a1 =   57721.7247139444d0  
       a2 =   9938.31388962037d0
       a3 =   1842.11088668000d0   
       a4 =   101.093806161906d0
       a5 =   5.03416184097568d0 
!
!      data b0, b1, b2, b3, b4 /

       b0 =   76537.3323337614d0
       b1 =   32597.1881290275d0 
       b2 =   6106.10794245759d0
       b3 =   635.419418378382d0 
       b4 =   37.2298352833327d0 
!
!      data c0, c1, c2, c3, c4, c5, c6 /

       c0 =  4.65627107975096d-7
       c1 =  0.999979577051595d0 
       c2 =  9.04161556946329d0
       c3 =  24.3784088791317d0    
       c4 =  23.0192559391333d0
       c5 =   6.90522522784444d0
       c6 =   0.430967839469389d0 
!
!      data d1, d2, d3, d4, d5, d6 /

       d1 =  10.0411643829054d0
       d2 =  32.4264210695138d0
       d3 =  41.2807841891424d0
       d4 =  20.4494785013794d0
       d5 =  3.31909213593302d0
       d6 =  0.103400130404874d0
!
!      data e0, e1, e2, e3, e4, e5, e6 /

       e0 =  -0.999999999998447d0
       e1 =  -26.6271060431811d0 
       e2 =  -241.055827097015d0
       e3 =  -895.927957772937d0 
       e4 =  -1298.85688746484d0
       e5 =  -545.374158883133d0
       e6 =  -5.66575206533869d0 
!
!      data f1, f2, f3, f4, f5, f6 /

       f1 =   28.6271060422192d0
       f2 =   292.310039388533d0 
       f3 =   1332.78537748257d0
       f4 =   2777.61949509163d0 
       f5 =   2404.01713225909d0
       f6 =   631.657483280800d0 
!*
!*------------------------------- EXECUTION -----------------------------
!*
      if (x .ne. x1) then
         x1 = x
         ex = dexp(-x1)
!*
         if (x1 .gt. 4.0d0)  then
            ex1 = (ex + ex * (e0 + (e1 + (e2 + (e3 + (e4 + (e5 + e6 /
     &              x1) / x1) / x1) / x1) / x1) / x1) /
     &       (x1 + f1 +(f2 + (f3 + (f4 + (f5 + f6 / x1) / x1) / x1)
     &               / x1) / x1)) / x1
!*
         else if (x1 .gt. 1.0d0) then
            ex1 = ex * (c6 + (c5 + (c4 + (c3 + (c2 + (c1 + c0 * x1) *
     &            x1) * x1) * x1) * x1) * x1) / (d6 + (d5 + (d4 +
     &            (d3 + (d2 + (d1 + x1) * x1) * x1) * x1) * x1) * x1)
!*
         else if (x1 .gt. 0.0d0) then
            ex1 = (a0 + (a1 + (a2 + (a3 + (a4 + a5 * x1) * x1) * x1) *
     &               x1) * x1) / (b0 + (b1 + (b2 + (b3 + (b4 + x1) *
     &               x1) * x1) * x1) * x1) -dlog(x1)
!*
         else
            ex1 = 0.0d0
         end if
!*
      end if
!*
      expi = ex1
!*
      if (n .gt. 1) then
!
         do i = 1, n - 1
            expi = (ex - x1 * expi) / dble(i)
         end do
!*
      end if
!*
      end function 
!*
!*********** E N D   O F   F U N C T I O N    E X P I *******************




      integer function map1 (xold, fold, nold, xnew, fnew, nnew)
      implicit none
!*
!*--------------------------- DUMMY ARGUMENTS ---------------------------
!*
      double precision  fnew(*), fold(*), xnew(*), xold(*)
      integer  nnew, nold
!*
!*--------------------------- LOCAL VARIABLES --------------------------
!*
      double precision  ab, af, bb, bf, cb, cf, d, db, df, wt
      integer  l, ll, k, lm1, lm2, lp1
!*
!*------------------------------- EXECUTION -----------------------------
!*
      l = 2
      ll = 0

      do k = 1, nnew

         do while (l .le. nold .and. xnew(k) .ge. xold(l))
            l = l + 1
         end do

         if (l .gt. nold) l = nold

         if (l .gt. 2 .and. l .lt. nold) then

!*.... PARABOLIC CASE
!*
            if (l .ne. ll) then

               if (l .gt. 3 .and. l .eq. ll+1) then
                  ab = af
                  bb = bf
                  cb = cf

               else
!*
!*.... MUST COMPUTE THE BACKWARD COEFFICIENTS
!*
                  lm1 = l - 1
                  lm2 = l - 2
                  d = (fold(lm1) - fold(lm2)) /
     &                (xold(lm1) - xold(lm2))
                  cb = ((fold(l) - fold(lm1)) /
     &                  (xold(l) - xold(lm1)) - d) /
     &                 (xold(l) - xold(lm2))
                  bb = d + cb * (xold(lm1) - xold(lm2))
                  ab = fold(lm1)
               end if
!*
               lp1 = l + 1
               lm1 = l - 1
               d = (fold(l) - fold(lm1)) / (xold(l) - xold(lm1))
               cf = ((fold(lp1) - fold(l)) / 
     &               (xold(lp1) - xold(l)) - d) /
     &              (xold(lp1) - xold(lm1))
               bf = d + cf * (xold(l) - xold(lm1))
               af = fold(l)
               wt = 0.0d0
               if (cf .ne. 0.0d0) wt = abs(cf) / (abs(cf) + abs(cb))
               ll = l
            end if
!*
            df = xnew(k) - xold(l)
            db = xnew(k) - xold(lm1)
            fnew(k) = (1.0d0 - wt) * (af + (bf + cf * df) * df) + wt * 
     &                (ab + (bb + cb * db) * db)
!*
         else
!*
            if (l .ne. ll) then
               ll = l
               lm1 = l - 1
               af = fold(lm1)
               bf = (fold(l) - fold(lm1)) / (xold(l) - xold(lm1))
            end if
!*
            fnew(k) = af + bf * (xnew(k) - xold(lm1))
         end if
!*
      end do
!*
      map1 = ll - 1
      end function 

      end module 
 
