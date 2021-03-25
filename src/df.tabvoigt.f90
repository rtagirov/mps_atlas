subroutine tabvoigt(vsteps, n)
!---- pre-tabulate voight profiles to speed up line profile calculations
  use types
  use dfcomm
  implicit none
  
  integer, intent(in):: n
  real(kind= 4), intent(in) :: vsteps
  integer:: i, ia, idum,  map1, iv 
  double precision, dimension(81) :: tabvi, tabh1
  real(kind=dp):: a, aa, v, vv, vvu, hh1, hh2, u, aau, hh3, hh4

!------------------------------------------------------------------
  tabvi = (/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,   &
           1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9, &
           3.,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6, &
           4.8,5.0,5.2,5.4,5.6, 5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,   &
           7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8, 9.0,9.2,9.4,9.6,9.8,   &
           10.0,10.2,10.4,10.6,10.8,11.0,11.2,11.4,11.6,11.8,12.0/)

  tabh1 = (/-1.12838,-1.10596,-1.04048,-.93703,-.80346,-.64945,     &
           -.48552,-.32192,-.16772,-.03012,.08594,.17789, .24537,  &
            .28981, .31394, .32130,.31573,.30094,.28027, .25648,   &
            .231726,.207528, .184882, .164341,.146128, .130236,    &
            .116515,.104739, .094653,.086005,.078565, .072129,     &
            .066526,.061615, .057281,.053430,.049988, .046894,     &
            .044098, .041561, .039250,.035195,.031762, .028824,    &
            .026288,.024081, .022146, .020441,.018929, .017582,    &
            .016375,.015291, .014312,.013426,.012620, .0118860,    &
            .0112145, .0105990,.0100332,.0095119, .0090306,        &
            .0085852, .0081722,.0077885,.0074314, .0070985,        &
            .0067875, .0064967,.0062243, .0059688, .0057287,       &
            .0055030, .0052903,.0050898,.0049006, .0047217,        &
            .0045526, .0043924,.0042405,.0040964, .0039595/)
!   PRETABULATE VOIGT FUNCTION
!     100 STEPS PER DOPPLER WIDTH GIVES 2 PER CENT ACCURACY
! debug
!  write(6,*)'tabvi=', tabvi
!  write(6,*) 'tabh1=', tabh1

  do i=1, n
    h0tab(i)=dble(i-1) / dble(vsteps) 
  end do

  idum = map1(tabvi, tabh1, 81, h0tab, h1tab, n)

  do i=1, n
    vv=(dble(i - 1) /dble(vsteps))**2
    h0tab(i) = exp(-vv)
    h2tab(i) = h0tab(i) - (vv + vv) * h0tab(i)
  end do


  do ia = 41, 281
    a = dble(ia - 1)/200.0d0 
    aa = a**2
    do iv = 1, 1001
      v = (dble(iv) - 1.0d0) / 200.0d0 
      vv = v**2
      if(a + v .gt.  3.2d0) then
        u =(aa + vv) * 1.4142d0 
        aau = aa / u
        vvu = vv / u
        atab(ia, iv) = a * .79788 / u *((((aau- 10. *vvu)*aau +5.*vvu*vvu) / u+ vvu- aau/3.) * 3./u + 1.)
      else
        hh1 = h1tab(iv) + h0tab(iv) * 1.12838
        hh2 = h2tab(iv) + hh1 * 1.12838 - h0tab(iv)
        hh3 = (1. - h2tab(iv)) * .37613 - hh1 * .66667 * vv + hh2 * 1.12838
        hh4 = (3. * hh3 - hh1) * .37613 + h0tab(iv) * .66667 * vv**2
        atab(ia, iv) = ((((hh4 * a + hh3) * a + hh2) * a + hh1) * a + h0tab(iv)) * &
             (((-.122727278 * a + .532770573) * a - .96284325) * a + .979895032)
      end if
    end do
  end do


  return

end subroutine tabvoigt
