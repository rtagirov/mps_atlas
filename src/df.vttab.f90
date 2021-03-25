subroutine vttab
!---- turbulent velocity broadening ---- ! 
  use types
  use dfcomm
  implicit none

  include 'common.turbpr'

  logical check
  integer iv, i
  real(kind=dp) v, sum

do iv=1, numvt 
  v = ivt(iv) 
  sum = 1.0d0
  check = .true.
  i = 0
  do while (check .and. i .lt. 120 )
    i = i + 1
    vtprof(i, iv) = exp(-(i / v *2.99792458e5 / 500000.)**2)
    sum = sum + vtprof(i,iv) * 2.
    if (vtprof(i, iv) < 1.e-5) check = .false.
  end do
  nvt(iv) = i
      do i = 1, nvt(iv)
        vtprof(i, iv) = vtprof(i, iv) / sum
      end do
   vtcenter(iv) = 1. / sum
   write(6, 7)iv, vtcenter(iv), (vtprof(i, iv), i = 1, nvt(iv))
7  format(i3,' km/s',10f10.7/(8x,10f10.7))
end do

return
end
