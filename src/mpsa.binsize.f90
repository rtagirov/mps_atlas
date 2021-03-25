 subroutine get_bin_size(nsize, wld, mode)
  use types
  use dfcomm

  implicit none 
  include 'common.ifopbl'
 
  integer mode, i
  integer, intent(out) :: nsize
  real(kind=dp), intent(out) :: wld(mode)

  if (ifkbin ) then 
     call alloc_wave(0, 12) 
  else
     print*, 'Somthing went wrong with the frequency grid!'
     stop
  end if 

  wlend = 0.0d0

  call def_binsize()


  if ( mode .eq. 328 ) then 
   nsize = nsizebig

   do i = 1, nsize
     wld(i) = wavebig(i+1)
   end do 

  else if (mode .eq. 1212) then 

   nsize = nsizelit

   do i = 1, nsize
     wld(i) = wavelit(i+1)
   end do 

  else
 
   print*, 'Something went wrong'
   stop
  end if

 
  call close_wave 
  return
 end

