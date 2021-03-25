subroutine str(i, s, n)
    implicit none
    character(len=n), intent(out) :: s
    integer, intent(in) :: i,n
    write (s, *) i
    s = adjustl(s)
    s = trim(s)
end subroutine

subroutine convits(i, s)
     implicit none
     character(len=6), intent(out) ::s
     integer, intent(in) :: i
     character(len=1) :: units, tens, hun, taus, ttaus, htaus

     htaus = char(i/100000+48)
     ttaus = char(mod(i,100000)/10000 +48)
     taus  = char(mod(i,10000)/1000 +48)
     
     hun   = char(mod(i,1000)/100 +48)
     tens  = char(mod(i,100)/10 + 48)
     units = char(mod(i,10) +48)

     s = htaus//ttaus//taus//hun//tens//units

     


end subroutine 
