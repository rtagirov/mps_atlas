subroutine bisort_filter(k, n, k1)
    use types
    use dfcomm
    implicit none
    integer n, i, j, j0, j1, jj0, mask, max_val, k(n), k1(n)
    double precision indexx1(lenbuff)

    do j = 1, n
        k(j) = k(j) + 30000
    end do

    mask = 1
    max_val = k(1)
    do j = 2, n
        max_val = max(max_val, k(j))
    end do

    do i = 1, 59
        if (i .ne. 1) then
            max_val = max_val / 2
            mask = mask * 2
        end if

        if (max_val .eq.  0) then
            do j = 1, n
                k(j) = k(j) - 30000
            end do
            return 
        endif

        j1 = 0
        j0 = 0
        do j = 1, n
! --- this condition first takes the bit-wise intersection of the mask and the number, and than compares if it less equal 0

            if ((k(j) .and. mask) .le. 0)  then
                j0 = j0 + 1
                k(j0) = k(j)
                indexx(j0)=indexx(j)


            else
               j1 = j1 + 1
                k1(j1) = k(j)
                indexx1(j1)=indexx(j)

            end if
        end do

        if(j1 .ne.  0) then
            do j = 1, j1
                jj0 = j + j0
                k(jj0) = k1(j)
                indexx(jj0)=indexx1(j)
            end do
        end if
    end do
end


