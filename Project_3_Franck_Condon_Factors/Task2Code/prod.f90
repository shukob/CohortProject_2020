subroutine dmatvecprod(prod, arg1, arg2, size)
    implicit none
    integer size, i
    real(kind=8) prod(size), arg1(size), arg2(size)

    do i = 1, size
        prod(i) = arg1(i) * arg2(i)
    end do
    return
end
subroutine dotprod(arg1, arg2, size, dot)
    implicit none
    integer size, i
    real(kind=8) dot, arg1(size), arg2(size), ddot
    !      dot=ddot(size,arg1,1,arg2,1)
    dot = 0.d0
    do i = 1, size
        dot = dot + arg1(i) * arg2(i)
    end do
    return
end
subroutine cdotprod(arg1, arg2, size, dot)
    implicit none
    integer size, i
    complex(kind=8) dot, arg1(size), arg2(size), zdotu
    !      dot=zdotu(size,arg1,1,arg2,1)
    dot = (0.d0, 0.d0)
    do i = 1, size
        dot = dot + dconjg(arg1(i)) * arg2(i)
    end do
    return
end