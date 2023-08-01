subroutine gaussq(ngaus, posgp, weigp)
! **********************************************************************
!
!  MIXDYN INPUT ROUTINE
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: ngaus,kgaus, igash, jgash

    real :: posgp(4), weigp(4)

    if(ngaus .le. 2) then
        posgp(1)=-0.577350269189626
        weigp(1)= 1.0
    else
        posgp(1)=-0.774596669241483
        posgp(2)=0.0
        weigp(1)=0.555555555555556
        weigp(2)=0.888888888888889
    endif

    kgaus=ngaus/2
    do igash = 1, kgaus
        jgash = ngaus + 1 - igash
        posgp(jgash) = -posgp(igash)
        weigp(jgash) = weigp(igash)
    end do

    return
end subroutine
