subroutine modps(dmatx, lprop, nmats, nstre, ntype, props)
! **********************************************************************
!
! *** THIS SUBROUTINE EVALUATES THE D-MATRIX
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: lprop, nmats, nstre, ntype, istre, jstre

    real :: young, poiss, const, conss
    real :: props(nmats, 1), dmatx(4,4)


    young = props(lprop, 1)
    poiss = props(lprop, 2)

    do istre=1, 4
        do jstre =1, 4
            dmatx(ISTRE,JSTRE)=0.0
        end do
    end do

! matrix for plane stress case
    if(nstre .eq. 1) then
        const = young/ (1.0-poiss*poiss)
        dmatx(1,1) = const
        dmatx(2,2) = const
        dmatx(1,2) = const*poiss
        dmatx(2,1) = const*poiss
        dmatx(3,3) = (1.0-poiss)*const/2.0
        return

    endif

! matrix for plane strain case
    if(nstre .eq. 2) then
        const =young*(1.0-poiss)/((1.0+poiss)*(1.0-2.0*poiss))
        dmatx(1,1) = const
        dmatx(2,2) = const
        dmatx(1,2) = const*poiss/(1.0-poiss)
        dmatx(2,1) = const*poiss/(1.0-poiss)
        dmatx(3,3) = (1.0-2.0*poiss)*const/(2.0*(1.0-poiss))
        return

    endif

! matrix for axisymmetricx case
    if(nstre .eq. 3) then
        const=young*(1.0-poiss)/((1.0+poiss) *(1.0-2.0*poiss))
        conss=poiss/(1.0-poiss)
        dmatx(1,1) = const
        dmatx(2,2) = const
        dmatx(3,3) = const*(1.0-2.0*poiss)/(2.0*(1.0-poiss))
        dmatx(1,2) = const*conss
        dmatx(1,4) = const*conss
        dmatx(2,1) = const*conss
        dmatx(2,4) = const*conss
        dmatx(4,1) = const*conss
        dmatx(4,2) = const*conss
        dmatx(4,4) = const
        return
    endif

    return
end
