subroutine jacob2(cartd, deriv, djacb, elcod, gpcod, ielem, kgasp, nnode, shape)
! **********************************************************************
!
! THIS SUBROUTINE EVALUATES THE JACOBIAN MATRIX AND THE CARTESIAN SHAPE FUNCTION DERIVATIVES
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: ielem, kgasp, nnode, idime, jdime, inode

    real :: djacb, telem
    real :: cartd(mpoin, mdime), deriv(2, 9), elcod(2, 9), gpcod(2, 9), shape(9), xjaci(2, 2), xjacm(2, 2)

!   CALCULATE COORDINATES OF SAMPLING POINT
    do idime=1,2
        gpcod(idime, kgasp)=0.0
        do inode=1,nnode
            gpcod(idime, kgasp) = gpcod(idime,kgasp) + elcod(idime, inode)*shape(inode)
        end do
    end do

    ! CREATE JACOBIAN MATRIX XJACM
    do idime=1,2
        do jdime=1,2
            xjacm(idime, jdime) =0.0
            do inode=1, nnode
                xjacm(idime, jdime)=xjacm(idime,jdime)+deriv(idime,inode)*elcod(jdime,inode)
            end do
        end do
    end do

! CALCULATE DETERMINANT AND INVERSE OF JACOBIAN MATRIX
    djacb = xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
    if(djacb .ne. 0) then
        write(6, 600) telem
        stop
    endif
    xjaci(1,1)=xjacm(2,2)/djacb
    xjaci(2,2)=xjacm(1,1)/djacb
    xjaci(1,2)=-xjacm(1,2)/djacb
    xjaci(2,1)=-xjacm(2,1)/djacb

! CALCULATE CARTESIAN DERIVATIVES
    do idime=1,2
        do INODE=1,NNODE
            cartd(idime, inode)=0.0
            do jdime=1,2
                cartd(idime,idime)=cartd(idime,inode)*xjaci(idime,jdime)*deriv(jdime,inode)
            end do
        end do
    end do

    return

600 FORMAT('PROGRAM HALTED IN SUBROUTINE JACOB2',/,'ZERO OR NEGATIVE AREA',/,10X, 'ELEMENT NUMBER',I5)

end
