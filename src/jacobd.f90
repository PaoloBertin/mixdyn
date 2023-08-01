subroutine jacobd(cartd, dlcod, djacm, ndime, nlaps, nnode)
! **********************************************************************
!
! DEFORMATION JACOBIAN
!
! **********************************************************************
    implicit none

    integer :: idime, ndime, jdime, inode, nlaps, nnode

    real :: cartd(2,9), dlcod(2,9), djacm(2,2)

    if(nlaps .gt. 1) goto 10
! for small displacement
    djacm(1, 1)=1.0
    djacm(2, 2)=1.0
    djacm(1, 2)=0.0
    djacm(2, 1)=0.0

! for large displacement
10  continue
    do 20 idime = 1, ndime
        do 25 jdime =1 , ndime
            djacm(idime, jdime) = 0.0
            do 30 inode=1, nnode
                djacm(idime, jdime) = djacm(idime, jdime) + DLCOD(idime, inode) * CARTD(jdime , inode)
30          continue
25      continue
20  continue

    return
end subroutine jacobd
