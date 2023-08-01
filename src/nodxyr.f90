subroutine nodxyr(coord, lnods, nelem, nnode, npoin, nrads, ntype)
! **********************************************************************
!
! INTERPOLATION OF MIDSIDE AND CENTER NODES
!
! **********************************************************************
implicit none

include 'param.inc'

integer :: ipoin, npoin, ntype, nelem, nnode, nrads, lnode, ielem, inode, nodst, igash, ndofn, midpt, nodmd, kount
integer :: lnods(melem, mnode)

real :: raddi, theta, total
real :: coord(mpoin, mdime)

if((ntype .ne. 3 .or. nrads .eq. 0) .and. nnode .eq. 4) return

! CHANGE POLAR COORDINATES TO CARTISIAN
do ipoin =1, npoin
  raddi = coord(ipoin, 1)
  theta = coord(ipoin, 2)
  theta = 0.017453292*theta
  coord(ipoin, 1) = raddi* sin(theta)
  coord(ipoin, 2) = raddi*cos(theta) 
end do

lnode = nnode - 1
do ielem=1, nelem
  do inode =1, nnode
    if(inode .eq. 9) cycle
!   COMPUTE THE NODE NUMBER OF THE FIRST NUMBER
    nodst = lnods(ielem, inode)
    igash = inode + 2
    if(igash .gt. lnode) igash=1

!   COMPUTE THE NODE NUMBER OF THE LAST NODE
    ndofn = lnods(ielem, igash)
    midpt = inode + 1

!   COMPUTE THE NODE NUMBER OF THE INTERMEDIATE NODE
    nodmd = lnods(IELEM, MIDPT)
    total = abs(coord(nodmd, 1) + coord(nodmd, 2))

!   IF THE COORDINATES OF THE INTERMEDIATE NODE ARE BOTH ZERO INTERPOLATE BY A STRAIGHT LINE
    if(total .gt. 0.0) cycle
    kount=1
10  coord(nodmd, kount)=(coord(nodst, kount) + COORD(ndofn, kount))/2.0 
    kount = kount + 1         
    if(kount .eq. 2) goto 10
  end do
end do

return
end
