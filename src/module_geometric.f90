module geometry
    implicit none
    integer, parameter :: mpoin=200, mdime=2, melem=50, mnode=9, mprop=7, mmats=10

    integer :: lnods(melem, mnode), matno(melem)

    real :: coord(mpoin, mdime), props(mmats, mprop)
    
contains
    
end module geometry