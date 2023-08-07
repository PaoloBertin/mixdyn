program mixdyn
! **********************************************************************
!
! TIME INTEGRATION IMPLICIT-EXPLICITY ALGORITHM
!
! Variables
!
!   NPOIN   Number of nodes
!   NELEM   Number of elements
!   NDOFN   Number of degrees of freedom per node
!   NVFIX   Number of bound nodes
!   NCONM   Number of concentrated masses
!   NRADS
!
!   COORD   Nodal coordinates
!   PROPS
!   NPOIN
! **********************************************************************
    use model

    implicit none

    integer :: istep, iiter, nchek

    real :: consd, consf

    integer :: dt(8)

    ! real :: dispi(mpoin*mdime), veloi(mpoin*mdime), accei(mpoin*mdime),                  &
    !     displ(mpoin*mdime), ymass(mpoin*mdime), accek(mpoin*mdime), accej(mpoin*mdime), resid(mpoin*mdime), dispt(mpoin*mdime),    &
    !     dispq(mpoin*mdime), tempe(mpoin), force(mpoin*mdime), acceh(mpoin*mdime), velot(mpoin*mdime), velol(mpoin*mdime),          &
    !     accel(mpoin*mdime), accev(mpoin*mdime), stiff(maxma), stifs(maxma), stifi(maxma), xmass(maxma), dampi(maxma),              &
    !     dampg(maxma), rload(melem, mnode*mdime), strin(4, 450), strag(4, 450), strsg(4, 450), epstn(450), effst(450), posgp(4),    &
    !     weigp(4)

    ! common stiff, xmass, dampg, stifi, stifs, dampi

    ! Start
    call date_and_time(values=dt)
    write(*, '(a7, i4, 5(a, i2.2), a, i3/)') 'Start ', dt(1), '-', dt(2), '-', dt(3), ' ', dt(5), ':', dt(6), ':', dt(7), '.', dt(8)

    ! Open the file I/O
    open(5, file='/home/paolo-bertin/Documenti/FortranProjets/mixdyn/data/input.dat',  status='OLD')
    open(6, file='/home/paolo-bertin/Documenti/FortranProjets/mixdyn/data/result.dat', status='OLD')

    call date_and_time(values=dt)
    write(*, '(a7, i4, 5(a, i2.2), a, i3/)') 'Start ', dt(1), '-', dt(2), '-', dt(3), ' ', dt(5), ':', dt(6), ':', dt(7), '.', dt(8)

    call contol()
    call inputd()
    call intime()
    goto 100
    call prevos()
    call loadpl()
    call lumass()
    call linkin()
    do istep = 1, nstep
        do iiter=1, miter
            call gstiff(istep)
            call impexp(consd, consf, iiter, istep)
            call resepl(rload, iiter, resid, istep)
            call itrate(consd, consf, iiter)
            if(nchek .eq. 1) exit
        end do
        call outdyn(iiter, istep)
    end do

!   Close file I/O
    close(5)
    close(6)

100 continue
    stop

end program mixdyn
