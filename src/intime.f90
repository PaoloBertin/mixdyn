subroutine intime(aalfa, acceh, accev, afact, azero, beeta, bzero, delta, dtime, dtend, gaama, ifixd, ifunc, intgr, kstep, miter,  &
    ndofn, nelem, ngrqs, noutd, noutp, npoin, nprqd, nreqd, nreqs, nstep, omega, dispi, toler, veloi, ipred)
! **********************************************************************
!
! INITIAL VALUES AND TIME 1NTEGRATION DATA
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: ifixd, ifunc, kstep, miter, ndofn, nelem, noutd, noutp, npoin, nreqd, nreqs, nstep, ipred, nacce, ireqd, ireqs,     &
        ielem, jpoin, ipoin, ngash, nposn, iacce, idofn

    integer :: intgr(melem), nprqd(mnode), ngrqs(mnode)

    real :: aalfa, afact, azero, beeta, bzero, delta, dtime, dtend, gaama, omega, toler
    real :: dtrec, xgash, ygash, dispi(mpoin*mdime), veloi(mpoin*mdime), acceh(mpoin*mdime), accev(mpoin*mdime)

! READ TIME STEPPING AND SELECTIVE OUTPUT PARAMETERS
    read(5, 900)
    read(5, 902) nstep, noutd, noutp, nreqd, nreqs, nacce, ifunc, ifixd, miter, kstep, ipred
    write(6, 950) nstep, noutd, noutp, nreqd, nreqs, nacce, ifunc, ifixd, miter, kstep, ipred

    read(5, 900)
    read(5, 190) dtime, dtend, dtrec, aalfa, beeta, delta, gaama, azero, bzero, omega, toler
    write(6, 960) dtime, dtend, dtrec, aalfa, beeta, delta, gaama, azero, bzero, omega, toler

! SELECTED NODES AND GAUSS POINTS FOR OUTPUT
    read(5, 900)
    read(5, 902) (nprqd(ireqd), ireqd=1, nreqd)
    read(5, 900)
    read(5, 902) (ngrqs(ireqs), ireqs=1, nreqs)
    write(6, 909)
    write(6, 910) (nprqd(ireqd), ireqd=1, nreqd)
    write(6, 911) (ngrqs(ireqs), ireqs=1, nreqs)

! READ THE INDICATOR FOR EXPLICIT OR IMPLICIT ELEMENT
    read(5,900)
    read(5, 902) (intgr(ielem), ielem=1 , nelem)
    write(6, 930)
    write(6, 902)(intgr(ielem), ielem=1 , nelem)

! INITIAL DISPLACEMENTS
    jpoin=0
    do ipoin=1, npoin
        do idofn =1, ndofn
            jpoin = jpoin + 1
            dispi(jpoin) = 0.0
            veloi(jpoin) = 0.0
        end do
    end do

    write(6, 903)
    read(5,900)
    do while(ngash .ne. npoin)
        read(5, 904) ngash, xgash, ygash
        nposn = (ngash - 1)*ndofn + 1
        dispi(nposn) = xgash
        nposn = nposn + 1
        dispi(nposn) = ygash
        write(6, 904) ngash, xgash, ygash
    end do

! INITIAL VELOCITIES
    write(6, 906)
    write(5,  900)
    do while(ngash .ne. npoin)
        read(5,  904) ngash, xgash , ygash
        nposn = (ngash - 1)*ndofn + 1
        veloi(nposn) = xgash
        nposn = nposn + 1
        veloi(nposn) = ygash
        write(6, 904) ngash, xgash, ygash
    end do

    if(ifunc .ne. 0) goto 250

! READ ACCELEREGRAM DATA , X-DIREC FROM TAPE 7 ,Y-DIREC FROM TAPE 12
    afact = dtrec/dtime
    if(ifixd .eq. 0) then
!        OPEN (7,  FILE='accex.dat', STATUS='OLD')
!        OPEN (12, FILE='accey.dat', STATUS='OLD')
        read(5, 907) (acceh(iacce), iacce=1, nacce)
        read(5, 907) (accev(iacce), iacce=1, nacce)
        write(6, 913) dtrec
        write(6, 907) (accev(iacce) , iacce=1, nacce)

        close(7)
        close(12)
    endif

    if(ifixd .eq. 1) then
        open(12, FILE='accex.dat', STATUS='OLD')
        read(12 , 907) (ACCEV(IACCE), IACCE=1, nacce)
        write(6 ,913) dtrec
        write(6, 907) (ACCEV(IACCE), IACCE=1 , nacce)
        close(12)
    endif

    if(ifixd .eq. 2) then
        open(12, FILE='accey.dat', STATUS='OLD')
        read(7, 907) (acceh(iacce), iacce=1, nacce)
        write(6 ,912)
        write(6, 907) (acceh(iacce), iacce=1, nacce)
        close(12)
    endif

250 continue

900 format()
902 format(16I5)
190 format(11F10.4)
950 format(//'T1ME STEPP1NG PARAMETERS',/ &
        'nstep=', I5, 5X, 'noutd=', I5, 5X, 'noutp=', I5,/ &
        'nreqd=', I5, 5X, 'nreqs=', I5, 5X, 'nacce=', I5,/ &
        'ifunc=', I5, 5X, 'ifixd=', I5, 5X, 'M1TER=', I5,/ &
        'kstep=', I5, 5X, 'ipred=', I5)
960 format(//, &
        'dtime = ', F10.7, 5X, ' dtend = ', F10.7,/ &
        'dtrec = ', F10.7, 5X, ' aalfa = ', F10.7,/ &
        'beeta = ', F10.7, 5X, ' delta = ', F10.7,/ &
        'gaama = ', F10.7, 5X, ' azero = ', F10.7,/ &
        'bzero = ', F10.7, 5X, ' omega = ', F10.7,/ &
        'toler = ', F10.7)
909 format(//'SELECTIVE OUTPUT REQUESTED FOR FOLLOWING')
910 format('NODES', 10I5)
911 format('G.P. ', 10I5)
930 format(//'TYPE OF ELEMENT , IMPLICIT=1 , EXPLICIT=2 ')
903 format(//' NODE', 7X, 'INITIAL X-DISP', 2X, 'INITIAL Y-DISP')
904 format(I5, F10.5, 5X, F10.4)
!905   format(I5, F10.4)
906 format(//' NODE ', 2X, 'INITIAL X-VEL.', 2X, 'INITIAL Y-VEL.')
907 format(7F10.3)
912 format('HORIZONTAL ACCELERATION ORDINATES AT', F9.4, 2X, 'SEC'/ )
913 format('VERTICAL ACCELERATION ORDINATES AT  ', F9.4, 2X, 'SEC'/ )

    RETURN
END
