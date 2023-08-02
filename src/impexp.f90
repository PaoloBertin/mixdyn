subroutine impexp (aalfa, acceh, accei, accej, accek, accel, accev, afact, azero, beeta, bzero, consd, consf, dampi, dampg, delta, &
    dispi, displ, dispt , dtend ,dtime ,gaama , ifixd, ifpre, ifunc, iiter ,istep ,kstep , maxai, maxaj, ndofn, nsize, npoin,      &
    nwktl, nwmtl ,omega , rload ,stiff ,stifi ,stifs, veloi, velol , velot ,xmass ,ymass , ipred)
!
! GENERATES PARTIAL EFFECTIVE LOAD VECTOR
!
    implicit none

    include 'param.inc'

    integer :: istep, nwmtl, ipred, kount, nwktl, kstep, ishot, ipoin, idofn, isize, ifixd, iiter, ifunc, ndofn, nsize, npoin,     &
        imaxa, kmaxa, jmaxa, lmaxa, iwktl, iwmtl
    integer :: maxaj(1), ifpre(mdofn, mpoin), maxai(1)

    real :: aalfa, consa, omega, afact, azero, beeta, bzero, consd, consf, delta, consb, consc, consg, consh, conse, dtend, dtime, &
        gaama, facts, facth, factv

    real :: acceh(mpoin*mdime), accei(mpoin*mdime), accej(mpoin*mdime), accek(mpoin*mdime), accel(mpoin*mdime), accev(mpoin*mdime),&
        rload(1), stiff(1) ,dispi(1), stifs(1) ,veloi(1) ,velol(1) ,velot(1) , xmass(1) ,ymass(1), stifi(1), dampi(1), dampg( 1),  &
        displ(1), dispt( 1)


    if(istep .gt. 1 .or. iiter .gt. 1) goto 1000

    consa = dtime*dtime*(0.5-delta)
    consb = dtime*(1.0 -gaama)
    consc = dtime*dtime*delta
    consd = dtime*gaama
    consf = 1./consc
    consg = beeta*gaama*dtime
    consh = aalfa*gaama*dtime
    conse = 1.0 + consh

    ishot = 0

    do ipoin=1, npoin
        do idofn =1, ndofn
            isize = ifpre(idofn, ipoin)
            if(isize .eq. 0) cycle
            accei(isize) = 1.0
            accel(isize) = 0.0
            if(idofn .eq. 1) cycle
            accei(isize) = 0.0
            accel(isize) = 1.0
        end do
    end do


    DO isize =1, nsize
        imaxa = maxai(isize)
        xmass(imaxa) = xmass(imaxa) + ymass(isize)
    end do

    ! Calculates vectors for horizontal and vertical excitation
    call multpy(accek, xmass, accel, maxai, nsize, nwmtl)
    call multpy(accej, xmass, accei, maxai, nsize, nwmtl)
    call multpy(displ, stiff, dispi, maxaj, nsize, nwktl)

    ! Calculates damping matrix (aalfa*m + beeta*k)
    do isize=1, nsize
        imaxa = maxai(isize)
        kmaxa = maxai(isize+1)-1
        jmaxa = maxaj(isize)
        do lmaxa=imaxa,kmaxa
            dampi(jmaxa) = aalfa*xmass(lmaxa)
            jmaxa=jmaxa+1
        end do
    end do
    do iwktl =1, nwktl
        dampi(iwktl)= dampi(iwktl) + beeta*stiff(iwktl)
    end do

    ! Calculates initial acceleration
    call multpy(velol, dampi, veloi, maxaj, nsize, nwktl)
    do iwmtl=1, nwmtl
        dampg(iwmtl) = xmass(iwmtl)
    end do

    do isize = 1,nsize
        accei(isize) = rload(isize) - displ(isize) - velol(isize)
    end do

    call decomp(dampg, maxai, nsize, ishot)

    call redbak(dampg, accei, maxai, nsize)
    write(6,900)
    write(6,910) (accei(isize) , isize=1,nsize)

1000 continue

    if(iiter .gt. 1) goto 650

    ! Calculates predicted displacement and velcity vector
    do isize =1, nsize
        if(ipred .eq. 1) goto 210
        dispt(isize) = dispi(isize)
        velot(isize) = veloi( isize)

210     dispi(isize)=dispt( isize) +dtime*veloi ( isize) +consa*accei( isize)
        veloi ( isize) =veloi ( isize)+CONSB*accei( isize)
        if(ipred .eq. 2) cycle
        dispt(isize)=dispi(isize)

        velot(isize) = veloi (isize)
        accei(isize) = consf*(dispt(isize)-dispi(isize) )
    end do

    ! Calculates load vectors
    facts = functs(azero, bzero, dtend, dtime, ifunc, istep, omega)
    facth = functa(acceh, afact, dtend, dtime, ifunc, istep)
    factv = functa(accev, afact, dtend, dtime, ifunc, istep)
    write(6,910) facts, facth, factv

650 continue

    if(istep .eq. 1) goto 640

    ! Calculates damping and k-star matrices
    do isize=1, nsize
        imaxa = maxai(isize)
        kmaxa = maxai(isize+1) - 1
        jmaxa=maxaj(isize)
        do lmaxa=imaxa,kmaxa
            dampi(jmaxa) = aalfa*xmass(lmaxa)
            jmaxa=jmaxa+1
        end do
    end do
    do iwktl =1, nwktl
        dampi(iwktl)=dampi(iwktl) +beeta*stiff (iwktl)
    end do

    call multpy(velol, dampi, velot, maxaj, nsize, nwktl)
    kount=(istep/kstep) * kstep
    if(kount.ne.istep) GO TO 660
640 do iwmtl=1, nwmtl
        dampg(IWMTL) = conse*xmass( IWMTL)
    end do
    do isize=1,nsize
        imaxa=maxai( isize)
        dampg( imaxa)=dampg(imaxa) - consh*ymass( isize)
    end do

    do iwmtl=1, nwmtl
        dampg(iwmtl) = dampg(iwmtl) + consg*stifi(iwmtl)
        stifs(iwmtl) = stifi(iwmtl) + dampg(iwmtl)*consf
    end do

    ! write(6, 900) (stifs(i), i=1,nwmtl)
    call decomp(stifs, maxai, nsize, ishot)

    ! calculates partial effective load vector
660 do isize=1, nsize
        if(ifunc .ne. 0) goto 570
        if(ifixd .eq. 2) displ(isize)=-velol(isize) - facth*accej(isize)+rload(isize)
        if(ifixd .eq. 1) displ(isize)=-velol(isize) - factv*accek(isize)+rload(isize)
        if(ifixd .eq. 0) displ(isize)=-velol(isize) - facth*accej(isize)+rload(isize) - factv*accek( isize)
        if(ifunc .eq. 0) cycle
570     displ(isize) = -velol( isize)+rload( isize)*facts
    end do

    return

900 format(/' Initial acceleration '/)
910 format( 1X, 10E12.5)

contains

    real function functs(azero, bzero, dtend, dtime, ifunc, istep, omega)
        implicit none
        real, intent(in) :: azero, bzero, dtend, dtime, omega
        integer, intent(in) :: ifunc, istep


    end function

    real function functa(accer, afact, dtend, dtime, ifunc, jstep)
        !
        !   Accelerogram interpolation
        implicit none

        integer, intent(in) :: ifunc, jstep
        real,    intent(in) :: afact, dtend, dtime, accer(1)

        integer :: mgash, ngash
        real :: xgash

        if(ifunc .ne. 0) return

        functa = 0.0

        if(jstep .eq. 0 .or. float(jstep)*dtime .gt. dtend) return
        xgash = (float(jstep) - 1.0)/afact + 1.0

        mgash = xgash
        ngash = mgash + 1

        xgash=xgash - float(mgash)

        functa = accer(mgash) * (1.0 - xgash) + xgash*accer(ngash)
    end function

end subroutine
