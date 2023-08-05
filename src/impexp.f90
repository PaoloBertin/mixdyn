subroutine impexp (consd, consf, iiter ,istep)
!
! Generates partial effective load vector
!
    use model

    implicit none

    integer :: istep, kount, ishot, ipoin, idofn, isize, iiter, nsize, imaxa, kmaxa, jmaxa, lmaxa, iwktl, iwmtl

    real :: consa, consd, consf, consb, consc, consg, consh, conse, facts, facth, factv
    real :: functa, functs

    if(istep .gt. 1 .or. iiter .gt. 1) goto 1000

    consa = dtime * dtime * (0.5 - delta)
    consb = dtime * (1.0 - gaama)
    consc = dtime * dtime * delta
    consd = dtime * gaama
    consf = 1./consc
    consg = beeta * gaama * dtime
    consh = aalfa * gaama * dtime
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


    do isize =1, nsize
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
        ! accei(isize) = rload(isize) - displ(isize) - velol(isize)
        accei(isize) = force(isize) - displ(isize) - velol(isize)
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

210     dispi(isize)=dispt(isize) + dtime*veloi (isize) +consa * accei( isize)
        veloi ( isize) =veloi (isize)+CONSB*accei(isize)
        if(ipred .eq. 2) cycle
        dispt(isize)=dispi(isize)

        velot(isize) = veloi (isize)
        accei(isize) = consf * (dispt(isize) - dispi(isize) )
    end do

    ! Calculates load vectors
    facts = functs(istep)
    facth = functa(acceh, istep)
    factv = functa(accev, istep)
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
        !if(ifixd .eq. 2) displ(isize) = -velol(isize) - facth*accej(isize) + rload(isize)
        !if(ifixd .eq. 1) displ(isize) = -velol(isize) - factv*accek(isize) + rload(isize)
        !if(ifixd .eq. 0) displ(isize) = -velol(isize) - facth*accej(isize) + rload(isize) - factv*accek(isize)
        if(ifixd .eq. 2) displ(isize) = -velol(isize) - facth*accej(isize) + force(isize)
        if(ifixd .eq. 1) displ(isize) = -velol(isize) - factv*accek(isize) + force(isize)
        if(ifixd .eq. 0) displ(isize) = -velol(isize) - facth*accej(isize) + force(isize) - factv*accek(isize)
        if(ifunc .eq. 0) cycle
!570     displ(isize) = -velol(isize) + rload(isize)*facts
570     displ(isize) = -velol(isize) + force(isize)*facts
    end do

    return

900 format(/' Initial acceleration '/)
910 format( 1X, 10E12.5)

end subroutine
