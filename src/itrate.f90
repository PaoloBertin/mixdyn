subroutine itrate(consd, consf, iiter)
!
! Calculates increment in displacement and applies convergence
!
    use model

    implicit none

    integer :: iiter, nchek, nsize, isize
    real :: consd, consf, sumpp, sumpq, dispp, dispq

    nchek=0
    call multpy(accel, xmass, accei, maxai, nsize, nwmtl)

    ! Calculates total effective load vector
    do isize =1, nsize
        accel(isize) = displ(isize) - accel(isize) - resid(isize)
    end do

    ! Calculates delta displacement
    call redbak(stifs, accel, maxai, nsize)

    ! Applies convergence
    sumpp = 0.0
    sumpq = 0.0
    do isize = 1, nsize
        dispp = accel(isize)
        dispq = dispt(isize) + dispp
        dispt(isize) = dispq
        sumpp = sumpp + dispp * dispp
        sumpq = sumpq + dispq*dispq
    end do

    do isize = 1, nsize
        accei(isize) = consf * (dispt(isize) - dispi(isize) )
        velot(isize) = veloi(isize) + consd * accei(isize)
    end do

    sumpp = sqrt(sumpp/sumpq)
    if(sumpp .gt. toler) goto 550
    nchek = 1
    goto 240

550 if(iiter .lt. miter) goto 230

240 do  isize = 1,nsize
        veloi(isize) = velot( isize)
        dispi(isize) = dispt(isize)
    end do
230 continue

    return

end subroutine itrate
