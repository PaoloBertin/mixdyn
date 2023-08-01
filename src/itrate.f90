subroutine itrate(accei, accel ,consd ,consf ,xmass ,dispi , resid, stifs ,toler ,veloi ,velol ,velot , iiter ,miter )
!
! CALCULATES INCREMENT IN DISPLACEMENT AND APPLIES CONVERGENCE
!
    DIMENSION dispi(1) ,veloi(1) ,accei(1) ,resid(1) ,MAXAI(1) , DISPL(1) ,velol(1) ,accel(1) ,stifs(1) ,disp(1) , xmass(1) , &
        velot(1)


    nchek=0
    call multpy (accel ,xmass ,accei ,MAXAI ,nsize ,NWMTL )

! Calculates total effective load vector
    do isize =1, nsize
        accel(isize) = displ(isize) - accel(isize) - resid(isize)
    end do

! Calculates delta displacement
    call redbak(stifs, accel, maxai, nsize)

! Applies convergence
    sumpp=0.0
    sumpq=0.0
    do isize=1,nsize
        dispp = accel(isize)
        dispq = disp( isize) + dispp
        disp(isize) = dispq
        sumpp = sumpp + dispp * dispp
        sumpq = sumpq + dispq*dispq
    end do

    DO isize=1,nsize
        accei( isize) =consf*(disp ( isize) -dispi( isize) )
        velot( isize) =veloi( isize) +consd*accei( isize)
    end do

    sumpp = sqrt(sumpp/sumpq)
    if(sumpp .gt. toler) goto 550
    nchek=1
    goto 240
550 if(iiter .lt. miter) goto 230
240 do  isize=1,nsize
        veloi( isize) =velot( isize)
        dispi( isize)=disp(iSIZE)
    end do
230 continue

    return

end subroutine
