subroutine yieldf(avect, devia, frict, ncrit, sint3, steff, theta, varj2)
!
! SELECTS YIELD FUNCTION AND CALCULATES VECTOR AVECT
!
    implicit none

    integer :: ncrit, nstr1, istr1

    real :: frict, sint3, steff, theta, varj2, tanth, sinth, costh, root3, cost3, cons1, abthe, cons2, cons3, plumi, tant3, snphi

    real :: AVECT(4) ,devia(4) ,veca1(4) ,veca2(4) ,veca3(4)

    if(steff .eq. 0.0) return
    nstr1=4
    tanth=tan(theta)
    sinth=SIN( theta)
    costh=COS( theta)
    cost3=COS(3.0*theta)
    root3=1.73205080757

    ! calculate vector A1
    veca1(1)=1.0
    veca1(2)=1.0
    veca1(3)=0.0
    veca1(4)=1.0

    ! calculate vector A2
    do istr1=1, nstr1
        veca2( istr1)=devia(istr1)/(2.0*steff)
    end do
    veca2(3)=devia(3)/steff
    
    ! calculate vector A3
    veca3(1)=devia(2)*devia(4)+varj2/3.0
    veca3(2)=devia( 1)*devia(4)+varj2/3.0
    veca3(3)=-2.0*devia(3)*devia(4)
    veca3(4)=devia(1)*devia(2)-devia(3)*devia(3)+varj2/3.0
    GO TO (1,2,3,4) ncrit

    ! TRESCA
1   cons1=0.0
    abthe=abs(theta*0.5729577951308)
    IF(abthe .lt. 29.0) GO TO 20
    cons2=root3
    cons3=0.0
    GO TO 40
20  cons2=2.0*(costh+sinth*TAN(3.0*theta) )
    cons3=root3*sinth/ ( varj2*cost3)
    GO TO 40

    ! VON MISES
2   cons1 = 0.0
    cons2 = root3
    cons3 = 0.0
    goto 40

    ! MOHR-COULOMB
3   cons1=sin(frict*0.017453292)/3.0
    abthe=abs(theta*57.29577951308)
    if(abthe .lt. 29.0) goto 30
    cons3=0.0
    plumi=1.0
    if(theta .gt. 0.0) plumi=-1.0
    cons2=0.5*( root3+plumi*cons1/ root3 )
    goto 40
30  tant3=tan(3.0*theta)
    cons2=costh*((1.0+tanth*tant3)+cons1*(tant3-tanth)/root3)
    cons3=(root3*sinth+cons1*costh)/(2.0*varj2*cost3)
    GO TO 40
    
    ! DRUCKER-PRAGER
4   snphi=sin(frict*0.017453292)
    cons1=2.0*snphi/( root3*(3.0-snphi) )
    cons2=1.0
    cons3=0.0
40  CONTINUE
    
    do istr1=1, nstr1
        avect(istr1)=cons1*veca1(istr1)+cons2*veca2(istr1)+cons3*veca3(istr1)
    end do

    return

    end subroutine
