subroutine invar(devia, lprop, ncrit, nmats, props, sint3, steff, stemp, theta, varj2, yield)
!
! STRESS INVARIANTS
!
    implicit none

    integer :: ncrit, nmats
    integer :: lprop

    real :: root3, smean, varj2, sint3, theta, varj3, yield, phira, snphi, steff
    real :: devia(4), stemp(4), props(nmats, 1)

    ! Invariants
    root3=1.73205080757
    smean=(stemp(1)+stemp(2)+stemp(4) )/3.0
    devia(1) = stemp(1)-smean
    devia(2) = stemp(2)-smean
    devia(3) = stemp(3)
    devia(4) = stemp(4)-smean
    varj2=devia(3)*devia(3)+0.5*(devia(1) *devia(1)+ devia(2)*devia(2)+devia(4) *devia(4))
    varj3=devia(4)*(devia(4) *devia(4)-varj2)
    steff =sqrt(varj2)
    if(varj2 .eq. 0.0 .or. steff .eq. 0.0) goto 5
    sint3=-2.5980762113*varj3/(varj2*steff)
    goto 6
5   sint3=0.0
6   continue
    IF(sint3.LT.-1.0) sint3 =-1.0

    IF(sint3.GT. 1.0) sint3 = 1.0
    theta=asin(sint3)/3.0
    goto (1,2,3,4) NCRIT

    ! TRESCA
1   yield=2.0*COS( theta) *steff
    return

    ! VON MISES
2   yield=root3*steff
    return

    ! MOHR-COULOMB
3   phira=props(LPROP, 8)*0.017453292
    snphi=SIN(phira)
    yield=smean*snphi+steff*(cos(theta)-sin(theta)*snphi/root3)
    return

    ! DRUCKER-PRAGER
4   phira=props(LPROP, 8)*0.017453292
    snphi=SIN(phira)
    yield=6.0*smean*snphi/(root3*(3.0-snphi))+steff
    return
    
end subroutine
