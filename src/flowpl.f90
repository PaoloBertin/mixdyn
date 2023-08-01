subroutine flowpl(avect, ABETA, dvect, ntype, props, lprop, nstr1, nmats)
!
!
! THIS SUBROUTINE EVALUATES THE PLASTIC D VECTOR
!
!
    implicit none

    include 'param.inc'

    integer :: nstr1, istr1, ntype, lprop, nmats

    real :: young, poiss, hards, fmul1, fmul2, fmul3, denom, abeta
    real :: avect(4), dvect(4), props(mmats, 7)

    young=props(lprop, 1)
    poiss=props(lprop, 2)
    hards=props(lprop, 6)
    fmul1=young/(1.0+poiss)
    if(ntype .eq. 1) goto 60
    fmul2=young*poiss*(avect( 1)+avect(2)+avect( 4))/((1.0+poiss)*(1.0-2.0*poiss))
    dvect(1) =fmul1 *avect (1) +fmul2
    dvect(2) =fmul1 *avect (2) +fmul2
    dvect(3)=0.5*avect(3)*young/(1.0+poiss)
    dvect(4) =fmul1*avect(4)+fmul2
    goto 70

60  fmul3=young*poiss*( avect( 1)+avect(2))/(1.0-poiss*poiss)
    dvect(1)=fmul1*avect(1)+fmul3
    dvect(2)=fmul1*avect (2) +fmul3
    dvect(3)=0.5*avect(3)*young/(1.0+poiss)
    dvect(4)=fmul1*avect(4)+fmul3

70  denom=hards
    do  istr1=1,nstr1
        denom=denom+avect( istr1)*dvect( istr1)
    end do

    abeta=1.0/denom

    RETURN

END


