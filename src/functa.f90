real function functa (accer, jstep)
!
!   Accelerogram interpolation
!
    use model
    
    implicit none
    
    integer, intent(in) :: jstep
    real,    intent(in) :: accer(1)
    
    integer :: ngash, mgash

    real :: xgash

    if(ifunc .ne. 0) return

    functa = 0.0

    if(jstep .eq. 0 .or. float(jstep)*dtime .gt. dtend) return
    xgash = (float(jstep) - 1.0)/afact + 1.0

    mgash = xgash
    ngash = mgash + 1

    xgash=xgash - float(mgash)

    functa = accer(mgash) * (1.0 - xgash) + xgash*accer(ngash)

    return

end function

