subroutine sfr2(deriv, etasp, exisp, nnode, shape)
! **********************************************************************
!
!      THIS SUBROUTINE EVALUATES shape FUNCTIONS AND THEIR derivATIVES
!      FOR LINEAR,QUADRATIC LAGRANGIAN AND SERENDIPITY
!      ISOPARAMETRIC 2-D ELEMENTS
!
! ***********************************************************************
    implicit none

    integer :: nnode

    real :: deriv(2,9), etasp, exisp, shape(9), s, t, st, s1, s2, s9, t9, t1, t2, tt, ss, sst, stt, st2

    s=exisp
    t=etasp
    if(nnode .eq. 4) then
        st=s*t
! shape functions for 4 nodes element
        shape(1)=(1-t-s+st)*0.25
        shape(2)=(1-t+s-st)*0.25
        shape(3)=(1+t+s+st)*0.25
        shape(4)=(1+t-s-st)*0.25

! shape function derivates
        deriv(1,1)=(-1+t)*0.25
        deriv(1,2)=(+1-t)*0.25
        deriv(1,3)=(+1+t)*0.25
        deriv(1,4)=(-1-t)*0.25
        deriv(2,1)=(-1+s)*0.25
        deriv(2,2)=(-1-s)*0.25
        deriv(2,3)=(+1+s)*0.25
        deriv(2,4)=(+1-s)*0.25
        return
    endif

    if(nnode .eq. 8) then
        s2=s*2.0
        t2=t*2.0
        ss=s*s
        tt=t*t
        st=s*t
        sst=s*s*t
        stt=s*t*t
        st2=s*t*2.0

! shape functions for 8 noded element
        shape(1)=(-1.0+ST+SS+TT-SST-STT)/4.0
        shape(2)=(1.0-T-SS+SST)/2.0
        shape(3)=(-1.0-ST+SS+TT-SST+STT)/4.0
        shape(4)=(1.0+S-TT-STT)/2.0
        shape(5)=(-1.0+ST+SS+TT+SST+STT)/4.0
        shape(6)=(1.0+T-SS-SST)/2.0
        shape(7) =(-1.0-ST+SS+TT+SST-STT)/4.0
        shape(8)=(1.0-S-TT+STT)/2.0

! shape function derivAtives
        deriv(1,1)=(T+S2-ST2-TT)/4.0
        deriv(1,2)=-S+ST
        deriv(1,3)=(-T+S2-ST2+TT)/4.0
        deriv(1,4)=(1.0-TT)/2.0
        deriv(1,5)=(T+S2+ST2+TT)/4.0
        deriv(1,6)=-S-ST
        deriv(1,7)=(-T+S2+ST2-TT)/4.0
        deriv(1,8)=(-1.0+TT)/2.0
        deriv(2,1)=(S+T2-SS-ST2)/4.0
        deriv(2,2)=(-1.0+SS)/2.0
        deriv(2,3)=(-S+T2-SS+ST2)/4.0
        deriv(2,4)=-T-ST
        deriv(2,5)=(S+T2+SS+ST2)/4.0
        deriv(2,6)=(1.0-SS)/2.0
        deriv(2,7)=(-S+T2*SS-ST2)/4.0
        deriv(2,8)=-T+ST
        return
    endif

    if(nnode .eq. 9) then
        ss=s*s
        st=s*t
        tt=t*t
        s1=s+1.0
        t1=t+1.0
        s2=s*2.0
        t2=t*2.0
        s9=s-1.0
        t9=t-1.0

! shape functions for 9 noded element
        shape(1)=0.25*S9*ST*T9
        shape(2)=0.5*(1.0-SS)*T*T9
        shape(3)=0.25*S1*ST*T9
        shape(4)=0.5*S*S1*(1.0-TT)
        shape(5)=0.25*S1*ST*T1
        shape(6)=0.5*(1.0-SS)*T*T1
        shape(7)=0.25*S9*ST*T1
        shape(8)=0.5*S*S9*(1.0-TT)
        shape(9)=(1.0-SS)*(1.0-TT)

! shape function derivatives
        deriv(1,1)=0.25*T*T9*(-1.0+S2)
        deriv(1,2)=-ST*T9
        deriv(1,3)=0.25*(1.0+S2)*T*T9
        deriv(1,4)=0.5*(1.0+S2)*(1.0-TT)
        deriv(1,5)=0.25*(1.0+S2)*T*T1
        deriv(1,6)=-ST*T1
        deriv(1,7)=0.25*(-1.0+S2) *T*T1
        deriv(1,8)=0.5*(-1.0+s2)*(1.0-TT)
        deriv(1,9)=-S2*(1.0-TT)
        deriv(2,1)=0.25*(-1.0+T2)*S*S9
        deriv(2,2)=0.5*(1.0-SS)*(-1.0+T2)
        deriv(2,3)=0.25*S*S1*(-1.0+T2)
        deriv(2,4)=-ST*S1
        deriv(2,5)=0.25*S*51*(1.0+T2)
        deriv(2,6)=0.5*(1.0-SS)*(1.0+T2)
        deriv(2,7)=0.25*S*S9*(1.0+T2)
        deriv(2,8)=-ST*S9
        deriv(2,9)=-T2*(1.0-SS)

        return
    endif

end subroutine
