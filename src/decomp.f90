subroutine decomp(stiff ,maxai ,neqns ,ishot )
!
! Factorises (L)#(D)#(L) transpose of stiffness matrix
!
    implicit none

    integer :: ieqns, imaxa, neqns, ishot, lower, kuper, khigh, ksize, icoun, juper, jhigh, kmaxa, ndiag, ncolm, icolm, jmaxa
    integer :: maxai(1)

    real :: count, bsumm, ratio
    real :: stiff(1)

    if(neqns .eq. 1) return

    do 200 ieqns=1,neqns
        imaxa = maxai(ieqns)
        lower = imaxa + 1
        kuper = maxai(ieqns+1) - 1
        khigh = kuper - lower
        if(khigh) 304,240,210
210     ksize = ieqns-khigh
        icoun = 0
        juper=kuper
        do jhigh=1,khigh
            icoun = icoun+1
            juper = juper-1
            kmaxa = maxai(ksize)
            ndiag = maxai(ksize+1)-kmaxa - 1
            if(ndiag .le. 0) then
                ksize=ksize+1
            else
                ncolm=min0(icoun, ndiag)
                count=0.0
                do icolm = 1, ncolm
                    count=count+stiff (kmaxa+icolm) *stiff(juper+icolm)
                end do
                stiff(juper) = stiff(juper) -count
                ksize=ksize+1
            endif
        end do
240     ksize=ieqns
        bsumm=0.0
        do icolm = lower, kuper
            ksize = ksize - 1
            jmaxa = maxai(ksize)
            ratio = stiff(icolm)/stiff(jmaxa)
            bsumm = bsumm + ratio*stiff(icolm)
            stiff(icolm)=ratio
        end do
        stiff(imaxa) = stiff(imaxa) - bsumm

304     if(stiff(imaxa)) 310,310,200

310     if(ishot.EQ.0) GO TO 320

        if(stiff(imaxa) .EQ.0) stiff( imaxa) =-1.E-16
        GO TO 200

320     WRITE(6,2000) ieqns,stiff( imaxa)
        STOP

200 CONTINUE
    RETURN

2000 FORMAT(//'STOP - stiffNESS MATRIX NOT POSITIVE DEFINITE ', &
        // 'NONPOSITIVE PIVOT FOR EQUATION ',I4,//' PIVOT = ',E20.12 )

end subroutine
