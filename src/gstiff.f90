subroutine gstiff(coord, epstn, intgr, istep, kstep, leqns, lnods , matno, maxai, maxaj, ncrit , ndime, ndofn , nelem , ngaus,     &
                  nlaps, nmats, nnode, npoin ,NSTRE ,NTYPE ,nwmtl ,nwktl ,POSGP, PROPS ,STIFF ,STIFI ,STRSG ,TDISP ,WEIGP)
! **********************************************************************
!
! EVALUATES GEOMETRICALLY NONLINEAR STIFFNESS MATRIX
! FOR 2-D PLANE STRESS/STRAIN 2-D ELEMENT
!
! **********************************************************************
   IMPLICIT NONE

   INCLUDE 'param.inc'

   integer :: istep, kstep, leqns, matno, maxaj, ncrit, ndime, ndofn, nelem, ngaus, nlaps, nmats, nnode,                           &
      npoin, NSTRE, NTYPE, nwmtl, nwktl, KOUNT, LSTEP, KGAUS, NEVAB, IWKTL, IELEM, IPOSN, INODE, IDIME, LNODE, NPOSN,              &
      IEVAB, JEVAB, KGASP, NSTR1, LPROP, ISIZE, IGAUS, JGAUS, ISTR1, ISTRE, JSTRE

   integer :: lnods(nelem, 1), maxai(1), intgr(1) 

   real :: coord(npoin,1), epstn(1), props, STIFI, STRSG, TDISP, WEIGP, AVECT, CARTD, DVECT, DBMAT,                  &
      DEVIA, DJACM, DERIV, SHAPE, STIFP, TWOPI, DISPT, YOUNG, THICK, HARDS, FRICT, exisp,                                          &
      etasp, DJACB, TELEM, DVOLU, SINT3, STEFF, THETA, VARJ2, YIELD, ABETA, POISS

   real :: dlcod(2,9), elcod(2, 9), estif(171), stres(4), dmatx(4,4), gpcod(2,9), bmatx(4, 18), stiff(1), posgp(1)

!   DIMENSION ELCOD(2,9), AVECT(4), BMATX(4,18), CARTD(2,9), DVECT(4),                &
!   +             PROPS(MMATS, MPROP) ,DBMAT(4,18) ,GPCOD(2,9) ,DEVIA(4), leqns(18, 1), STRSG(4, 1), DLCOD(2,9), STRES(4),         &
!   +             ESTIF(171), DJACM(2, 2), DERIV(2,9), SHAPE(9)

   DIMENSION maxaj(1), TDISP(1), STIFI(1), STIFP(1), WEIGP(1), matno(1)

   IF(istep .EQ. 1) GO TO 200
   KOUNT=(LSTEP/kstep)*kstep
   IF(KOUNT .NE. istep) RETURN
200 CONTINUE
   TWOPI=6.283185307179586
   KGAUS=0

   !     LOOP OVER EACH ELEMENT
   NSTR1=0
   NEVAB=ndofn*nnode
   DO 500 IWKTL=1, nwktl
      STIFP(IWKTL) = 0.0
      STIFI(IWKTL) = 0.0
500 CONTINUE

   DO 70 IELEM=1,nelem
      LPROP=matno(IELEM)

      !       EVALUATE THE coordINATES OF THE ELEMENT NODAL POINTS
      IPOSN=0
      DO 10 INODE=1,nnode
         LNODE=lnods( IELEM, INODE)
         DO 15 IDIME=1,ndime
            IPOSN=IPOSN+1
            NPOSN=leqns( IPOSN, IELEM)
            IF(NPOSN .EQ. 0) DISPT=0.
            IF(NPOSN .NE. 0) DISPT=TDISP(NPOSN)
            DLCOD(IDIME, INODE) = coord(LNODE, IDIME)+DISPT
            ELCOD(IDIME, INODE) = coord(LNODE, IDIME)
15       CONTINUE
10    CONTINUE
      YOUNG=PROPS(LPROP, 1)
      POISS=PROPS(LPROP, 2)
      THICK=PROPS(LPROP, 3)
      HARDS=PROPS(LPROP, 7)
      FRICT=PROPS(LPROP, 8)

      !       INITIALIZE THE ELEMENT STIFFNESS MATRIX 171=NEVAB*(NEVAB+1)/2
      DO 20 ISIZE=1,171
         ESTIF(ISIZE)=0.0
20    CONTINUE
      KGASP=0

      !       ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
        DO 50 IGAUS=1,ngaus
         exisp=POSGP( IGAUS)
         DO 55 JGAUS=1,ngaus
            etasp=POSGP( JGAUS)
            KGASP=KGASP+1
            KGAUS=KGAUS+1
            CALL MODPS(dmatx,LPROP,nmats, NSTRE,NTYPE, PROPS)
            CALL SFR2(DERIV,nnode, SHAPE, exisp, etasp)
            CALL JACOB2(CARTD, DERIV, DJACB, ELCOD, GPCOD, TELEM, KGASP, nnode, SHAPE)
            CALL JACOBD(CARTD, DLCOD, DJACM, ndime, nlaps, nnode)
            DVOLU=DJACB*WEIGP( IGAUS)*WEIGP(JGAUS)
            IF(NTYPE .EQ. 3) DVOLU=DVOLU*TWOPI*GPCOD(1,KGASP)
            IF(NTYPE .EQ. 1) DVOLU=DVOLU*THICK

            ! EVALUATE THE B AND DB MATRICES
            CALL BLARGE(BMATX, CARTD, DJACM, DLCOD, GPCOD, KGASP, nlaps, nnode, NTYPE, SHAPE)
            IF(nlaps .EQ. 2 .OR. nlaps .EQ. 0) GO TO 80
            IF(istep .EQ. 1) GO TO 80
            IF(epstn(KGAUS) .EQ.0.0) GO TO 80
            DO 90 ISTR1=1,NSTR1
               STRES(ISTR1)=STRSG(ISTR1,KGAUS)
90          CONTINUE
            CALL INVAR(DEVIA, LPROP, ncrit, nmats, PROPS, SINT3, STEFF, STRES, THETA, VARJ2, YIELD)
            CALL YIELDF(AVECT, DEVIA, FRICT, ncrit, SINT3, STEFF, THETA, VARJ2)
            CALL FLOWPL(AVECT, ABETA, DVECT, HARDS, NTYPE, POISS, YOUNG)
            DO 100 ISTRE=1,NSTRE
               DO 105 JSTRE=1,NSTRE
                  dmatx(ISTRE, JSTRE)=dmatx(ISTRE,JSTRE)- ABETA*DVECT(ISTRE)*DVECT(JSTRE)
105            CONTINUE
100         CONTINUE
80          CONTINUE
            CALL DINTOB(BMATX,DBMAT,dmatx, NEVAB, NSTRE)

            !           EVALUATE GEOMETRIC STIFFNESS TERMS
            IF(nlaps.LT.2) GO TO 85
            CALL GEOMST(CARTD, DVOLU, ESTIF, KGAUS, ndofn, nnode, STRSG, SHAPE, NTYPE, GPCOD, KGASP)

            !           CALCULATE THE ELEMENT STIFFNESSES
85          KOUNT=0
            DO 30 IEVAB=1,NEVAB
               DO 35 JEVAB=IEVAB, NEVAB
                  KOUNT=KOUNT+1
                  DO 40 ISTRE=1,NSTRE
                     ESTIF(KOUNT)=ESTIF(KOUNT)+BMATX(ISTRE,IEVAB)*DBMAT(ISTRE,JEVAB)*DVOLU
40                CONTINUE
35             CONTINUE
30          CONTINUE
55       CONTINUE
50    CONTINUE

      ! GENERATES GLOBAL STIFFNSS MATRIX IN COMPACTED COLUMN FORM
      IF(intgr(IELEM) .EQ. 2) GO TO 210
      CALL ADDBAN(STIFI, maxai, ESTIF, leqns(1, IELEM), NEVAB)
210   CALL ADDBAN(STIFF, maxaj, ESTIF, leqns(1, IELEM), NEVAB)
70 CONTINUE
   WRITE(6,900) (STIFI(IWKTL) ,IWKTL=1,nwmtl)
900 FORMAT(10E12.4)

   RETURN
END

