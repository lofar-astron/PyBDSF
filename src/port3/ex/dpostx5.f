C$TEST DPT5
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DPT5
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM DPOST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(4000)
      COMMON /TIME/ T
      DOUBLE PRECISION T
      COMMON /PARAM/ VC, X
      DOUBLE PRECISION VC(3), X(3)
      EXTERNAL DEE, HANDLE, UOFX, BC, AF
      INTEGER NDX, IDLUMB, ISTKGT, K, IU, IS(1000)
      INTEGER NU, NV, IMMMD, IMESH, NMESH
      REAL ERRPAR(2), RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION TSTART, V(3), DT, XB(3), WS(500), TSTOP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST DPOST ON
C      U SUB T = ( K(T,X) * U SUB X ) SUB X + G      ON (-1,+2) * (0,+1)
C WITH A MOVING FRONT X(T) CHARACTERIZED BY U(X(T),T) == 1 AND
C    JUMP ACROSS X(T) OF K(T,X) U SUB X = - 3 * X'(T).
C WHERE K(T,X) IS PIECEWISE CONSTANT, SAY
C            1 FOR X < X(T)
C   K(T,X) =
C            2 FOR X > X(T)
C AND G IS CHOSEN SO THAT THE SOLUTION IS
C               EXP(X-X(T))  FOR X < X(T)
C      U(X,T) = 
C               EXP(X(T)-X)  FOR X > X(T)
C AND X(1,T) = T. THE MOVING FRONT IS TRACKED
C IMPLICITLY BY FORCING U(X(1,T),T) = 1 AS A PSEUDO-RANKINE-HEUGONIOT RE
CLATION.
C V(1,2,3) GIVES THE MOVING MESH.
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(4000, 4)
      CALL ENTER(1)
      NU = 1
      NV = 3
      ERRPAR(1) = 0
C ABSOLUTE ERROR.
      ERRPAR(2) = 1E-2
      TSTART = 0
      TSTOP = 1
      DT = 0.1
      K = 4
C NDX UNIFORM MESH POINTS ON EACH INTERVAL OF XB ARRAY.
      NDX = 6
      XB(1) = 0
      XB(2) = 1
      XB(3) = 2
C GET MESH ON PORT STACK.
      IMESH = IDLUMB(XB, 3, NDX, K, NMESH)
C MAKE 1 OF MULTIPLICITY K-1.
      IMESH = IMMMD(IMESH, NMESH, 1D0, K-1)
      X(1) = -1
      X(2) = 0
      X(3) = 2
C INITIAL VALUES FOR V.
      CALL DLPLMG(3, X, VC)
C GET U ON THE PORT STACK.
      IU = ISTKGT(NMESH-K, 4)
C UOFX NEEDS TIME.
      T = TSTART
C UOFX NEEDS V FOR MAPPING.
      CALL MOVEFD(NV, VC, V)
C INITIAL CONDITIONS FOR U.
      CALL DL2SFF(UOFX, K, WS(IMESH), NMESH, WS(IU))
C OUTPUT THE ICS.
      CALL HANDLE(T-1D0, WS(IU), V, T, WS(IU), V, NU, NMESH-K, NV, K, 
     1   WS(IMESH), NMESH, DT, TSTOP)
      CALL DPOST(WS(IU), NU, K, WS(IMESH), NMESH, V, NV, TSTART, TSTOP
     1   , DT, AF, BC, DEE, ERRPAR, HANDLE)
      CALL LEAVE
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, XI, NX, U, UX, UT, UTX, NU, V, VT, NV, A, 
     1   AU, AUX, AUT, AUTX, AV, AVT, F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NV, NX
      DOUBLE PRECISION T, XI(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), 
     1   UTX(NX, NU)
      DOUBLE PRECISION V(NV), VT(NV), A(NX, NU), AU(NX, NU, NU), AUX(NX,
     1   NU, NU), AUT(NX, NU, NU)
      DOUBLE PRECISION AUTX(NX, NU, NU), AV(NX, NU, NV), AVT(NX, NU, NV)
     1   , F(NX, NU), FU(NX, NU, NU), FUX(NX, NU, NU)
      DOUBLE PRECISION FUT(NX, NU, NU), FUTX(NX, NU, NU), FV(NX, NU, NV)
     1   , FVT(NX, NU, NV)
      COMMON /DPOSTF/ FAILED
      LOGICAL FAILED
      INTEGER I
      DOUBLE PRECISION KAY, XXI(99), XTV(99), XVV(99), X(99), DEXP
      DOUBLE PRECISION XXIV(99), AX(99), FX(99), XT(99), XV(99)
      LOGICAL TEMP
      TEMP = V(2) .LE. V(1)
      IF (.NOT. TEMP) TEMP = V(2) .GE. V(3)
      IF (.NOT. TEMP) GOTO 1
         FAILED = .TRUE.
         RETURN
C MAP XI INTO X.
   1  CALL DLPLM(XI, NX, V, 3, X, XXI, XXIV, XV, XVV, XT, XTV)
C MAP U INTO X SYSTEM.
      CALL DPOSTU(XI, X, XT, XXI, XV, VT, NX, 3, UX, UT, NU, AX, FX)
      DO  7 I = 1, NX
         IF (XI(I) .GT. 1D0) GOTO 2
            KAY = 1
            GOTO  3
   2        KAY = 2
   3     A(I, 1) = KAY*UX(I, 1)
         AUX(I, 1, 1) = KAY
         IF (XI(I) .GT. 1D0) GOTO 4
            A(I, 1) = A(I, 1)-3D0*VT(2)
            AVT(I, 1, 2) = -3
   4     F(I, 1) = UT(I, 1)
         FUT(I, 1, 1) = 1
         IF (XI(I) .GT. 1D0) GOTO 5
            F(I, 1) = F(I, 1)+2D0*DEXP(X(I)-T)
            FX(I) = 2D0*DEXP(X(I)-T)
            GOTO  6
   5        F(I, 1) = F(I, 1)+DEXP(T-X(I))
            FX(I) = -DEXP(T-X(I))
   6     CONTINUE
   7     CONTINUE
C MAP A AND F INTO XI SYSTEM.
      CALL DPOSTI(XI, X, XT, XXI, XV, XTV, XXIV, XVV, NX, UX, UT, NU, V,
     1   VT, NV, 1, 3, A, AX, AU, AUX, AUT, AUTX, AV, AVT, F, FX, FU, 
     2   FUX, FUT, FUTX, FV, FVT)
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV, B, BU,
     1   BUX, BUT, BUTX, BV, BVT)
      INTEGER NU, NV
      DOUBLE PRECISION T, L, R, U(NU, 2), UX(NU, 2), UT(NU, 2)
      DOUBLE PRECISION UTX(NU, 2), V(NV), VT(NV), B(NU, 2), BU(NU, NU, 2
     1   ), BUX(NU, NU, 2)
      DOUBLE PRECISION BUT(NU, NU, 2), BUTX(NU, NU, 2), BV(NU, NV, 2), 
     1   BVT(NU, NV, 2)
      DOUBLE PRECISION DEXP
      B(1, 1) = U(1, 1)-DEXP((-1D0)-T)
      B(1, 2) = U(1, 2)-DEXP(T-2D0)
      BU(1, 1, 1) = 1
      BU(1, 1, 2) = 1
      RETURN
      END
      SUBROUTINE DEE(T, K, X, NX, U, UT, NU, NXMK, V, VT, NV, D, 
     1   DU, DUT, DV, DVT)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      DOUBLE PRECISION T, X(NX), U(NXMK, NU), UT(NXMK, NU), V(NV), VT(
     1   NV)
      DOUBLE PRECISION D(NV), DU(NV, NXMK, NU), DUT(NV, NXMK, NU), DV(
     1   NV, NV), DVT(NV, NV)
      INTEGER INTRVD, I, ILEFT
      DOUBLE PRECISION BX(10), XX(1)
      INTEGER TEMP
      D(1) = V(1)+1D0
C X(0,V) = -1.
      DV(1, 1) = 1
      XX(1) = 1
C FIND 1 IN THE MESH.
      ILEFT = INTRVD(NX, X, XX(1))
C GET THE B-SPLINE BASIS AT XX.
      CALL DBSPLN(K, X, NX, XX, 1, ILEFT, BX)
C U(X(1,V),T) = 1.
      D(2) = -1
      DO  1 I = 1, K
         TEMP = ILEFT+I-K
         D(2) = D(2)+U(TEMP, 1)*BX(I)
         TEMP = ILEFT+I-K
         DU(2, TEMP, 1) = BX(I)
   1     CONTINUE
      D(3) = V(3)-2D0
C X(2,V) = +2.
      DV(3, 3) = 1
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, V0, T, U, V, NU, NXMK, NV, K, X, 
     1   NX, DT, TSTOP)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      DOUBLE PRECISION T0, U0(NXMK, NU), V0(NV), T, U(NXMK, NU), V(NV)
      DOUBLE PRECISION X(NX), DT, TSTOP
      COMMON /PARAM/ VC, XX
      DOUBLE PRECISION VC(3), XX(3)
      COMMON /TIME/ TT
      DOUBLE PRECISION TT
      EXTERNAL UOFX
      INTEGER I1MACH
      DOUBLE PRECISION DEESFF, EU, EV(3)
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .NE. T) GOTO 2
         TEMP = I1MACH(2)
         WRITE (TEMP,  1) T
   1     FORMAT (16H RESTART FOR T =, 1PE10.2)
         RETURN
   2  TT = T
C UOFX NEEDS V FOR MAPPING.
      CALL MOVEFD(NV, V, VC)
      EU = DEESFF(K, X, NX, U, UOFX)
      EV(1) = V(1)+1D0
      EV(2) = V(2)-T
      EV(3) = V(3)-2D0
      TEMP = I1MACH(2)
      WRITE (TEMP,  3) T, EU, EV
   3  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 1PE10.2, 6H   V =, 3(
     1   1PE10.2))
      RETURN
      END
      SUBROUTINE UOFX(XI, NX, U, W)
      INTEGER NX
      DOUBLE PRECISION XI(NX), U(NX), W(NX)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PARAM/ VC, X
      DOUBLE PRECISION VC(3), X(3)
      COMMON /TIME/ T
      DOUBLE PRECISION T
      INTEGER IXV, IXX, ISTKGT, I, IS(1000)
      REAL RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION DEXP, WS(500), XOFXI
      INTEGER TEMP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
      IXX = ISTKGT(NX, 4)
C SPACE FOR X AND XV.
      IXV = ISTKGT(3*NX, 4)
C MAP INTO USER SYSTEM.
      CALL DLPLMX(XI, NX, VC, 3, WS(IXX), WS(IXV))
      DO  3 I = 1, NX
         TEMP = IXX+I
         XOFXI = WS(TEMP-1)
         IF (XI(I) .GT. 1D0) GOTO 1
            U(I) = DEXP(XOFXI-T)
            GOTO  2
   1        U(I) = DEXP(T-XOFXI)
   2     CONTINUE
   3     CONTINUE
      CALL LEAVE
      RETURN
      END
