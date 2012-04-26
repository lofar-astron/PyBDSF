C$TEST PDEW
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PDEW
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM POSTU
C
C***********************************************************************
C  THE PORT STACK
C
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(4000)
      REAL WS(1000)
      EQUIVALENCE (DS(1),WS(1))
C
C  TIME FOR THE FUNCTION UOFX.
C
      COMMON /TIME/ T
      REAL T
C
C  MAPPING PARAMETERS FOR UOFX.
C
      COMMON /PARAM/ VC, X
      REAL VC(4), X(3)
      EXTERNAL DEE, HANDLE, UOFX, BC, AF
      INTEGER NDX, K, IMMM, ISTKGT
      INTEGER NU, NV, IMESH, ILUMB, NMESH
      REAL ERRPAR(2), TSTART, TSTOP, V(4), DT, XB(3), U(1000)
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(4000, 4)
      CALL ENTER(1)
      NU = 1
      NV = 4
      ERRPAR(1) = 0
C ABSOLUTE ERROR.
      ERRPAR(2) = 1E-2
      TSTART = 0
      TSTOP = 3.14
      DT = 0.4
      K = 4
C NDX UNIFORM MESH POINTS ON EACH INTERVAL OF XB.
      NDX = 6
      XB(1) = 0
      XB(2) = 1
      XB(3) = 2
C GET MESH ON PORT STACK.
      IMESH = ILUMB(XB, 3, NDX, K, NMESH)
C MAKE 1 OF MULTIPLICITY K-1.
      IMESH = IMMM(IMESH, NMESH, 1E0, K-1)
      X(1) = -3.14
      X(2) = 3.14/2.
      X(3) = 3.14
C INITIAL VALUES FOR V.
      CALL LPLMG(3, X, VC)
C GET U ON THE PORT STACK.
      IU = ISTKGT(NMESH-K, 3)
C UOFX NEEDS TIME.
      T = TSTART
C THE INITIAL HEIGHT OF THE JUMP.
      VC(4) = 1
C UOFX NEEDS V FOR MAPPING.
      CALL MOVEFR(NV, VC, V)
C INITIAL CONDITIONS FOR U.
      CALL L2SFF(UOFX, K, WS(IMESH), NMESH, U)
C OUTPUT ICS.
      CALL HANDLE(T-1., U, V, T, U, V, NU, NMESH-K, NV, K, WS(
     1   IMESH), NMESH, DT, TSTOP)
      CALL POST(U, NU, K, WS(IMESH), NMESH, V, NV, TSTART, TSTOP,
     1   DT, AF, BC, DEE, ERRPAR, HANDLE)
      CALL LEAVE
      CALL WRAPUP
      STOP
      END
      SUBROUTINE AF(T, XI, NX, U, UX, UT, UTX, NU, V, VT, NV,
     *              A, AU, AUX, AUT, AUTX, AV, AVT,
     *              F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NV, NX
      REAL T, XI(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), UTX(NX, NU)
      REAL V(NV), VT(NV),
     *     A(NX,NU),AU(NX,NU,NU),AUX(NX,NU,NU),AUT(NX,NU,NU),
     *     AUTX(NX,NU,NU),AV(NX,NU,NV),AVT(NX,NU,NV),
     *     F(NX,NU),FU(NX,NU,NU),FUX(NX,NU,NU),FUT(NX,NU,NU),
     *     FUTX(NX,NU,NU),FV(NX,NU,NV),FVT(NX,NU,NV)
      COMMON /POSTF/ FAILED
      LOGICAL FAILED
      INTEGER I
      REAL COS, SIN, XXI(99), XTV(99), XVV(99), X(99)
      REAL XXIV(99), AX(99), FX(99), XT(99), XV(99)
      LOGICAL TEMP
      TEMP = V(2) .LE. V(1)
      IF (.NOT. TEMP) TEMP = V(2) .GE. V(3)
      IF (.NOT. TEMP) GOTO 1
         FAILED = .TRUE.
         RETURN
C MAP XI INTO X.
   1  CALL LPLM(XI, NX, V, 3, X, XXI, XXIV, XV, XVV, XT, XTV)
C MAP U INTO X SYSTEM.
      CALL POSTU(XI, X, XT, XXI, XV, VT, NX, 3, UX, UT, NU, AX, FX)
      DO  4 I = 1, NX
         A(I, 1) = -U(I, 1)
         AU(I, 1, 1) = -1
         F(I, 1) = UT(I, 1)
         FUT(I, 1, 1) = 1
         IF (XI(I) .GT. 1.) GOTO 2
            F(I, 1) = F(I, 1)-2.*COS(X(I)+T)
            FX(I) = 2.*SIN(X(I)+T)
            GOTO  3
   2        F(I, 1) = F(I, 1)-VT(4)
            FVT(I, 1, 4) = -1
            F(I, 1) = F(I, 1)+2.*SIN(X(I)+T)
            FX(I) = 2.*COS(X(I)+T)
   3     CONTINUE
   4     CONTINUE
C MAP A AND F INTO XI SYSTEM.
      CALL POSTI(XI, X, XT, XXI, XV, XTV, XXIV, XVV, NX, UX, UT, NU, V
     1   , VT, NV, 1, 3, A, AX, AU, AUX, AUT, AUTX, AV, AVT, F, FX, FU
     2   , FUX, FUT, FUTX, FV, FVT)
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV,
     *              B, BU, BUX, BUT, BUTX, BV, BVT)
      INTEGER NU, NV
      REAL T,L,R,U(NU,2),UX(NU,2),UT(NU,2),UTX(NU,2),V(NV),VT(NV)
      REAL B(NU,2),BU(NU,NU,2),BUX(NU,NU,2),BUT(NU,NU,2),BUTX(NU,NU,2),
     *     BV(NU,NV,2),BVT(NU,NV,2)
      B(1, 1) = U(1, 1)-SIN(T-3.14)
C U(-PI,T) = SIN(-PI+T).
      BU(1, 1, 1) = 1
      RETURN
      END
      SUBROUTINE DEE(T, K, X, NX, U, UT, NU, NXMK, V, VT, NV,
     *               D, DU, DUT, DV, DVT)
      INTEGER NXMK, NU, NV, NX, K
      REAL T, X(NX), U(NXMK, NU), UT(NXMK, NU), V(NV), VT(NV)
      REAL D(NV),DU(NV,NXMK,NU),DUT(NV,NXMK,NU),DV(NV,NV),DVT(NV,NV)
      INTEGER INTRVR, I, ILEFT
      REAL BX(10), XX(1), R1MACH
      INTEGER TEMP
      D(1) = V(1)+3.14
C X(0,V) = -PI.
      DV(1, 1) = 1
C XX(1) = 1 + A ROUNDING ERROR.
      XX(1) = R1MACH(4)+1.
      ILEFT = INTRVR(NX, X, XX(1))
C GET THE B-SPLINE BASIS AT XX.
      CALL BSPLN(K, X, NX, XX, 1, ILEFT, BX)
      D(2) = -V(4)
C U(X(T)+,T) - JUMP = 0.
      DV(2, 4) = -1
      DO  1 I = 1, K
         TEMP = ILEFT+I-K
         D(2) = D(2)+U(TEMP, 1)*BX(I)
         TEMP = ILEFT+I-K
         DU(2, TEMP, 1) = BX(I)
   1     CONTINUE
      D(3) = V(3)-3.14
C X(2,V) = +PI.
      DV(3, 3) = 1
C JUMP + D( X(1,V(T)) )/DT = 0.
      D(4) = VT(2)+V(4)
      DVT(4, 2) = 1
      DV(4, 4) = 1
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, V0, T, U, V, NU, NXMK, NV,
     *                  K, X, NX, DT, TSTOP)
      INTEGER NXMK, NU, NV, NX, K
      REAL T0, U0(NXMK, NU), V0(NV), T, U(NXMK, NU), V(NV),
     *     X(NX), DT, TSTOP
      COMMON /PARAM/ VC, XX
      REAL VC(4), XX(3)
      COMMON /TIME/ TT
      REAL TT
      EXTERNAL UOFX
      INTEGER I1MACH
      REAL EU, EESFF, EV(2)
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .NE. T) GOTO 2
         TEMP = I1MACH(2)
         WRITE (TEMP,  1) T, DT
   1     FORMAT (16H RESTART FOR T =, 1PE10.2, 7H   DT =, 1PE10.2)
         RETURN
   2  TT = T
C UOFX NEEDS V FOR MAPPING.
      CALL MOVEFR(NV, V, VC)
      EU = EESFF(K, X, NX, U, UOFX)
C ERROR IN POSITION OF SHOCK.
      EV(1) = V(2)-(3.14/2.-T)
C ERROR IN HEIGHT OF SHOCK.
      EV(2) = V(4)-1.
      TEMP = I1MACH(2)
      WRITE (TEMP,  3) T, EU, EV
   3  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 1PE10.2, 6H   V =, 2(
     1   1PE10.2))
      RETURN
      END
      SUBROUTINE UOFX(XI, NX, U, W)
      INTEGER NX
      REAL XI(NX), U(NX), W(NX)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(4000)
      COMMON /PARAM/ VC, X
      REAL VC(4), X(3)
      COMMON /TIME/ T
      REAL T
      INTEGER IXV, IXX, ISTKGT, I, IS(1000)
      REAL EWE, RS(1000), WS(1000)
      LOGICAL LS(1000)
      INTEGER TEMP
      EQUIVALENCE (DS(1), WS(1), RS(1), IS(1), LS(1))
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
      IXX = ISTKGT(NX, 3)
C SPACE FOR X AND XV.
      IXV = ISTKGT(3*NX, 3)
C MAP INTO USER SYSTEM.
      CALL LPLMX(XI, NX, VC, 3, WS(IXX), WS(IXV))
      DO  1 I = 1, NX
         TEMP = IXX+I
         U(I) = EWE(T, WS(TEMP-1), VC(2))
         IF (XI(I) .GT. 1.) U(I) = U(I)+1.
   1     CONTINUE
      CALL LEAVE
      RETURN
      END
      REAL FUNCTION EWE(T, X, XBREAK)
      REAL T, X, XBREAK
      REAL COS, SIN
      IF (X .GE. XBREAK) GOTO 1
         EWE = SIN(X+T)
         RETURN
   1     IF (X .LE. XBREAK) GOTO 2
            EWE = COS(X+T)
            RETURN
C/6S
   2        CALL SETERR(17HEWE - X == XBREAK, 17, 1, 2)
C/7S
C  2        CALL SETERR('EWE - X == XBREAK', 17, 1, 2)
C/
   3  CONTINUE
   4  STOP
      END
