C$TEST PST7
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PST7
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM POST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(4000)
      COMMON /TIME/ T
      REAL T
      COMMON /PARAM/ VC, X, XI0
      REAL VC(4), X(3), XI0
      EXTERNAL DEE, HANDLE, UOFX, BC, AF
      INTEGER NDX, ISTKGT, K, IMMM, IU, IS(1000)
      INTEGER NU, NV, IMESH, ILUMB, NMESH
      REAL ERRPAR(2), TSTART, D, V(4), DT, XB(3)
      REAL RS(1000), WS(1000), TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST  POST ON
C      U SUB T = U SUB XX + F      ON (20,10**6)
C WHERE F IS CHOSEN SO THAT THE SOLUTION IS
C      U(X,T) = EXP(-X*T),
C AND X(1,T) IS CHOSEN SO THAT THE BOUNDARY-LAYER IS TRACKED
C IMPLICITLY BY FORCING U(X(1,T)/2.3/D,T) = 1/E.
C THIS IS THE SAME AS REQUIRING THE EXACT SOLUTION TO HAVE
C U(X(1,T),T) = 10 ** -D.
C V(1,2,3) GIVES THE MOVING MESH, V(4) IS TIME.
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(4000, 4)
      CALL ENTER(1)
      NU = 1
      NV = 4
      ERRPAR(1) = 1E-2
C MIXED RELATIVE AND ABSOLUTE ERROR.
      ERRPAR(2) = 1E-2
      D = 3
C W(XI0,T) = 1/E.
      XI0 = 1./2.3/D
      TSTART = 20
      TSTOP = 1E+6
      DT = 1E-2
      K = 4
C NDX UNIFORM MESH POINTS ON EACH INTERVAL OF XB.
      NDX = 6
      XB(1) = 0
      XB(2) = 1
      XB(3) = 2
C GET MESH ON PORT STACK.
      IMESH = ILUMB(XB, 3, NDX, K, NMESH)
C MAKE 1D0 OF MULTIPLICITY K-1.
      IMESH = IMMM(IMESH, NMESH, 1E0, K-1)
      X(1) = 0
      X(2) = 2.3*D/TSTART
      X(3) = 1
C INITIAL VALUES FOR V.
      CALL LPLMG(3, X, VC)
C GET U ON PORT STACK.
      IU = ISTKGT(NMESH-K, 3)
C UOFX NEEDS TIME.
      T = TSTART
      VC(4) = TSTART
C UOFX NEEDS V FOR MAPPING.
      CALL MOVEFR(NV, VC, V)
C INITIAL CONDITIONS FOR U.
      CALL L2SFF(UOFX, K, WS(IMESH), NMESH, WS(IU))
C OUTPUT ICS.
      CALL HANDLE(T-1., WS(IU), V, T, WS(IU), V, NU, NMESH-K, NV, K, WS(
     1   IMESH), NMESH, DT, TSTOP)
      CALL POST(WS(IU), NU, K, WS(IMESH), NMESH, V, NV, TSTART, TSTOP, 
     1   DT, AF, BC, DEE, ERRPAR, HANDLE)
      CALL LEAVE
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, XI, NX, U, UX, UT, UTX, NU, V, VT, NV, A, 
     1   AU, AUX, AUT, AUTX, AV, AVT, F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NV, NX
      REAL T, XI(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), UTX(NX, NU)
      REAL V(NV), VT(NV), A(NX, NU), AU(NX, NU, NU), AUX(NX, NU, NU), 
     1   AUT(NX, NU, NU)
      REAL AUTX(NX, NU, NU), AV(NX, NU, NV), AVT(NX, NU, NV), F(NX, NU),
     1   FU(NX, NU, NU), FUX(NX, NU, NU)
      REAL FUT(NX, NU, NU), FUTX(NX, NU, NU), FV(NX, NU, NV), FVT(NX, 
     1   NU, NV)
      COMMON /POSTF/ FAILED
      LOGICAL FAILED
      INTEGER I
      REAL XXI(99), XTV(99), XVV(99), X(99), EXPL, XXIV(99)
      REAL AX(99), FX(99), XT(99), XV(99)
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
      DO  2 I = 1, NX
         A(I, 1) = -UX(I, 1)
         AUX(I, 1, 1) = -1
         F(I, 1) = (-UT(I, 1))-EXPL((-X(I))*V(4))*(X(I)+V(4)**2)
         FUT(I, 1, 1) = -1
         FV(I, 1, 4) = (-EXPL((-X(I))*V(4)))*(2.*V(4)+(X(I)+V(4)**2)*(-X
     1      (I)))
         FX(I) = (-EXPL((-X(I))*V(4)))*(1.-V(4)*X(I)-V(4)**3)
   2     CONTINUE
C MAP A AND F INTO XI SYSTEM.
      CALL POSTI(XI, X, XT, XXI, XV, XTV, XXIV, XVV, NX, UX, UT, NU, V
     1   , VT, NV, 1, 3, A, AX, AU, AUX, AUT, AUTX, AV, AVT, F, FX, FU
     2   , FUX, FUT, FUTX, FV, FVT)
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV, B, BU,
     1   BUX, BUT, BUTX, BV, BVT)
      INTEGER NU, NV
      REAL T, L, R, U(NU, 2), UX(NU, 2), UT(NU, 2)
      REAL UTX(NU, 2), V(NV), VT(NV), B(NU, 2), BU(NU, NU, 2), BUX(NU, 
     1   NU, 2)
      REAL BUT(NU, NU, 2), BUTX(NU, NU, 2), BV(NU, NV, 2), BVT(NU, NV, 2
     1   )
      REAL EXPL
C U(0,T) = 1
      B(1, 1) = U(1, 1)-1.
C U(1,T) = EXP(-T)
      B(1, 2) = U(1, 2)-EXPL(-V(4))
      BU(1, 1, 1) = 1
      BU(1, 1, 2) = 1
      BV(1, 4, 2) = EXPL(-V(4))
      RETURN
      END
      SUBROUTINE DEE(T, K, X, NX, U, UT, NU, NXMK, V, VT, NV, D, 
     1   DU, DUT, DV, DVT)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      REAL T, X(NX), U(NXMK, NU), UT(NXMK, NU), V(NV), VT(NV)
      REAL D(NV), DU(NV, NXMK, NU), DUT(NV, NXMK, NU), DV(NV, NV), DVT(
     1   NV, NV)
      COMMON /PARAM/ VC, XC, XI0
      REAL VC(4), XC(3), XI0
      INTEGER INTRVR, I, ILEFT
      REAL EXP, BX(10), XX(1)
      INTEGER TEMP
      D(1) = V(1)
C X(0,V) = 0.
      DV(1, 1) = 1
      XX(1) = XI0
      ILEFT = INTRVR(NX, X, XX(1))
C GET THE B-SPLINE BASIS AT XX.
      CALL BSPLN(K, X, NX, XX, 1, ILEFT, BX)
      D(2) = -EXP(-1E0)
C D(2) = W(XI0,T) - EXP(-1).
      DO  1 I = 1, K
         TEMP = ILEFT+I-K
         D(2) = D(2)+U(TEMP, 1)*BX(I)
         TEMP = ILEFT+I-K
         DU(2, TEMP, 1) = BX(I)
   1     CONTINUE
      D(3) = V(3)-1.
C X(2,V) = 1.
      DV(3, 3) = 1
      D(4) = VT(4)-1.
      DVT(4, 4) = 1
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, V0, T, U, V, NU, NXMK, NV, K, X, 
     1   NX, DT, TSTOP)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      REAL T0, U0(NXMK, NU), V0(NV), T, U(NXMK, NU), V(NV)
      REAL X(NX), DT, TSTOP
      COMMON /PARAM/ VC, XX, XI0
      REAL VC(4), XX(3), XI0
      COMMON /TIME/ TT
      REAL TT
      EXTERNAL UOFX
      INTEGER I1MACH
      REAL EU, EESFF, EV, LPLMT
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .NE. T) GOTO 2
         TEMP = I1MACH(2)
         WRITE (TEMP,  1) T, DT
   1     FORMAT (16H RESTART FOR T =, 1PE10.2, 7H   DT =, 1PE10.2)
         RETURN
C LET DT CARRY V(2) DOWN BY NO MORE THAN A FACTOR OF 10.
   2  DT = LPLMT(T, V, NV, T0, V0, 1E-1, DT)
      TT = T
C UOFX NEEDS V FOR MAPPING.
      CALL MOVEFR(NV, V, VC)
      EU = EESFF(K, X, NX, U, UOFX)
C ERROR IN POSITION OF BOUNDARY LAYER.
      EV = V(2)-1./XI0/T
      TEMP = I1MACH(2)
      WRITE (TEMP,  3) T, EU, EV
   3  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 1PE10.2, 6H   V =, 1P
     1   E10.2)
      RETURN
      END
      SUBROUTINE UOFX(XI, NX, U, W)
      INTEGER NX
      REAL XI(NX), U(NX), W(NX)
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /PARAM/ VC, X, XI0
      REAL VC(4), X(3), XI0
      COMMON /TIME/ T
      REAL T
      INTEGER IXV, IXX, ISTKGT, I, IS(1000)
      REAL EXPL, RS(1000), WS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      INTEGER TEMP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
      IXX = ISTKGT(NX, 3)
C SPACE FOR X AND XV.
      IXV = ISTKGT(3*NX, 3)
C MAP INTO USER SYSTEM.
      CALL LPLMX(XI, NX, VC, 3, WS(IXX), WS(IXV))
      DO  1 I = 1, NX
         TEMP = IXX+I
         U(I) = EXPL((-WS(TEMP-1))*T)
   1     CONTINUE
      CALL LEAVE
      RETURN
      END
