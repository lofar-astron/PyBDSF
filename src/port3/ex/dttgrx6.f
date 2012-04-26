C$TEST DTG6
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DTG6
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM DTTGR
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(350000)
      EXTERNAL HANDLE, BC, AF
      INTEGER IUE, NDX, NDY, IUR, IXR, IYR
      INTEGER NXR, NYR, ISTKGT, I, IS(1000), IU
      INTEGER IX, IY, NU, KX, NX, KY
      INTEGER NY, I1MACH
      REAL ERRPAR(2), RS(1000), FLOAT
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION TSTART, DBLE, DABS, EERR, ERRE, ERRR
      DOUBLE PRECISION DMAX1, DT, LX, LY, RX, RY
      DOUBLE PRECISION WS(500), TSTOP
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO GET ERROR ESTIMATES FOR LAPLACES EQUATION WITH REAL ( Z*LOG(Z) ) AS
C SOLUTION.
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(350000, 4)
      CALL ENTER(1)
      NU = 1
      LX = 0
      RX = 1
      LY = 0
      RY = 1
      KX = 4
      KY = 4
      NDX = 2
      NDY = 2
      TSTART = 0
      TSTOP = 1
      DT = 1
      ERRPAR(1) = 1E-2
      ERRPAR(2) = 1E-4
      NX = NDX+2*(KX-1)
C SPACE FOR X MESH.
      IX = ISTKGT(NX, 4)
      DO  1 I = 1, KX
         TEMP = IX+I
         WS(TEMP-1) = 0
         TEMP = IX+NX-I
         WS(TEMP) = RX
   1     CONTINUE
C 0 AND RX MULT = KX.
      TEMP = NDX-1
      DO  2 I = 1, TEMP
         TEMP2 = IX+KX-2+I
         WS(TEMP2) = RX*(DBLE(FLOAT(I-1))/(DBLE(FLOAT(NDX))-1D0))**KX
   2     CONTINUE
      NY = NDY+2*(KY-1)
C SPACE FOR Y MESH.
      IY = ISTKGT(NY, 4)
      DO  3 I = 1, KY
         TEMP = IY+I
         WS(TEMP-1) = 0
         TEMP = IY+NY-I
         WS(TEMP) = RY
   3     CONTINUE
C 0 AND RY MULT = KY.
      TEMP = NDY-1
      DO  4 I = 1, TEMP
         TEMP2 = IY+KY-2+I
         WS(TEMP2) = RY*(DBLE(FLOAT(I-1))/(DBLE(FLOAT(NDY))-1D0))**KY
   4     CONTINUE
C SPACE FOR THE SOLUTION.
      IU = ISTKGT(NU*(NX-KX)*(NY-KY), 4)
      CALL SETD(NU*(NX-KX)*(NY-KY), 0D0, WS(IU))
      TEMP = I1MACH(2)
      WRITE (TEMP,  5) 
   5  FORMAT (23H SOLVING ON CRUDE MESH.)
      CALL DTTGR(WS(IU), NU, KX, WS(IX), NX, KY, WS(IY), NY, TSTART, 
     1   TSTOP, DT, AF, BC, ERRPAR, HANDLE)
      DT = 1
      NDX = 2*NDX-1
C REFINE MESH.
      NDY = 2*NDY-1
      NXR = NDX+2*(KX-1)
C SPACE FOR X MESH.
      IXR = ISTKGT(NXR, 4)
      DO  6 I = 1, KX
         TEMP = IXR+I
         WS(TEMP-1) = 0
         TEMP = IXR+NXR-I
         WS(TEMP) = RX
   6     CONTINUE
C 0 AND RX MULT = KX.
      TEMP = NDX-1
      DO  7 I = 1, TEMP
         TEMP2 = IXR+KX-2+I
         WS(TEMP2) = RX*(DBLE(FLOAT(I-1))/(DBLE(FLOAT(NDX))-1D0))**KX
   7     CONTINUE
      NYR = NDY+2*(KY-1)
C SPACE FOR Y MESH.
      IYR = ISTKGT(NYR, 4)
      DO  8 I = 1, KY
         TEMP = IYR+I
         WS(TEMP-1) = 0
         TEMP = IYR+NYR-I
         WS(TEMP) = RY
   8     CONTINUE
C 0 AND RY MULT = KY.
      TEMP = NDY-1
      DO  9 I = 1, TEMP
         TEMP2 = IYR+KY-2+I
         WS(TEMP2) = RY*(DBLE(FLOAT(I-1))/(DBLE(FLOAT(NDY))-1D0))**KY
   9     CONTINUE
C SPACE FOR THE SOLUTION.
      IUR = ISTKGT(NU*(NXR-KX)*(NYR-KY), 4)
      CALL SETD(NU*(NXR-KX)*(NYR-KY), 0D0, WS(IUR))
      TEMP = I1MACH(2)
      WRITE (TEMP,  10) 
  10  FORMAT (25H SOLVING ON REFINED MESH.)
      CALL DTTGR(WS(IUR), NU, KX, WS(IXR), NXR, KY, WS(IYR), NYR, 
     1   TSTART, TSTOP, DT, AF, BC, ERRPAR, HANDLE)
      DT = 1
      ERRPAR(1) = ERRPAR(1)/10.
      ERRPAR(2) = ERRPAR(2)/10.
C SPACE FOR THE SOLUTION.
      IUE = ISTKGT(NU*(NX-KX)*(NY-KY), 4)
      CALL SETD(NU*(NX-KX)*(NY-KY), 0D0, WS(IUE))
      TEMP = I1MACH(2)
      WRITE (TEMP,  11) 
  11  FORMAT (24H SOLVING WITH ERRPAR/10.)
      CALL DTTGR(WS(IUE), NU, KX, WS(IX), NX, KY, WS(IY), NY, TSTART, 
     1   TSTOP, DT, AF, BC, ERRPAR, HANDLE)
      ERRR = EERR(KX, IX, NX, KY, IY, NY, WS(IU), NU, IXR, NXR, IYR, 
     1   NYR, WS(IUR), TSTOP)
      ERRE = 0
      TEMP = NU*(NX-KX)*(NY-KY)
      DO  12 I = 1, TEMP
         TEMP2 = IU+I
         TEMP1 = IUE+I
         ERRE = DMAX1(ERRE, DABS(WS(TEMP2-1)-WS(TEMP1-1)))
  12     CONTINUE
      TEMP = I1MACH(2)
      WRITE (TEMP,  13) ERRE
  13  FORMAT (24H U ERROR FROM U AND UE =, 1PE10.2)
      TEMP = I1MACH(2)
      WRITE (TEMP,  14) ERRR
  14  FORMAT (24H U ERROR FROM U AND UR =, 1PE10.2)
      CALL LEAVE
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, XI, NX, YI, NY, NU, U, UT, UX, UY, UXT, 
     1   UYT, A, AU, AUT, AUX, AUY, AUXT, AUYT, F, FU, FUT, FUX, FUY, 
     2   FUXT, FUYT)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, XI(NX), YI(NY), U(NX, NY, NU), UT(NX, NY, NU),
     1   UX(NX, NY, NU)
      DOUBLE PRECISION UY(NX, NY, NU), UXT(NX, NY, NU), UYT(NX, NY, NU),
     1   A(NX, NY, NU, 2), AU(NX, NY, NU, NU, 2), AUT(NX, NY, NU, NU, 2)
      DOUBLE PRECISION AUX(NX, NY, NU, NU, 2), AUY(NX, NY, NU, NU, 2), 
     1   AUXT(NX, NY, NU, NU, 2), AUYT(NX, NY, NU, NU, 2), F(NX, NY, NU)
     2   , FU(NX, NY, NU, NU)
      DOUBLE PRECISION FUT(NX, NY, NU, NU), FUX(NX, NY, NU, NU), FUY(NX,
     1   NY, NU, NU), FUXT(NX, NY, NU, NU), FUYT(NX, NY, NU, NU)
      INTEGER I, P, Q
      DO  3 I = 1, NU
         DO  2 Q = 1, NY
            DO  1 P = 1, NX
               A(P, Q, I, 1) = UX(P, Q, I)
               A(P, Q, I, 2) = UY(P, Q, I)
               AUX(P, Q, I, I, 1) = 1
               AUY(P, Q, I, I, 2) = 1
   1           CONTINUE
   2        CONTINUE
   3     CONTINUE
      RETURN
      END
      SUBROUTINE BC(T, X, NX, Y, NY, LX, RX, LY, RY, U, UT, UX, 
     1   UY, UXT, UYT, NU, B, BU, BUT, BUX, BUY, BUXT, BUYT)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, X(NX), Y(NY), LX, RX, LY
      DOUBLE PRECISION RY, U(NX, NY, NU), UT(NX, NY, NU), UX(NX, NY, NU)
     1   , UY(NX, NY, NU), UXT(NX, NY, NU)
      DOUBLE PRECISION UYT(NX, NY, NU), B(NX, NY, NU), BU(NX, NY, NU, 
     1   NU), BUT(NX, NY, NU, NU), BUX(NX, NY, NU, NU), BUY(NX, NY, NU
     2   , NU)
      DOUBLE PRECISION BUXT(NX, NY, NU, NU), BUYT(NX, NY, NU, NU)
      INTEGER I, J
      DOUBLE PRECISION R, DCOS, DLOG, DSIN, DATAN, THETA
      DOUBLE PRECISION DSQRT
      DO  6 J = 1, NY
         DO  5 I = 1, NX
            IF (Y(J) .NE. LY) GOTO 1
               B(I, J, 1) = UY(I, J, 1)
C NEUMANN DATA ON BOTTOM.
               BUY(I, J, 1, 1) = 1
               GOTO  4
   1           R = DSQRT(X(I)**2+Y(J)**2)
C DIRICHLET DATA.
               IF (X(I) .LE. 0D0) GOTO 2
                  THETA = DATAN(Y(J)/X(I))
                  GOTO  3
   2              THETA = 2D0*DATAN(1D0)
   3           B(I, J, 1) = U(I, J, 1)-R*(DCOS(THETA)*DLOG(R)-THETA*
     1            DSIN(THETA))
               BU(I, J, 1, 1) = 1
   4        CONTINUE
   5        CONTINUE
   6     CONTINUE
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, T, U, NV, DT, TSTOP)
      INTEGER NV
      DOUBLE PRECISION T0, U0(NV), T, U(NV), DT, TSTOP
      COMMON /D7TGRP/ ERRPAR, NU, MXQ, MYQ
      INTEGER NU, MXQ, MYQ
      REAL ERRPAR(2)
      COMMON /D7TGRM/ KX, IX, NX, KY, IY, NY
      INTEGER KX, IX, NX, KY, IY, NY
      IF (T0 .NE. T) GOTO 2
         WRITE (6,  1) T
   1     FORMAT (16H RESTART FOR T =, 1PE10.2)
         RETURN
   2  CALL GERR(KX, IX, NX, KY, IY, NY, U, NU, T)
      RETURN
      END
      SUBROUTINE GERR(KX, IX, NX, KY, IY, NY, U, NU, T)
      INTEGER KX, IX, NX, KY, IY, NY
      INTEGER NU
      DOUBLE PRECISION U(1), T
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IFA, ITA(2), IXA(2), NTA(2), NXA(2), IDLUMD
      INTEGER IXS, IYS, NXS, NYS, ISTKGT, I
      INTEGER IEWE, KA(2), MA(2), IS(1000), I1MACH
      REAL RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION DABS, ERRU, DMAX1, WS(500)
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO GET AND PRINT THE ERROR AT EACH TIME-STEP.
C U(NX-KX,NY-KY,NU).
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
C FIND THE ERROR IN THE SOLUTION AT 2*KX * 2*KY POINTS / MESH RECTANGLE.
C X SEARCH GRID.
      IXS = IDLUMD(WS(IX), NX, 2*KX, NXS)
C Y SEARCH GRID.
      IYS = IDLUMD(WS(IY), NY, 2*KY, NYS)
C U SEARCH GRID VALUES.
      IEWE = ISTKGT(NXS*NYS, 4)
C THE EXACT SOLUTION.
      CALL EWE(T, WS(IXS), NXS, WS(IYS), NYS, WS(IEWE), NU)
      KA(1) = KX
      KA(2) = KY
      ITA(1) = IX
      ITA(2) = IY
      NTA(1) = NX
      NTA(2) = NY
      IXA(1) = IXS
      IXA(2) = IYS
      NXA(1) = NXS
      NXA(2) = NYS
      MA(1) = 0
C GET SOLUTION.
      MA(2) = 0
C APPROXIMATE SOLUTION VALUES.
      IFA = ISTKGT(NXS*NYS, 4)
C EVALUATE THEM.
      CALL DTSD1(2, KA, WS, ITA, NTA, U, WS, IXA, NXA, MA, WS(IFA))
C ERROR IN SOLUTION VALUES.
      ERRU = 0
      TEMP = NXS*NYS
      DO  1 I = 1, TEMP
         TEMP2 = IEWE+I
         TEMP1 = IFA+I
         ERRU = DMAX1(ERRU, DABS(WS(TEMP2-1)-WS(TEMP1-1)))
   1     CONTINUE
      TEMP = I1MACH(2)
      WRITE (TEMP,  2) T, ERRU
   2  FORMAT (14H ERROR IN U(.,, 1PE10.2, 3H) =, 1PE10.2)
      CALL LEAVE
      RETURN
      END
      DOUBLE PRECISION FUNCTION EERR(KX, IX, NX, KY, IY, NY, U, 
     1   NU, IXR, NXR, IYR, NYR, UR, T)
      INTEGER KX, IX, NX, KY, IY, NY
      INTEGER NU, IXR, NXR, IYR, NYR
      DOUBLE PRECISION U(1), UR(1), T
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IFA, ITA(2), IXA(2), NTA(2), NXA(2), IDLUMD
      INTEGER IXS, IYS, NXS, NYS, ISTKGT, I
      INTEGER IFAR, KA(2), MA(2), IS(1000), I1MACH
      REAL RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION DABS, ERRU, DMAX1, WS(500)
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO GET AND PRINT THE ERROR ESTIMATE AT EACH TIME-STEP.
C U(NX-KX,NY-KY,NU), UR(NXR-KX,NYR-KY,NU).
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
C FIND THE ERROR IN THE SOLUTION AT 2*KX * 2*KY POINTS / FINE MESH RECTA
CNGLE.
C X SEARCH GRID.
      IXS = IDLUMD(WS(IXR), NXR, 2*KX, NXS)
C Y SEARCH GRID.
      IYS = IDLUMD(WS(IYR), NYR, 2*KY, NYS)
      KA(1) = KX
      KA(2) = KY
      ITA(1) = IX
      ITA(2) = IY
      NTA(1) = NX
      NTA(2) = NY
      IXA(1) = IXS
      IXA(2) = IYS
      NXA(1) = NXS
      NXA(2) = NYS
      MA(1) = 0
C GET SOLUTION.
      MA(2) = 0
C APPROXIMATE SOLUTION VALUES.
      IFA = ISTKGT(NXS*NYS, 4)
C EVALUATE THEM.
      CALL DTSD1(2, KA, WS, ITA, NTA, U, WS, IXA, NXA, MA, WS(IFA))
      KA(1) = KX
      KA(2) = KY
      ITA(1) = IXR
      ITA(2) = IYR
      NTA(1) = NXR
      NTA(2) = NYR
      IXA(1) = IXS
      IXA(2) = IYS
      NXA(1) = NXS
      NXA(2) = NYS
      MA(1) = 0
C GET SOLUTION.
      MA(2) = 0
C APPROXIMATE SOLUTION VALUES.
      IFAR = ISTKGT(NXS*NYS, 4)
C EVALUATE THEM.
      CALL DTSD1(2, KA, WS, ITA, NTA, UR, WS, IXA, NXA, MA, WS(IFAR))
C ERROR IN SOLUTION VALUES.
      ERRU = 0
      TEMP = NXS*NYS
      DO  1 I = 1, TEMP
         TEMP2 = IFAR+I
         TEMP1 = IFA+I
         ERRU = DMAX1(ERRU, DABS(WS(TEMP2-1)-WS(TEMP1-1)))
   1     CONTINUE
      CALL LEAVE
      EERR = ERRU
      RETURN
      END
      SUBROUTINE EWE(T, X, NX, Y, NY, U, NU)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, X(NX), Y(NY), U(NX, NY, NU)
      INTEGER I, J, P
      DOUBLE PRECISION R, DCOS, DLOG, DSIN, DATAN, THETA
      DOUBLE PRECISION DSQRT
C THE EXACT SOLUTION.
      DO  7 P = 1, NU
         DO  6 I = 1, NX
            DO  5 J = 1, NY
               R = DSQRT(X(I)**2+Y(J)**2)
               IF (X(I) .LE. 0D0) GOTO 1
                  THETA = DATAN(Y(J)/X(I))
                  GOTO  2
   1              THETA = 2D0*DATAN(1D0)
   2           IF (R .LE. 0D0) GOTO 3
                  U(I, J, P) = R*(DCOS(THETA)*DLOG(R)-THETA*DSIN(THETA))
                  GOTO  4
   3              U(I, J, P) = 0
   4           CONTINUE
   5           CONTINUE
   6        CONTINUE
   7     CONTINUE
      RETURN
      END
