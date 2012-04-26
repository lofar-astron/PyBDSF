C$TEST DTG2
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DTG2
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM DTTGR
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(350000)
      EXTERNAL HANDLE, BC, AF
      INTEGER NDX, NDY, ISTKGT, IS(1000), IU, IX
      INTEGER IY, NU, KX, NX, KY, NY
      INTEGER IDUMB
      REAL ERRPAR(2), RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION TSTART, DT, LX, LY, RX, RY
      DOUBLE PRECISION WS(500), TSTOP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO SOLVE TWO COUPLED, NONLINEAR HEAT EQUATIONS.
C   U1 SUB T = DIV . ( U1X, U1Y ) - U1*U2 + G1
C   U2 SUB T = DIV . ( U2X, U2Y ) - U1*U2 + G2
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(350000, 4)
      CALL ENTER(1)
      NU = 2
      LX = 0
      RX = 1
      LY = 0
      RY = 1
      KX = 4
      KY = 4
      NDX = 3
      NDY = 3
      TSTART = 0
      TSTOP = 1
      DT = 1E-2
      ERRPAR(1) = 1E-2
      ERRPAR(2) = 1E-4
C UNIFORM GRID.
      IX = IDUMB(LX, RX, NDX, KX, NX)
C UNIFORM GRID.
      IY = IDUMB(LY, RY, NDY, KY, NY)
C SPACE FOR THE SOLUTION.
      IU = ISTKGT(NU*(NX-KX)*(NY-KY), 4)
      CALL SETD(NU*(NX-KX)*(NY-KY), 1D0, WS(IU))
      CALL DTTGR(WS(IU), NU, KX, WS(IX), NX, KY, WS(IY), NY, TSTART, 
     1   TSTOP, DT, AF, BC, ERRPAR, HANDLE)
      CALL LEAVE
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, X, NX, Y, NY, NU, U, UT, UX, UY, UXT, UYT
     1   , A, AU, AUT, AUX, AUY, AUXT, AUYT, F, FU, FUT, FUX, FUY, FUXT,
     2   FUYT)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, X(NX), Y(NY), U(NX, NY, NU), UT(NX, NY, NU), 
     1   UX(NX, NY, NU)
      DOUBLE PRECISION UY(NX, NY, NU), UXT(NX, NY, NU), UYT(NX, NY, NU),
     1   A(NX, NY, NU, 2), AU(NX, NY, NU, NU, 2), AUT(NX, NY, NU, NU, 2)
      DOUBLE PRECISION AUX(NX, NY, NU, NU, 2), AUY(NX, NY, NU, NU, 2), 
     1   AUXT(NX, NY, NU, NU, 2), AUYT(NX, NY, NU, NU, 2), F(NX, NY, NU)
     2   , FU(NX, NY, NU, NU)
      DOUBLE PRECISION FUT(NX, NY, NU, NU), FUX(NX, NY, NU, NU), FUY(NX,
     1   NY, NU, NU), FUXT(NX, NY, NU, NU), FUYT(NX, NY, NU, NU)
      INTEGER P, Q
      DOUBLE PRECISION DEXP
      DO  2 Q = 1, NY
         DO  1 P = 1, NX
            A(P, Q, 1, 1) = UX(P, Q, 1)
            AUX(P, Q, 1, 1, 1) = 1
            A(P, Q, 1, 2) = UY(P, Q, 1)
            AUY(P, Q, 1, 1, 2) = 1
            F(P, Q, 1) = UT(P, Q, 1)+U(P, Q, 1)*U(P, Q, 2)
            FU(P, Q, 1, 1) = U(P, Q, 2)
            FU(P, Q, 1, 2) = U(P, Q, 1)
            FUT(P, Q, 1, 1) = 1
            A(P, Q, 2, 1) = UX(P, Q, 2)
            AUX(P, Q, 2, 2, 1) = 1
            A(P, Q, 2, 2) = UY(P, Q, 2)
            AUY(P, Q, 2, 2, 2) = 1
            F(P, Q, 2) = UT(P, Q, 2)+U(P, Q, 1)*U(P, Q, 2)
            FU(P, Q, 2, 1) = U(P, Q, 2)
            FU(P, Q, 2, 2) = U(P, Q, 1)
            FUT(P, Q, 2, 2) = 1
            F(P, Q, 1) = F(P, Q, 1)-(DEXP(T*(X(P)-Y(Q)))*(X(P)-Y(Q)-2D0*
     1         T*T)+1D0)
            F(P, Q, 2) = F(P, Q, 2)-(DEXP(T*(Y(Q)-X(P)))*(Y(Q)-X(P)-2D0*
     1         T*T)+1D0)
   1        CONTINUE
   2     CONTINUE
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
      DOUBLE PRECISION DEXP
      DO  2 J = 1, NY
         DO  1 I = 1, NX
            BU(I, J, 1, 1) = 1
            B(I, J, 1) = U(I, J, 1)-DEXP(T*(X(I)-Y(J)))
            BU(I, J, 2, 2) = 1
            B(I, J, 2) = U(I, J, 2)-DEXP(T*(Y(J)-X(I)))
   1        CONTINUE
   2     CONTINUE
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, T, U, NV, DT, TSTOP)
      INTEGER NV
      DOUBLE PRECISION T0, U0(NV), T, U(NV), DT, TSTOP
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      COMMON /D7TGRP/ ERRPAR, NU, MXQ, MYQ
      INTEGER NU, MXQ, MYQ
      REAL ERRPAR(2)
      COMMON /D7TGRM/ KX, IX, NX, KY, IY, NY
      INTEGER KX, IX, NX, KY, IY, NY
      INTEGER IFA, ITA(2), IXA(2), NTA(2), NXA(2), IDLUMD
      INTEGER IXS, IYS, NXS, NYS, ISTKGT, I
      INTEGER J, IEWE, KA(2), MA(2), IS(1000), I1MACH
      REAL RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION DABS, ERRU, DMAX1, WS(500)
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C THE PORT LIBRARY STACK AND ITS ALIASES.
      IF (T0 .NE. T) GOTO 2
         WRITE (6,  1) T
   1     FORMAT (16H RESTART FOR T =, 1PE10.2)
         RETURN
   2  CALL ENTER(1)
C FIND THE ERROR IN THE SOLUTION AT 2*KX * 2*KY POINTS / MESH RECTANGLE.
C X SEARCH GRID.
      IXS = IDLUMD(WS(IX), NX, 2*KX, NXS)
C Y SEARCH GRID.
      IYS = IDLUMD(WS(IY), NY, 2*KY, NYS)
C U SEARCH GRID VALUES.
      IEWE = ISTKGT(NU*NXS*NYS, 4)
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
      DO  5 J = 1, NU
C EVALUATE THEM.
         TEMP = (J-1)*(NX-KX)*(NY-KY)
         CALL DTSD1(2, KA, WS, ITA, NTA, U(TEMP+1), WS, IXA, NXA, MA, 
     1      WS(IFA))
C ERROR IN SOLUTION VALUES.
         ERRU = 0
         TEMP = NXS*NYS
         DO  3 I = 1, TEMP
            TEMP2 = IEWE+I-1+(J-1)*NXS*NYS
            TEMP1 = IFA+I
            ERRU = DMAX1(ERRU, DABS(WS(TEMP2)-WS(TEMP1-1)))
   3        CONTINUE
         TEMP = I1MACH(2)
         WRITE (TEMP,  4) T, J, ERRU
   4     FORMAT (14H ERROR IN U(.,, 1PE10.2, 1H,, I2, 3H) =, 1PE10.2)
   5     CONTINUE
      CALL LEAVE
      RETURN
      END
      SUBROUTINE EWE(T, X, NX, Y, NY, U, NU)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, X(NX), Y(NY), U(NX, NY, NU)
      INTEGER I, J, P
      REAL FLOAT
      DOUBLE PRECISION DBLE, DEXP
C THE EXACT SOLUTION.
      DO  3 P = 1, NU
         DO  2 I = 1, NX
            DO  1 J = 1, NY
               U(I, J, P) = DEXP(DBLE(FLOAT((-1)**(P+1)))*T*(X(I)-Y(J)))
   1           CONTINUE
   2        CONTINUE
   3     CONTINUE
      RETURN
      END
