C$TEST MFTG
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE MFTG
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM MFTCR
C
C***********************************************************************
      REAL A(200,100)
      REAL AA(200,100)
      REAL TX(200),TY(200)
      REAL RSTAK(3212)
      REAL FN2,RMSERR,SIGN,SUM
      INTEGER IFX(25),IFY(25)
      INTEGER I,J,N,N2,NP2,NP3,N2MK,NNS0,NNS,NSEGS
      DOUBLE PRECISION DSTAK(1606)
      COMMON /CSTAK/DSTAK
      EQUIVALENCE (RSTAK(1),DSTAK(1))
C
C
      CALL ISTKIN(3212,3)
C
      N   = 100
      NP2 = 102
      N2  = 200
C
C   THE SEGMENT SIZE IS ARBITRARILY CHOSEN TO BE 16 X N - I.E. USING
C   MFTCC, MFTRC, MFTCR TO COMPUTE UP TO 16 INDEPENDENT VECTORS AT A
C   TIME.
C
      NSEGS = (N-1)/16 + 1
      NNS0  = MOD(N-1,16) + 1
C
C   EXAMPLE USES RANDOM INPUT DATA
C
      DO 1 J = 1,N
         DO 2 I = 1,N
            A(I,J)  = UNI(0)
            AA(I,J) = A(I,J)
 2       CONTINUE
 1    CONTINUE
C
      SIGN = 1.0E0
      CALL MFTRI(N,IFX,TX)
      CALL MFTCI(N,IFY,TY)
C
C  X-DIMENSION
C
      NNS = NNS0
      L   = 1
      DO 3 LL = 1,NSEGS
         CALL MFTRC(N,NNS,A(1,L),1,N2,A(1,L),A(2,L),2,N2,IFX,TX,SIGN)
         L = L + NNS
         NNS = 16
 3    CONTINUE
C
C  FILL-IN FROM CONJUGATION OF TERMS
C
      NP3  = N+3
      N2MK = N-1
      DO 4 I = NP3,N2,2
         DO 5 J = 1,N
            A(I,J)   = A(N2MK,J)
            A(I+1,J) = - A(N2MK+1,J)
 5       CONTINUE
         N2MK = N2MK - 2
 4    CONTINUE
C
C  DO COMPLEX PART IN Y-DIRECTION
C
      NNS = NNS0
      L   = 1
      DO 6 LL = 1,NSEGS
         CALL MFTCC(N,NNS,A(L,1),A(L+1,1),N2,2,
     *              A(L,1),A(L+1,1),N2,2,IFY,TY,SIGN)
         L   = L + 2*NNS
         NNS = 16
 6    CONTINUE
C
C  NOW GO BACKWARDS, COMPLEX TO COMPLEX FIRST
C
      SIGN = -1.0E0
      NNS  = NNS0
      L    = 1
      DO 7 LL = 1,NSEGS
         CALL MFTCC(N,NNS,A(L,1),A(L+1,1),N2,2,
     *              A(L,1),A(L+1,1),N2,2,IFY,TY,SIGN)
         L   = L + 2*NNS
         NNS = 16
 7    CONTINUE
C
C  AND BACK TO REAL
C
      NNS  = NNS0
      L    = 1
      DO 8 LL = 1,NSEGS
         CALL MFTCR(N,NNS,A(1,L),A(2,L),2,N2,A(1,L),1,N2,IFX,TX,SIGN)
         L   = L + NNS
         NNS = 16
 8    CONTINUE
C
C  COMPARE TO INPUT
C
      FN2 = 1./FLOAT(N*N)
      SUM = 0.0E0
      DO 9 J = 1,N
         DO 10 I = 1,N
            SUM = SUM + (AA(I,J) - FN2*A(I,J))**2
 10      CONTINUE
 9    CONTINUE
C
C   PRINT ROOT MEAN SQUARE ERROR
C
      RMSERR = SQRT(SUM*FN2)
      IWRITE   = I1MACH(2)
      WRITE(IWRITE,1000) N,N,RMSERR
 1000 FORMAT(1X,5H FOR ,I3,1HX,I3,20H ARRAY, RMS ERROR = ,1PE12.3)
      STOP
      END
