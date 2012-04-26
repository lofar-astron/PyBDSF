      SUBROUTINE D4SQR(K, M, N, Q, R, B, RHS)
      INTEGER M
      INTEGER N
      REAL Q(K, 1), R( 1), B(1), RHS(K)
      INTEGER I
      REAL BETA, ALPHA,U, X
C
C THIS SUBROUTINE UPDATES THE QR DECOMPOSTION WHENE A NEW
C ROW CONTAINED IN B IS ADDED TO THE MATRIX
C
       M=M+1
       MM1=M-1
C
C ZERO OUT ROW AND COLUMN OF Q MATRIX
C
       Q(M,M)=1.
       IF(M.EQ.1)RETURN
       DO 10 II=1,MM1
          Q(M,II)=0.0
          Q(II,M)=0.0
 10       CONTINUE
       X=RHS(M)
         IF (N.EQ.0) RETURN
         IS=1
         DO 20 I=1,N
         CALL SROTG(R(IS), B(I), ALPHA, BETA)
         CALL SROT(M, Q(I, 1), K, Q(M, 1), K, ALPHA, BETA)
         U=RHS(I)
         RHS(I)=ALPHA*U+BETA*X
         X=-BETA*U+ALPHA*X
         IS=IS+I+1
         IF (N-I.GE.1)
     1    CALL SROT2(N-I,R(IS-1),I+1,B(I+1),-1,ALPHA,BETA)
 20     CONTINUE
      RHS(M)=X
      RETURN
      END
