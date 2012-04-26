C$TEST LYMK
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LYMK
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SYNM
C
C***********************************************************************
        INTEGER I, J, L, N, I1MACH, IWRITE
        REAL C(1300), CC(1300), B(50), X(50)
        REAL RELERR, RELRES, XNORM, RNORM, ERR, R(50)
        REAL SYNM, SAMAX
        L=0
C
C GENERATE MATRIX
C
        N=50
        DO 20 I=1,N
           DO 10 J=I,N
              L=L+1
              C(L)=J-I
              CC(L)=C(L)
 10        CONTINUE
           B(I)=I
 20     CONTINUE
C
C GENERATE RIGHT HAND SIDE
C
        CALL SYML(N,C,B,X)
C
C MAKE COPY OF RIGHT HAND SIDE
C
        CALL MOVEFR(N,X,B)
C
C SOLVE THE SYSTEM
C
        CALL SYLE(N,C,B,N,1)
C
C COMPUTE THE RELATIVE ERROR AND THE RELATIVE RESIDUAL
C
        CALL SYML(N,CC,B,R)
        ERR=0.0
        DO 30 I=1,N
           ERR=AMAX1(ERR,ABS(B(I)-FLOAT(I)))
           R(I)=R(I)-X(I)
 30     CONTINUE
        XNORM=SAMAX(N,X,1)
        RNORM=SAMAX(N,R,1)
        RELERR=ERR/XNORM
        RELRES=RNORM/(XNORM*SYNM(N,CC))
        IWRITE=I1MACH(2)
        WRITE(IWRITE,31)RELERR,RELRES
 31     FORMAT(16H RELATIVE ERROR=,E15.5,19H RELATIVE RESIDUAL=,
     1   E15.5)
        STOP
        END
