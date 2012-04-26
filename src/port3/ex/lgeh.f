C$TEST LGEH
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LGEH
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM GENM
C
C***********************************************************************
        INTEGER I, J, L, N, IA, IWRITE, I1MACH
        REAL A(50, 50), AA(50, 50), B(50), X(50)
        REAL RELERR, RELRES, XNORM, RNORM, ERR, R(50)
        REAL GENM, SAMAX
        IA = 50
C
C GENERATE MATRIX
C
        N=50
        DO 20 I=1,N
           DO 10 J=I,N
              A(I,J)=J-I
              A(J,I)=J-I + 1
              AA(I,J)=A(I,J)
              AA(J,I)=A(J,I)
 10        CONTINUE
           B(I)=I
 20     CONTINUE
C
C GENERATE RIGHT HAND SIDE
C
        CALL GEML(N,A,IA,B,X)
C
C MAKE COPY OF RIGHT HAND SIDE
C
        CALL MOVEFR(N,X,B)
C
C SOLVE THE SYSTEM
C
        CALL GELE(N,A,IA,B,N,1)
C
C COMPUTE THE RELATIVE ERROR AND THE RELATIVE RESIDUAL
C
        CALL GEML(N,AA,IA,B,R)
        ERR=0.0
        DO 30 I=1,N
           ERR=AMAX1(ERR,ABS(B(I)-FLOAT(I)))
           R(I)=R(I)-X(I)
 30     CONTINUE
        XNORM=SAMAX(N,X,1)
        RNORM=SAMAX(N,R,1)
        RELERR=ERR/XNORM
        RELRES=RNORM/(XNORM*GENM(N,AA,IA))
        IWRITE=I1MACH(2)
        WRITE(IWRITE,31)RELERR,RELRES
 31     FORMAT(16H RELATIVE ERROR=,E15.5,19H RELATIVE RESIDUAL=,
     1   E15.5)
        STOP
        END
