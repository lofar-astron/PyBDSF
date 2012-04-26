C$TEST NTLT
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE NTLT
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SNSFB
C
C***********************************************************************
        INTEGER N,P,L
        INTEGER LP1, IINC,INC(4,2)
        EXTERNAL OSBA
        REAL Y(10),T(10),X(2),C(3),B(2,2)
        COMMON /TT/T
C GENERATE DATA FOR PROBLEM
        DATA Y(1)/8.44E-1/, Y(2) /9.36E-1/, Y(3) /8.81E-1/
     1  Y(4)/7.84E-1/, Y(5)/ 6.85E-1/, Y(6)/6.03E-1/,
     2  Y(7) /5.38E-1/ , Y(8) /4.90E-1/, Y(9)/4.57E-1/
        P=2
        N=9
        L=3
        DO 10 I=1,9
           T(I)=-30.E0*FLOAT(I-1)
 10     CONTINUE
C INITIALIZE X
        X(1)=.01
        X(2)=.03
C GENERATE THE INCIDENCE MATRIX
        LP1=L+1
        DO 30 J=1,P
           DO 20 I=1,LP1
              INC(I,J)=0
 20        CONTINUE
 30     CONTINUE
        INC(2,1)=1
        INC(3,2)=1
        IINC=LP1
C SUPPLY BOUNDS
        B(1,1)=-R1MACH(2)
        B(2,1)=0.125
        B(1,2)=.03
        B(2,2)=R1MACH(2)
C
C SOLVE THE PROBLEM
C
        CALL SNSFB(N, P, L, X, B, C, Y, OSBA, INC, IINC, 100, 1.E-4)
C       PRINT RESULTS ON STANDARD OUTPUT UNIT
        IWRITE = I1MACH(2)
        WRITE(IWRITE, 40)(X(I),I=1,P)
 40     FORMAT(22H NONLINEAR PARAMETERS-,2E15.5)
        WRITE(IWRITE, 50)(C(I),I=1,L)
 50     FORMAT(19H LINEAR PARAMETERS-, 3E15.5)
        STOP
        END
        SUBROUTINE OSBA(N,P,L,X,NF,A)
C THIS SUBROUTINE COMPUTES THE MODEL
        INTEGER P, N, NF, L
        REAL X(P), A(N,L)
        REAL  T(10)
        COMMON /TT/  T
        DO 10 I=1,N
           A(I,1)=1.0
           A(I,2)=EXP(X(1)*T(I))
           A(I,3)=EXP(X(2)*T(I))
 10     CONTINUE
        RETURN
        END
