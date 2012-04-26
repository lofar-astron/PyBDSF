C$TEST NP2E
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE NP2E
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SN2G
C
C***********************************************************************
        INTEGER N,P
        EXTERNAL OSBN, OSBNJ
        REAL Y(10),YY(10),T(10),X(5)
        COMMON /YT/YY,T
C GENERATE DATA FOR PROBLEM
        DATA Y(1)/8.44E-1/, Y(2) /9.36E-1/, Y(3) /8.81E-1/
     1  Y(4)/7.84E-1/, Y(5)/ 6.85E-1/, Y(6)/6.03E-1/,
     2  Y(7) /5.38E-1/ , Y(8) /4.90E-1/, Y(9)/4.57E-1/
        P=5
        N=9
        DO 10 I=1,9
           YY(I) = Y(I)
           T(I)=-30.E0*FLOAT(I-1)
 10     CONTINUE
C INITIALIZE X
        X(1)=0.5
        X(2)=1.5
        X(3)=-1.
        X(4)=.01
        X(5)=.02
C
C SOLVE THE PROBLEM
C
        CALL SN2G(N, P, X, OSBN, OSBNJ, 100, 1.E-4)
C       PRINT RESULTS ON STANDARD OUTPUT UNIT
        IWRITE = I1MACH(2)
        WRITE(IWRITE, 20)(X(I),I=1,P)
 20     FORMAT(10H SOLUTION-,5E15.5)
        STOP
        END
        SUBROUTINE OSBN(N,P,X,NF,R)
C THIS SUBROUTINE COMPUTES THE MODEL
        INTEGER P, N, NF
        REAL X(P), R(N)
        REAL Y(10), T(10)
        COMMON /YT/ Y, T
        DO 10 I=1,N
           R(I)=Y(I)-(X(1)+X(2)*EXP(X(4)*T(I))+X(3)*EXP(X(5)*T(I)))
 10     CONTINUE
        RETURN
        END
        SUBROUTINE OSBNJ(N,P,X,NF,J)
C THIS SUBROUTINE COMPUTES THE JACOBIAN OF THE MODEL
        INTEGER P, N, NF
        REAL X(P), J(N,P)
        REAL Y(10), T(10)
        COMMON /YT/ Y, T
        DO 10 I=1,N
           J(I,1)=-1.0E0
           J(I,2)=-EXP(X(4)*T(I))
           J(I,3)=-EXP(X(5)*T(I))
           J(I,4)=-T(I)*X(2)*EXP(X(4)*T(I))
           J(I,5)=-T(I)*X(3)*EXP(X(5)*T(I))
 10     CONTINUE
        RETURN
        END
