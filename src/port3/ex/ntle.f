C$TEST NTLE
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE NTLE
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SMNFB
C
C***********************************************************************
        INTEGER N
        EXTERNAL ROSN
        REAL X(2), B(2,2)
        N=2
C INITIALIZE X
        X(1)=-1.2
        X(2)=1.0
C SET UP THE BOUND ARRAY
C R1MACH(2) CONTAINS THE LARGEST NUMBER IN THE MACHINE
        B(1,1)=-R1MACH(2)
        B(2,1)=0.5
        B(1,2)=0.0
        B(2,2)=1.0
C
C SOLVE THE PROBLEM
C
        CALL SMNFB(N, X, B, ROSN, 100, 1.E-4)
C       PRINT RESULTS ON STANDARD OUTPUT UNIT
        IWRITE=I1MACH(2)
        WRITE(IWRITE,10)(X(I),I=1,N)
 10     FORMAT(10H SOLUTION-,5E15.5)
        STOP
        END
        SUBROUTINE ROSN(N,X,NF,F)
C THIS SUBROUTINE COMPUTES THE  FUNCTION
        INTEGER N, NF
        REAL X(N), F
        F=100.0*(X(2)-X(1)*X(1))**2 + (1.0 - X(1))**2
        RETURN
        END
