C$TEST LBAF
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LBAF
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BALE
C
C***********************************************************************
         INTEGER IG, IWRITE, I1MACH, N, ML, II, MP1, I, K
         INTEGER IB, NB, IT, ILAPSZ
         REAL G(19, 100), B(100, 10), BB(100, 10), GG(19, 100)
         REAL COND, TIME1, TIME2, UNI
C
C THIS PROGRAM SOLVES BANDED SYSTEMS USING BALE AND
C BASS AND COMPARES THE TIME FOR EACH OF THEM. THE
C SYSTEMS HAVE VARIOUS BANDWIDTHS,DIMENSIONS, AND
C NUMBERS OF RIGHT-HAND SIDES
         DOUBLE PRECISION D(600)
         COMMON /CSTAK/ D
C MAKE SURE THE STACK MECHANISM HAS SUFFICIENT SPACE
C FOR BASS
         CALL ISTKIN(1200,3)
         IG=19
         IWRITE=I1MACH(2)
         IB=100
         DO 70 N=50,100,50
            DO 60 ML=2,10,8
               M=2*ML - 1
               MP1=M+1
               DO 50 NB=1,10,9
                  WRITE(IWRITE,1)N,M,NB
  1               FORMAT(/5H N IS,I4,6H M IS ,I3,7H NB IS ,I3)
C
C CONSTRUCT THE MATRIX A(I,J)=ABS(I-J) AND PACK IT INTO G
C AND MAKE A COPY OF THE MATRIX SO THE SYSTEM CAN BE
C SOLVED WITH BOTH BALE AND BASS
C
                  K=ML - 1
                  DO 20 I=1,ML
                     II=MP1 - I
                     DO 10 J=1,N
                        G(I,J)=K
                        GG(I,J)=K
                        G(II,J)=K
                        GG(II,J)=K
  10                 CONTINUE
                     K=K - 1
  20              CONTINUE
C
C CONSTRUCT RANDOM RIGHT-HAND SIDES
C AND MAKE A COPY
C
                  DO 40 I=1,NB
                     DO 30 II=1,N
                        B(II,I)=UNI(0)
                        BB(II,I)=B(II,I)
  30                 CONTINUE
  40              CONTINUE
C
C SOLVE THE SYSTEM USING BOTH BASS AND BALE
C
                  IT=ILAPSZ(0)
                  CALL BASS(N,ML,M,G,IG,B,IB,NB,COND)
                  TIME1=(ILAPSZ(0)-IT)/64.0
                  WRITE(IWRITE,41)TIME1
  41              FORMAT(34H TIME FOR BASS IN MILLISECONDS IS ,F10.1)
                  IT=ILAPSZ(0)
                  CALL BALE(N,ML,M,GG,IG,BB,IB,NB)
                  TIME2=(ILAPSZ(0)-IT)/64.0
                  WRITE(IWRITE,42)TIME2
  42              FORMAT(34H TIME FOR BALE IN MIILISECONDS IS ,F10.1)
  50           CONTINUE
  60        CONTINUE
  70     CONTINUE
         STOP
         END
      INTEGER FUNCTION ILAPSZ(N)
      INTEGER N
      ILAPSZ = 0
      RETURN
      END
