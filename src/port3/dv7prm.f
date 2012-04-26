      SUBROUTINE DV7PRM(N, IP, X)
C
C     PERMUTE X SO THAT X.OUTPUT(IP(I)) = X.INPUT(I).
C     IP IS UNCHANGED ON OUTPUT.
C
      INTEGER N
      INTEGER IP(N)
      DOUBLE PRECISION X(N)
C
      INTEGER I, J, K
      DOUBLE PRECISION S, T
      DO 30 I = 1, N
         J = IP(I)
         IF (J .EQ. I) GO TO 30
         IF (J .GT. 0) GO TO 10
            IP(I) = -J
            GO TO 30
 10      T = X(I)
 20      S = X(J)
         X(J) = T
         T = S
         K = J
         J = IP(K)
         IP(K) = -J
         IF (J .GT. I) GO TO 20
         X(J) = T
 30      CONTINUE
 999  RETURN
C  ***  LAST LINE OF DV7PRM FOLLOWS  ***
      END
