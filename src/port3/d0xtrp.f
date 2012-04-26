      SUBROUTINE D0XTRP(TM,M,NVAR,NG,KMAX,XPOLY,T,ERROR,EBEST,RHG,EMAG,
     1                  ESAVE)
C
      DOUBLE PRECISION TM(NVAR),NG(M),T(NVAR,KMAX),RHG(1)
C     DOUBLE PRECISION RHG(MIN(M-1,KMAX))
      REAL ERROR(NVAR,1),EBEST(NVAR),EMAG(1)
C     REAL ERROR(NVAR,MIN(M-1,KMAX)),EMAG(MIN(M-1,KMAX))
      LOGICAL XPOLY,ESAVE
C
      DOUBLE PRECISION U,V,TI,TV,TEMP
      REAL ERR
C
      IF (M.GT.1) GO TO 20
C
C ... INITIALIZE T.
C
      DO 10 I=1,NVAR
 10      T(I,1)=TM(I)
C
      GO TO 80
C
 20   MR=MIN0(M-1,KMAX)
C
      DO 30 J=1,MR
         MMJ=M-J
         RHG(J)=NG(M)/NG(MMJ)
         EMAG(J)=1.0D0+1.0D0/(RHG(J)-1.0D0)
         IF (XPOLY) RHG(J)=RHG(J)-1.0D0
 30      CONTINUE
C
      DO 70 I=1,NVAR
C
         V=0.0D0
         U=T(I,1)
         TI=TM(I)
         T(I,1)=TI
C
         DO 60 J=1,MR
C
C ......... OBTAIN SIGNED ERROR ESTIMATE.
C
            ERR=(T(I,J)-U)*EMAG(J)
            IF (ESAVE) ERROR(I,J)=ERR
            ERR=ABS(ERR)
            IF (J.EQ.1) EBEST(I)=ERR
            EBEST(I)=AMIN1(EBEST(I),ERR)
            IF (EBEST(I).EQ.ERR) JBEST=J
C
            IF (J.EQ.KMAX) GO TO 60
C
            IF (XPOLY) GO TO 40
C
C ......... RATIONAL EXTRAPOLATION.
C
            TV=TI-V
            TEMP=RHG(J)*(U-V)-TV
            IF (TEMP.NE.0.0D0) TI=TI+(TI-U)*(TV/TEMP)
            V=U
            GO TO 50
C
C ......... POLYNOMIAL EXTRAPOLATION.
C
 40         TI=TI+(TI-U)/RHG(J)
C
 50         U=T(I,J+1)
            T(I,J+1)=TI
 60         CONTINUE
C
 70      TM(I)=T(I,JBEST)
C
 80   RETURN
C
      END
