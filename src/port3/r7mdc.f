      REAL FUNCTION  R7MDC(K)
C
C  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  ***
C
      INTEGER K
C
C  ***  THE CONSTANT RETURNED DEPENDS ON K...
C
C  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS.
C  ***        K = 2... SQUARE ROOT OF ETA.
C  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH
C  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1.
C  ***        K = 4... SQUARE ROOT OF MACHEP.
C  ***        K = 5... SQUARE ROOT OF BIG (SEE K = 6).
C  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS.
C
      REAL BIG, ETA, MACHEP
C/+
      REAL SQRT
C/
      REAL R1MACH, ZERO
      EXTERNAL R1MACH
      DATA BIG/0.E+0/, ETA/0.E+0/, MACHEP/0.E+0/, ZERO/0.E+0/
      IF (BIG .GT. ZERO) GO TO 1
         BIG = R1MACH(2)
         ETA = R1MACH(1)
         MACHEP = R1MACH(4)
 1    CONTINUE
C
C-------------------------------  BODY  --------------------------------
C
      GO TO (10, 20, 30, 40, 50, 60), K
C
 10    R7MDC = ETA
      GO TO 999
C
 20    R7MDC = SQRT(256.E+0*ETA)/16.E+0
      GO TO 999
C
 30    R7MDC = MACHEP
      GO TO 999
C
 40    R7MDC = SQRT(MACHEP)
      GO TO 999
C
 50    R7MDC = SQRT(BIG/256.E+0)*16.E+0
      GO TO 999
C
 60    R7MDC = BIG
C
 999  RETURN
C  ***  LAST CARD OF  R7MDC FOLLOWS  ***
      END
