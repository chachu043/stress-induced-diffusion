
 PROGRAM PGDEM3

      INTEGER PGBEG
INTEGER I,J
      REAL F(40,40),FMIN,FMAX,ALEV(1),TR(1,6)
      WRITE (*,'(A)') ' Demonstration of PGPLOT contouring routines'

      IF (PGBEG(0,'?',1,1) .NE. 1) STOP

   TR(1,1) = 0.0
      TR(1,2) = 1.0
      TR(1,3) = 0.0
      TR(1,4) = 0.0
      TR(1,5) = 0.0
      TR(1,6) = 1.0

      FMIN = 0.0
      FMAX = 0.0
      DO 20 I=1,40
          DO 10 J=1,40
              F(I,J) = COS(0.3*SQRT(I*2.)-0.4*J/3.)*COS(0.4*I/3)+(I-J)/40.0
              FMIN = MIN(F(I,J),FMIN)
              FMAX = MAX(F(I,J),FMAX)
   10     CONTINUE
   20 CONTINUE

      CALL PGPAGE
      CALL PGSVP(0.05,0.95,0.05,0.95)
      CALL PGSWIN(1.0,40.0,1.0,40.0)
      CALL PGBOX('bcts',0.0,0,'bcts',0.0,0)
      CALL PGMTXT('t',1.0,0.0,0.0,'Contouring using PGCONX with arrows')

       CALL PGBBUF
      DO 50 I=2,21,2
          ALEV(1) = FMIN + (I-1)*(FMAX-FMIN)/20.0
          WRITE (LABEL,'(I2)') I
C         WRITE (LABEL,'(F8.2)') ALEV
          IF (I.LT.10) THEN
              CALL PGSCI(2)
          ELSE
              CALL PGSCI(3)
          END IF
          CALL PGCONL(F,40,40,1,40,1,40,ALEV,TR,LABEL,16,8)
 50   CONTINUE
      CALL PGSCI(1)
      CALL PGEBUF
    
CALL PGEND
end 
