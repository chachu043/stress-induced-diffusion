!======================================================================
! STEADY STATE SOLUTION IN 2-D WITH 4 TEMPERATURE BOUNDARY CONDITIONS
! First and last elements are used for boundary conditions
!======================================================================

program solve
  implicit none
  real, dimension(0:11,0:11) :: C,D,F ! D is a matrix declared to store the immediate previous values of C to calculate error
  integer :: i,j,counter,IER, PGBEG,k
  real :: e,error,ALEV(1), TR(1,6),FMIN,FMAX
  integer, dimension(200) :: X
  real,dimension(200) :: Y

! Boundary conditions
  C(0,0:11) = 100
  C(1:10,0) = 0
  C(1:10,11) =50
  C(11,1:11) =75

  counter=0 ! To check the number of steps taken

  e=10
  error = 0.00001

  do i=0,11
    do j=0,11
      D(i,j) = C(i,j)
    end do
  end do

do i =1,200
  X(i) = 0
  Y(i) = 0
end do

! NOT NECESSARY -- Initialising all the elements to 0
  do i = 1,10
    do j=1,10
      C(i,j) = 0
    end do
  end do

do while (e>error)
counter = counter +1
! Storing old values of C in D
  do i=0,11
    do j=0,11
      D(i,j) = C(i,j)
    end do
  end do

! Solving the generic equation
  do i = 1, 10
    do j = 1,10
      C(i,j) = (C(i+1,j)+C(i,j+1)+C(i-1,j)+C(i,j-1))/4.0
    end do
  end do

  e = ABS(C(5,5)-D(5,5))*100/C(5,5)
  Y(counter) = e

end do
do i=0,11
write(*,900) (C(j,i),j=0,11)
900 format(12f7.3)
end do
!print*,C(11,4)

      IF (PGBEG(0,'?',1,1) .NE. 1) STOP

      
      TR(1,1) = 0.0
      TR(1,2) = 1.0
      TR(1,3) = 0.0
      TR(1,4) = 0.0
      TR(1,5) = 0.0
      TR(1,6) = 1.0



      FMIN = 0.0
      FMAX = 0.0
      DO 20 I=0,11
          DO 10 J=0,11
              F(I,J) = C(I,J)
              FMIN = MIN(F(I,J),FMIN)
              FMAX = MAX(F(I,J),FMAX)
   10     CONTINUE
   20 CONTINUE
     


      CALL PGPAGE
      CALL PGSVP(0.0,1.0,0.0,1.0)
      CALL PGSWIN(1.0,12.0,1.0,12.0)
      CALL PGBOX('bcts',0.0,0,'bcts',0.0,0)
      CALL PGMTXT('t',1.0,0.0,0.0,'Contouring using PGCONT')
CALL PGBBUF
DO 30 I=1,50
          ALEV(1) = FMIN + (I-1)*(FMAX-FMIN)/49.0      
!IF (MOD(I,5).EQ.0) THEN
              !CALL PGSLW(3)
          !ELSE
              !CALL PGSLW(1)
          !END IF
          !IF (F(i,j).LT.50) THEN
              !CALL PGSCI(4)
              !CALL PGSLS(2)
!do k=0,11
!do j=0,11
        IF(maxval(F(:,j)).LT.50) THEN
              CALL PGSCI(4)
              CALL PGSLS(1)
       ELSE If(maxval(F(:,j)).GT.50 .and. maxval(F(:,j)).LT.75) THEN  
              CALL PGSCI(3)
              CALL PGSLS(1)
ELSE If(maxval(F(:,j)).GT.75 .and. maxval(F(:,j)).LT.100) THEN  
              CALL PGSCI(2)
              CALL PGSLS(1)

          END IF

      CALL PGBBUF
      
          
          CALL PGCONT(F,12,12,1,12,1,12,ALEV,-1,TR)
!end do 
!end do
   30 CONTINUE
       CALL PGSLW(1)
      CALL PGSLS(1)
      CALL PGSCI(1)
      CALL PGEBUF
call PGEND

end program
