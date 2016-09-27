program composite1
real, dimension (:,:), allocatable ::r,k,rho,c,T,x,xJ,y1,y2,J,temp400sec,temp200sec,temp80sec,temp40sec

real  L, r1,c1,c2,k1,k2,rho1,rho2,dx,dt,s, area, pi,k3,rho3,c3
integer IER, PGBEG, nx, nt,i,m

pi= 3.141592  
nx = 150
nt= 130000
L=0.25
r1=0.01
allocate(r(1,nx))
r=r1
!print*, r


k1 = 0.72
k2 = 0.034
k3 =1.33
allocate(k(1,nx))
k = k1
do i=70,90
k(1,i) = k2
end do
do i=91,150
k(1,i) = k3
end do

    
c1 = 1000  
c2 = 1000
c3=1000
allocate(c(1,nx))
c= c1

  
rho1 = 1
rho2 = 1
rho3=  1
allocate(rho(1,nx))
rho = rho1
allocate(T(nt,nx))
T=0.0

!boundary conditions
do i=1,nt
T(i,1) = 1200.0
T(i,nx) = 50.0
end do

dx= L/nx
dt = 0.25* (dx**2) * minval(rho) * minval(c)
dt= dt/maxval(k)
allocate(x(1,nx))
x= dx/2+ dx

!print*, x(1,150)

allocate(xJ(1,nx))
do i=1,nx
xJ(1,i)= s
s=s+dx
end do
!print*, xJ(1,2)

!allocate(t(1,nt-1))
!t= dt

allocate(y1(1,nx))
do i=1,nx
y1(1,i) = real(k(1,i))/real(dx)
end do

allocate(y2(1,nx))
do i=1,nx
y2(1,i) = rho(1,i)*c(1,i)*dx
y2(1,i) = real(dt)/real(y2(1,i))
end do

area= pi*(r(1,1)**2)
allocate(J(nt,nx+1))

do i=1,nt-1
do m=2,nx
 J(i+1,m)= y1(1,m) * ( T(i,m-1) - T(i,m) )
    J(i+1,1)= J(i+1,2)
    J(i+1,nx) = J(i+1,nx-1)
end do  

do m = 2,nx-1
    T(i+1,m) = T(i,m) + y2(1,m) * ( J(i+1,m) - J(i+1,m+1) ) 
end do
end do

! energy flux dQ/dt = JA
    dQ_dt = J(nt,nx) * area

! temperature gradient dT/dx
    dT_dx = (T(nt,nx) - T(nt,1)) / (x(1,nx) - x(1,1))
allocate(temp400sec(1,nx))
allocate(temp200sec(1,nx))
allocate(temp80sec(1,nx))
allocate(temp40sec(1,nx))

do i=1,nx
temp400sec(1,i) = T(130000,i)
end do
print*, temp400sec


IER = PGBEG(0,'?',1,1)
IF (IER.NE.1) STOP
CALL PGENV(0.0,0.25,0.0,1200.0,0,1) 
CALL PGLAB('temperature', 'position', 'final temperature distribution')
CALL PGLINE(150,xJ,temp400sec)
!CALL PGLINE(150,xJ,temp200sec)
!CALL PGLINE(150,xJ,temp80sec)
!CALL PGLINE(150,xJ,temp40sec)
CALL PGEND
    






end
