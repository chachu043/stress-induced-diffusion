program composite1
real, dimension (:,:), allocatable ::r,k,rho,c,T,x,xJ,y1,y2,J,temp400sec,temp200sec,temp80sec,temp40sec
real, dimension (:,:), allocatable ::temp_400sec,temp_200sec,temp_80sec,temp_40sec

real  L, r1,c1,c2,k1,k2,rho1,rho2,dx,dt,s, area, pi
integer IER, PGBEG, nx, nt,i,m

pi= 3.141592  
nx = 150
nt= 120000
L=0.2
r1=0.01
allocate(r(1,nx))
r=r1
!print*, r

!thermal conductivity (400 W/(m.degC)
    k1 = 400
    k2 = 50
allocate(k(1,nx))
    k = k1
do i=nx/2,nx
    k(1,i) = k2
end do
    !k(100:end) = k3
! specific heat capacity (380 J/(kg.degC)
    c1 = 380  
    c2 = 450
   allocate(c(1,nx))
    c= c1
do i=nx/2,nx
    c(1,i) = c2
end do
    !k(100:end) = k3
! density  (8900 kg/m^3)
    rho1 = 8900
    rho2 = 7900
    allocate(rho(1,nx))
    rho = rho1
do i=nx/2,nx
    rho(1,i) = rho2
    !k(100:end) = k3
end do

allocate(T(nt,nx))
T=0.0

!boundary conditions
do i=1,nt
T(i,1) = 100.0
T(i,nx) = 0.0
end do

dx= L/nx
dt = 0.25* (dx**2) * minval(rho) * minval(c)
dt= dt/maxval(k)
print*, dt
allocate(x(1,nx))
x= dx/2+ dx

!print*, x(1,150)

allocate(xJ(1,nx))
do i=1,nx
s=s+dx
xJ(1,i)= s

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
 J(i+1,m)= y1(1,m) * ( T(i,m-1) - T(i,m) );
    J(i+1,1)= J(i+1,2);
    J(i+1,nx) = J(i+1,nx-1);
end do  

do m = 2,nx-1
    T(i+1,m) = T(i,m) + y2(1,m) * ( J(i+1,m) - J(i+1,m+1) ) ;
end do
end do

! energy flux dQ/dt = JA
    dQ_dt = J(nt,nx) * area;

! temperature gradient dT/dx
    dT_dx = (T(nt,nx) - T(nt,1)) / (xJ(1,nx) - xJ(1,1));
allocate(temp400sec(1,nx))
allocate(temp200sec(1,nx))
allocate(temp80sec(1,nx))
allocate(temp40sec(1,nx))
allocate(temp_400sec(1,nx))
allocate(temp_200sec(1,nx))
allocate(temp_80sec(1,nx))
allocate(temp_40sec(1,nx))

do i=1,nx
temp400sec(1,i) = J(119921,i)
end do
do i=1,nx
temp200sec(1,i) = J(59960,i)
end do
!print*, temp200sec
do i=1,nx
temp80sec(1,i) = J(23984,i)
end do
do i=1,nx
temp40sec(1,i) = J(11992,i)
end do
do i=1,nx
temp_400sec(1,i) = T(119921,i)
end do
do i=1,nx
temp_200sec(1,i) = T(59960,i)
end do
!print*, temp200sec
do i=1,nx
temp_80sec(1,i) = T(23984,i)
end do
do i=1,nx
temp_40sec(1,i) = T(11992,i)
end do

!print*, dT_dx
IER = PGBEG(0,'?',1,1)
IF (IER.NE.1) STOP
CALL PGSVP(0.0,0.5,0.5,1.0)
!CALL PGENV(0.0,0.2,0.0,300000.0,0,1) 
!CALL PGLAB('position (m)','FLux density (W/m^2)','Flux density vs position graph')
!CALL PGSCI(2)
!CALL PGLINE(150,xJ,temp400sec)
!CALL PGTEXT(0.16,290000.0,'flux at 400sec')
!CALL PGSCI(3)
!CALL PGLINE(150,xJ,temp200sec)
!CALL PGTEXT(0.16,280000.0,'flux at 200sec')
!CALL PGSCI(4)
!CALL PGLINE(150,xJ,temp80sec)
!CALL PGTEXT(0.16,270000.0,'flux at 80sec')
!CALL PGSCI(5)
!CALL PGLINE(150,xJ,temp40sec)
!CALL PGTEXT(0.16,260000.0,'flux at 40sec')
!CALL PGSVP(0.5,1.0,0.0,1.0)
CALL PGENV(0.0,0.2,0.0,100.0,0,1) 
CALL PGLAB('position (m)','temperature (k)','temperature vs position graph')
CALL PGSCI(2)
CALL PGLINE(200,xJ,temp_400sec)
CALL PGTEXT(0.17,95.0,'temp at 400sec')
CALL PGSCI(3)
CALL PGLINE(200,xJ,temp_200sec)
CALL PGTEXT(0.17,91.0,'temp at 200sec')
CALL PGSCI(4)
CALL PGLINE(200,xJ,temp_80sec)
CALL PGTEXT(0.17,87.0,'temp at 80sec')
CALL PGSCI(5)
CALL PGLINE(200,xJ,temp_40sec)
CALL PGTEXT(0.17,83.0,'temp at 40sec')
CALL PGEND
    






end
