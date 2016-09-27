
! Author: chanakya yadav mm14b043
! using an implicit (Backward Euler method) for the radioactive decay of a nuclear rod for 0 < r < 100cm, 0 < t < 100years
! solving T_t = kappa*(T_rr + T_r/r) + S(r,t)

program cylinder
implicit none
real, dimension (:,:), allocatable :: A,B,C,T_tk,temp,ti,temp_sol,T_tk1
real, dimension (:,:), allocatable ::source_vector,source,temp1year,temp10year,temp50year,temp100year, r, D
real  dr, dt, kappa, s, isol,r1,r2,rc,a_source,Trod,Tau_0,t1,T,p,omega,eps
integer n, m,i, j,k,x,y, IER, PGBEG,iter
omega=1.2
eps=0.0001
r1=0
rc=100
r2=rc
T=100. !length of time
n=99
m=1000.
dr= real(rc)/real(n+1) !discretization in r spacial direction.
dt= real(T)/real(m) !discretization in time
t1=0.0              !initial time    
kappa= 2*(10**7)    ! thermal diffusivity of ground
s= real(kappa*dt)/real(dr**2) 
!print *, kappa,s, dt


! build the A matrix to march finite difference solution forward in time
! check first with n=3 to test that you are getting the right A matrix and b vector:
allocate(A(n,n))
allocate(B(n,1))
B=0.0
do i=1,n
do j=1,n
if (i .eq. j) then
A(i,j) = 2*s
A(i,j) = A(i,j)+ 1.0
else if (i .eq. j+1) then
A(i,j) = -s +real(s)/real(2*i)
else if (i .eq. j-1) then
A(i,j) = -s - real(s)/real(2*i)
else 
A(i,j) =0.0
end if
end do
end do

! specify boundary conditions through vector D and by changing any rows in A matrix needed - for Neumann boundary conditions

! for Neumann boundary conditions at r=0: dT(0,t)/dr = 0 is approximated by T_0^k = T_1^k
allocate(D(n,1))
D=0
A(1,1) = s+ (real(s)/real(2))
A(1,1) = A(1,1)+ real(1.0000) !neumann T_0=T_1, 1+2s-s-s/2
!for Dirichlet boundary conditions at r=rc: b[n] = (-s-s/(2*n))*300. 
!T(r=rc, t) = 300.
D(n,1)= -(-s - real(s)/real(2*n))*300.
print *,A(1,1), A(2,2), A(2,1), A(1,2), D(n,1)

allocate(T_tk1(n,1))
T_tk1=300.
!print *,D
allocate(temp1year(n,1))
temp1year=0
allocate(temp10year(n,1))
temp10year=0
allocate(temp50year(n,1))
temp50year=0
allocate(temp100year(n,1))
temp100year=0

!a mesh in r direction
r1= r1+dr
!print*, r1
allocate(r(n,1))
do y=1,n
r(y,1) = r1
r1=r1+dr
end do
!print*, r
!Set up vector T_tk: T(r,0) = 300K 
allocate(T_tk(n,1))
T_tk = 300.9
!print*, T_tk

!store solution for each time in temp
allocate(temp(n,m))
do i=1,n
temp(i,1) = T_tk(i,1)
end do
allocate(ti(m,1))
ti=t1 !time vector

!define source term due to radioactive decay of nuclear rod
!source vector is only nonzero for r < a=25cm
!print *,D
allocate(source_vector(n,1))
source_vector = 0
allocate(source(n,1))
source= 0

a_source = 25.   !a = 25cm radius of rod is 25cm
Trod = 1.  !initial temperature change due to nuclear rod is 1K
tau_0 = 100.  !half-life of rod is 100 years

do i = 1,n
    if (r(i,1) .lt. a_source) then
       source_vector(i,1) = 1.0
    end if
end do
!write(*,10) (source_vector(i,1),i=1,n)
!10 format(5f4.2)

isol = 1
allocate(temp_sol(4,1))
temp_sol = 0 

!now march solution forward in time
do k= 1,m
t1= t1 + dt
p=real(Trod*exp(-t1/tau_0))/real(a_source**2)
do i=1,n
source(i,1) = p*source_vector(i,1)
source(i,1) = kappa*dt*source(i,1)
end do
do i=1,n
B(i,1)= T_tk(i,1)+ D(i,1)+ source(i,1)
end do
!call gauss(A,B,T_tk1,n)
call gs_sor(A,B,T_tk1,omega,eps,n,iter)

!print*,T_tk1
!if (k.eq.10)then
!print*,T_tk1
!end if
do i=1,n
temp(i,k) = T_tk1(i,1)
end do
!for next time step
do i=1,n
T_tk(i,1) = T_tk1(i,1)
end do
T_tk1=300.0
ti(k,1) = t1

if (k.eq.10) then
       temp_sol(1,1) = t1
       temp1year = T_tk
elseif (k.eq.100) then
       temp_sol(2,1) = t1
       
       temp10year = T_tk
elseif (k.eq.500) then
       temp_sol(3,1) = t1
      
       temp50year = T_tk
elseif (k.eq.1000) then
       temp_sol(4,1) = t1
       
       temp100year = T_tk
end if
end do
!print*, T_tk
print*, temp1year
print*, temp10year
print*, temp50year
print*, temp100year
!print*, temp_sol

IER = PGBEG(0,'?',1,1)
IF (IER.NE.1) STOP
CALL PGENV(0.,100.,300.,301.,0,1) 
CALL PGLAB('r (cms)','temperature (k)','temperature distribution near nuclear rod at different time intervals')
CALL PGSCI(2)
CALL PGLINE(n,r,temp1year)
CALL PGTEXT(80.0,300.9,'temp at 1years')
CALL PGSCI(3)
CALL PGLINE(n,r,temp10year)
CALL PGTEXT(80.0,300.8,'temp at 10years')
CALL PGSCI(4)
CALL PGLINE(n,r,temp50year)
CALL PGTEXT(80.0,300.7,'temp at 50years')
CALL PGSCI(5)
CALL PGLINE(n,r,temp100year)
CALL PGTEXT(80.0,300.6,'temp at 100years')
CALL PGEND
end


 subroutine gs_sor(a,b,x,omega,eps,n,iter)
implicit none 
integer, parameter::kmax=5000
integer n
real a(n,n), b(n), x(n)
real c, omega, eps, delta, conv, sum
integer i, j, k, iter, flag

! check if the system is diagonally dominant
flag = 0
do i=1,n
  sum = 0.0
  do j=1,n
    if(i == j) cycle
    sum = sum+abs(a(i,j))
  end do
  if(abs(a(i,i)) < sum) flag = flag+1
end do
if(flag >0) write(*,*) 'The system is NOT diagonally dominant'    

do k=1,kmax
  conv = 0.0
  do i=1,n
    delta = b(i)
    do j=1,n
      delta = delta - a(i,j)*x(j)
    end do
    x(i) = x(i)+omega*delta/a(i,i)
    if(abs(delta) > conv) conv=abs(delta)
  end do
  if(conv < eps) exit
end do
iter = k
if(k == kmax) write (*,*)'The system failed to converge'

end subroutine gs_sor


