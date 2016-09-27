program ADI
real, dimension (:,:), allocatable ::A,B,T_tk,T_tk1,T_tk2,D, bound1,bound2,E,F,G,P,Q,X,y,H
real dt, dx, T,l,s,kappa,ti,omega,eps
integer m,n,i,j,k,u,iter
omega=1.2
eps=0.00001
m=100
n=3
T=0.0 !k
l=40.0 !cms
dt=10
ti=0
print*,dt
dx=real(l)/real(n+1)
print*,dx**2
kappa= 0.835 !cm^2/sec
s= (kappa*dt)/(dx**2) !stability criterion says lamda<=0.25 for 2d problems
print*,s

allocate(A(n,n))
A=0.
allocate(B(n,n))
B=0.
allocate(G(n,n))
G=0.
allocate(G(n,n))
G=0.

do i=1,n
do j=1,n
if (i .eq. j) then
A(i,j) = 2*(1.+ s)
else if (i .eq. j+1) then
A(i,j) = -s
else if (i .eq. j-1) then
A(i,j) = -s
end if
end do
end do
do i=1,n
write(*,10) (A(i,j),j=1,n)
10 format(3f7.3)
end do


do i=1,n
do j=1,n
if (i .eq. j) then
B(i,j) = 2*(1.-s)
else if (i .eq. j+1) then
B(i,j) = s
else if (i .eq. j-1) then
B(i,j) = s
end if
end do
end do
do i=1,n
!write(*,11) (B(i,j),j=1,n)
!11 format(3f7.3)
end do
allocate(X(n,1))
allocate(E(1,n))
E=0.
allocate(F(n,1))
F=0.
allocate(D(1,n))
D=0.
print*,D
allocate(P(n,1))
P=0.
allocate(y(n,1))
y=0.
allocate(Q(n,1))
Q=0.
allocate(T_tk(n,n))
T_tk=T
do i=1,n
!write(*,18) (T_tk(i,j),j=1,n)
!18 format(3f7.3)
end do
allocate(T_tk1(n,n))
T_tk1=T
allocate(T_tk2(n,n))
T_tk2=T

allocate(bound1(n,n))
bound1=0
do i=1,n
j=1
bound1(i,j)= (75.*s)
end do
do i=1,n
j=n
bound1(i,j)= (50.*s)
end do
do i=1,n
!write(*,12) (bound1(i,j),j=1,n)
!12 format(3f7.3)
end do
allocate(bound2(n,n))
bound2=0
do j=1,n
i=1
bound2(i,j)= (0.*s)
end do
do j=1,n
i=n
bound2(i,j)= (100.*s)
end do
do i=1,n
!write(*,13) (bound2(i,j),j=1,n)
!13 format(3f8.3)
end do
bound1= bound1+bound2
do i=1,n
write(*,20) (bound1(i,j),j=1,n)
20 format(3f8.3)
end do
do i=1,n
do j=1,n
bound2(i,j) = bound1(j,i)
end do
end do

!do k= 1,10000
ti=ti+dt

print*,D
do i=1,n
do j=1,n
X(j,1)= bound1(j,i)
end do
!print*,X
!print*,'--------------'
do j=1,n
E(1,j) = B(i,j)
end do
D = matmul(E,T_tk)
do j=1,n
P(j,1)=D(1,j)
end do
P = P+X
call gs_sor(A,P,Q,omega,eps,n,iter)
do j=1,n
T_tk(i,j) = Q(j,1)
T_tk1(i,j) = Q(j,1)
end do
!print*, T_tk
!print*, '--------------'
end do
!print*, T_tk1
do i=1,n
write(*,311) (T_tk1(i,j),j=1,n)
311 format(3f8.3)
end do


do i=1,n
do j=1,n
X(j,1) = bound2(j,i)
end do
!print *,X
print*,'-------'
do j=1,n
E(1,j)=B(i,j)
end do
do j=1,n
F(j,1)=B(j,i)
end do
!print*, F
!print*, P
D= matmul(T_tk1,F)
!print*,D
!print*,D(2,1)
do j=1,n
y(j,1) = D(1,j)
end do
y(2,1)=D(1,2)
!print*,D(1,2)
!print*,'-------'
D=D+X
print*,D
!print*,'-----------'
call gs_sor(A,D,Q,omega,eps,n,iter)
do j=1,n
T_tk1(j,i)=Q(j,1)
T_tk2(i,j)=Q(j,1)
end do
end do
do i=1,n
write(*,312) (T_tk2(i,j),j=1,n)
312 format(3f8.3)
end do


end program

 subroutine gs_sor(a,b,x,omega,eps,n,iter)
!==========================================================
! Solutions to a system of linear equations A*x=b
! Method: The successive-over-relaxation (SOR)
! Alex G. (November 2009)
!----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - array of the right hand coefficients b
! x(n)   - solutions (initial guess)
! n      - number of equations (size of matrix A)
! omega  - the over-ralaxation factor
! eps    - convergence tolerance 
! output ...
! x(n)   - solutions
! iter   - number of iterations to achieve the tolerance
! coments ...
! kmax   - max number of allowed iterations
!==========================================================
implicit none 
integer, parameter::kmax=10000
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

  subroutine gauss_1(a,b,x,n)
!============================================================
! Solutions to a system of linear equations A*x=b
! Method: the basic elimination (simple Gauss elimination)
! Alex G. November 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - vector of the right hand coefficients b
! n      - number of equations
! output ...
! x(n)   - solutions
! comments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real a(n,n), b(n), x(n)
real c
integer i, j, k

!step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do

!step 2: back substitution
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*x(j)
   end do 
   x(i) = (b(i)- c)/a(i,i)
end do
end subroutine gauss_1
