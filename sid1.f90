program sid1
real, dimension (:,:), allocatable ::J,C,D,P,om,rho,s,x,y1,y2
real b,t,dt,dx,l,c1,cn,f
integer i,k,n
b=1.38*(10^-23)
print*, "discretization in x"
read*, nx
print*, "discretization in t"
read*, nt
print*, "length"
read*, L
allocate(D(1,nx))
print*,"diffusivity"
read*, D
allocate(rho(1,nx))
print*,"density"
read*, rho
allocate(s(1,nx))
print*,"specific heat"
read*, s
print*,"temperature"
read*, T

allocate(c(nt,nx))
print*,"initial concentration"
read*, c

!boundary conditions
print*, "boundary cond"
read*, c1,cn
do i=1,nt
c(i,1) = c1
c(i,nx) =cn
end do

dx= L/n
dt = 0.25* (dx**2) * minval(rho) * minval(s)
dt= dt/maxval(D)

allocate(x(1,nx))
do i=1,nx
f=f+dx
x(1,i)=f
end do

allocate(y1(1,nx))
do i=1,nx
y1(1,i) = real(c(1,i))/real(dx)
end do

allocate(y2(1,nx))
do i=1,nx
y2(1,i) = rho(1,i)*s(1,i)*dx
y2(1,i) = real(dt)/real(y2(1,i))
end do

allocate(p(nt,nx))
print*, "stress"
read*,p

allocate(J(nt,nx+1))
allocate(om(nt,nx))
print*, "omega"
read*,om


do i=1,nt-1
do m=2,nx
 J(i+1,m)= y1(1,m) * ( c(i,m-1) - c(i,m) - ((c(i,m)*om(i,m))/(b*T))*(P(i,m-1)-P(i,m)))
    J(i+1,1)= J(i+1,2)
    J(i+1,nx) = J(i+1,nx-1)
end do  

do m = 2,nx-1
    c(i+1,m) = c(i,m) + y2(1,m) * ( J(i+1,m) - J(i+1,m+1) ) 
end do
end do

end




