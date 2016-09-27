program try2
!two blocks each of concentration 0.75 (0<x<50, 0<y<100)  and 0.25 (0<x<50, 0<y<100) and diffusivity 10^(-9). 
! stability criteria for 2d is 0.25=Fo. ergo, when Dx= 10^3 Dt<=250
real A(10,10), B(10,10), Dx, Dt
integer i,j, t 
Dx= 0.001
Dt = 250
do i=1,10
do j=1,5
A(i,j)=0.75
end do
end do
do i=1,10
do j=6,10
A(i,j)=0.25
end do
end do 
do i=1,10
write (*,11) (A(i,j),j=1,10)
11 format(10f5.2)
end do
do t= 1,100000
B=A
do i=2,9
do j=2,9
A(i,j) = 0.25*(B(i-1,j)+ B(i+1,j) + B(i,j-1) + B(i,j+1))
end do
end do
do j=2,9
A(1,j)=0.5*B(2,j)+ 0.25*(B(1,j+1)+B(1,j-1))
A(10,j)=0.5*B(9,j)+ 0.25*(B(10,j+1)+B(10,j-1))
end do
end do
do i=1,10
write (*,10) (A(i,j),j=1,10)
10 format(10f6.3)
end do
end
