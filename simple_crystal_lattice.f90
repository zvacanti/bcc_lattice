program simple_crystal_lattice

implicit none

integer :: a !a = n, the number of ions on a side
integer :: i
integer :: j
integer :: k
integer :: N  !number of particles
integer :: counter
real :: s     !side lengths of crystal
real :: rho
real :: x
real :: y
real :: z
real :: x_cen
real :: y_cen
real :: z_cen

!N = 54	!N = particles on the vertices of each cube + particles at the center of each cube

print  *, 'Enter an integer number of ions, N'
read  *, N
a = ceiling((N/2.0)**(1.0/3.0))
counter = 0
print  *, 'The number of ions is', N
print  *, 'The value of a is', a

s = 1
!print  *, 'Enter a real density value, rho'
!read  *, rho
!s = ((2)/(rho))**(1.0/3.0)
!print  *, 'The value of s is', s

open(unit = 2, file = "simple_crystal.txt")

!write(2,*) N
!write(2,*)

!!!particles on the vertices of each cube 
!!!# of particles = # of loops performed by i * # of loops performedd by j * # of loops performed by k
!!!e.g. i=0,2 j=0,2 k=0,2 which is  3*3*3 = 27 particles

	do i=0,a-1	        	!change the end term of each do while loop to get larger figures
		x = i*s + s/4		!for perfect cubes i=j=k, e.g. for a 2x2x2 cube let i=j=k=2
		x_cen = i*s + s/2 + s/4
		do j=0,a-1
			y = j*s + s/4
			y_cen = j*s + s/2 + s/4
			do k=0,a-1
				if (counter .EQ. N) EXIT
				z = k*s + s/4
				z_cen = k*s + s/2 + s/4
				write(2,*) x, y, z
				counter = counter + 1
				if (counter .EQ. N) EXIT
				write(2,*) x_cen, y_cen, z_cen
				counter = counter + 1
				if (counter .EQ. N) EXIT
			end do
		
		end do	
	
	end do

print  *, counter

end  program simple_crystal_lattice
