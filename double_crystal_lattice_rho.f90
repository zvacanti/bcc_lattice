program simple_crystal_lattice

implicit none

integer :: a !a = n, the number of ions on a side in crystal A
!integer :: b !the number of ions on a side in crystal B
integer :: i
integer :: j
integer :: k
!integer :: p
!integer :: q
!integer :: r
integer :: N  !number of particles
integer :: counter1
integer :: counter2
real :: s	!distance between each ion in crystal A
!real :: l	!distance between each ion in crystal B
real :: rho
real :: x1
real :: y1
real :: z1
real :: x1_cen
real :: y1_cen
real :: z1_cen
real :: x2
real :: y2
real :: z2
real :: x2_cen
real :: y2_cen
real :: z2_cen
real :: d	!entire side length of crystal A

!N = 54	!N = particles on the vertices of each cube + particles at the center of each cube

print  *, 'Enter an integer number of ions, N'
read  *, N
a = ceiling((N/2.0)**(1.0/3.0))
counter1 = 0
counter2 = 0
print  *, 'The number of ions per lattice is', N
print  *, 'The number of ions per side of crystal A (a) is', a

!s = 1
print  *, 'Enter a real density value, rho'
read  *, rho
s = (2.0/rho)**(1.0/3.0)
print  *, 'The distance between each ion in crystal A (s) is', s

open(unit = 2, file = "double_crystal.txt")

d = s*a

!write(2,*) N
!write(2,*)

!!!particles on the vertices of each cube 
!!!# of particles = # of loops performed by i * # of loops performedd by j * # of loops performed by k
!!!e.g. i=0,2 j=0,2 k=0,2 which is  3*3*3 = 27 particles

	do i=0,a-1	        	!change the end term of each do while loop to get larger figures
		x1 = i*s + s/4		!for perfect cubes i=j=k, e.g. for a 2x2x2 cube let i=j=k=2
		x1_cen = i*s + s/2 + s/4
		do j=0,a-1
			y1 = j*s + s/4
			y1_cen = j*s + s/2 + s/4
			do k=0,a-1
				if (counter1 .EQ. N) EXIT
				z1 = k*s + s/4
				z1_cen = k*s + s/2 + s/4
				if (counter2 .EQ. N) EXIT
				z2 = z1 + d
				z2_cen = z1_cen + d
				write(2,*) x1, y1, z1
				write(2,*) x1, y1, z2
				counter1 = counter1 + 1
				counter2 = counter2 + 1
				if (counter1 .EQ. N) EXIT
				if (counter2 .EQ. N) EXIT
				write(2,*) x1_cen, y1_cen, z1_cen
				write(2,*) x1_cen, y1_cen, z2_cen
				counter1 = counter1 + 1
				counter2 = counter2 + 1
				if (counter1 .EQ. N) EXIT
				if (counter2 .EQ. N) EXIT
			end do
		
		end do	
	
	end do

print  *, counter1
print  *, counter2

!!Creates the second crystal

	!do p=0,b-1
		!x = p*l + l/4
		!x_cen = p*l + l/2 + l/4
		!do q=0,b-1
			!y = q*l + l/4
			!y_cen = q*l + l/2 + l/4			
			!do r=0,b-1
				!if (counter .EQ. N) EXIT
				!z = k*l + l/4 + 
				!z_cen = k*l + l/2 + l/4 + 
				!write(2,*) x, y, z
				!counter = counter + 1
				!if (counter .EQ. N) EXIT
				!write(2,*) x_cen, y_cen, z_cen
				!counter = counter + 1
				!if (counter .EQ. N) EXIT
			!end do			

		!end do
		
	!end do 

end  program simple_crystal_lattice
