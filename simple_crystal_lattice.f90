program simple_crystal_lattice

implicit none

integer :: a
integer :: i
integer :: j
integer :: k
integer :: p
integer :: q
integer :: r
integer :: N  !number of particles
real :: s     !side lengths of crystal
real :: x
real :: y
real :: z
real :: x_cen
real :: y_cen
real :: z_cen

N = 35	!N = particles on the vertices of each cube + particles at the center of each cube
s = 1

open(unit = 2, file = "simple_crystal.xyz")

write(2,*) N
write(2,*)

!!!MAKE SURE THAT p,q,r are all one less than i,j,k

!!!particles on the vertices of each cube 
!!!# of particles = # of loops performed by i * # of loops performedd by j * # of loops performed by k
!!!e.g. i=0,2 j=0,2 k=0,2 which is  3*3*3 = 27 particles

	do i=0,2	!change the end term of each do while loop to get larger figures
		x = i*s	!for perfect cubes i=j=k, e.g. for a 2x2x2 cube let i=j=k=2

		do j=0,2
			y = j*s

			do k=0,2
				z = k*s
				write(2,*) "C", x, y, z
			end do
		
		end do	

	end do

write(2,*)

!!!particles at the center of each cube
!!!# of particles = # of loops performed by p * # of loops performedd by q * # of loops performed by r
!!!e.g. p=0,1 q=0,1 r=0,1 which is  2*2*2 = 8 particles

	do p=0,1
		x_cen = p*s + s/2

		do q=0,1
			y_cen = q*s + s/2

			do r=0,1
				z_cen = r*s + s/2
				write(2,*) "C", x_cen, y_cen, z_cen
			end do

		end do

	end do

end  program simple_crystal_lattice
