program simple_crystal_lattice
!Further work: Make it so we can specify the box size.

implicit none

integer :: a 															!a = n, the number of ions on a side
integer :: i
integer :: j
integer :: k
integer :: N  															!number of particles
integer :: counter
real :: s, length														!side lengths of crystal
real :: rho
real :: x, y, z														!Ion locations
real :: x_cen, y_cen, z_cen										!Center ions
real :: x_copy, y_copy, z_copy									!Needed for loops
real :: theta, phi, psi												!Rotation angles
real, dimension (3,1) :: xyz_prime

N = 54	!N = particles on the vertices of each cube + particles at the center of each cube

!print  *, 'Enter an integer number of ions, N'
!read  *, N
a = ceiling((N/2.0)**(1.0/3.0))
counter = 0
print  *, 'The number of ions is', N
print  *, 'The number of ions per side (a) is', a

print *, 'Enter a rotation angle around x in radians'
read *, theta
print *, 'Enter a rotation angle around y in radians'
read *, phi
print *, 'Enter a rotation angle around z in radians'
read *, psi

s = 1
!print  *, 'Enter a real density value, rho'
!read  *, rho
!s = (2.0/rho)**(1.0/3.0)
!print  *, 'The value of s is', s

length = s*a !is the total side length. This is relevant later. 

open(unit = 2, file = "simple_crystal.txt")

!write(2,*) N
!write(2,*)

!!!particles on the vertices of each cube 
!!!# of particles = # of loops performed by i * # of loops performedd by j * # of loops performed by k
!!!e.g. i=0,2 j=0,2 k=0,2 which is  3*3*3 = 27 particles

!Note that I make copies of everything; If we didn't, it would throw the loops off completely. 
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
				x_copy = x
				y_copy = y
				z_copy = z
				call rotation(x_copy,y_copy,z_copy,theta,phi,psi)
				!call boundary_conditions(x_copy, y_copy, z_copy, length)
				write(2,*) x_copy, y_copy, z_copy
				counter = counter + 1
				if (counter .EQ. N) EXIT
				x_copy = x_cen
				y_copy = y_cen
				z_copy = z_cen
				call rotation(x_copy,y_copy,z_copy,theta,phi,psi)
				!call boundary_conditions(x_copy, y_copy, z_copy, length)
				write(2,*) x_copy, y_copy, z_copy
				counter = counter + 1
				if (counter .EQ. N) EXIT
			end do
		
		end do	
	
	end do

print  *, counter
print *, x, y, z
call rotation(x,y,z,theta, phi,psi)
!call boundary_conditions(x,y,z,length)
print *, x_copy, y_copy, z_copy
end  program simple_crystal_lattice

!---------------------------------------------------------------------------------------
!This makes the rotation matrix and rotates each particle. 
subroutine rotation (x, y, z, theta, phi, psi)
implicit none
real, intent ( inout ) :: x, y, z, theta, phi, psi

real, dimension (3,3) :: R
real, dimension (3,1) :: MX
real, dimension (3,1) :: MX_PRIME

MX(1,1) = x
MX(2,1) = y
MX(3,1) = z

!Need to use radians, because that's the output of cos in FORTRAN.
!This subroutine is off. The crystal breaks up, for some reason.
R(1,1) =  cos(phi)*cos(psi) - sin(phi)*cos(theta)*sin(phi)
R(1,2) = -cos(phi)*sin(psi) - sin(phi)*cos(theta)*cos(phi)
R(1,3) =  cos(phi)*sin(theta)

R(2,1) =  sin(phi)*cos(psi) + cos(phi)*cos(theta)*sin(phi)
R(2,2) = -sin(phi)*sin(psi) + cos(phi)*cos(theta)*cos(phi)
R(2,3) = -cos(phi)*sin(theta)

R(3,1) = sin(theta)*sin(psi)
R(3,2) = sin(theta)*cos(psi)
R(3,3) = cos(theta)

MX_PRIME = MATMUL(R, MX)

x = MX_PRIME(1,1)
y = MX_PRIME(2,1)
z = MX_PRIME(3,1)

end subroutine rotation


!----------------------------------------------------------------------------
!This makes it so we don't violate the periodic boundary conditions. 
!It's simple, because we can't possibly rotate it past -s, 2s. So 
!we don't have to do much math really. 

subroutine boundary_conditions(x, y, z, s)
implicit none
real, intent( inout) :: x, y, z, s

if (x < 0) then
	x = abs(x)
else if (x > s) then
	x = x - s
else
	x = x
endif

if (y < 0) then
	y = abs(y)
else if (y > s) then
	y = y - s
else
	y = y
endif

if (z < 0) then
	z = abs(z)
else if (z > s) then
	z = z - s
else
	z = z
endif

end subroutine boundary_conditions





















