program simple_crystal_lattice
implicit none

integer :: a, i, j, k, N, counter 					!a = n, the number of ions on a side
									!N = total num ions
real :: s, length, rho							!side lengths of crystal
real :: x, y, z								!Ion locations
real :: x_cen, y_cen, z_cen						!Center ions
real :: x_copy, y_copy, z_copy						!Needed for loops
real :: theta, phi, psi							!Rotation angles
real, dimension (3,1) :: xyz_prime

real :: x1, x2, y1, y2, z1, z2, distance
real, allocatable :: positions(:,:)					!Positions that get written.
real, allocatable :: posit_rot(:,:)					!Positions BEFORE rotating. 
real, allocatable :: distances(:,:)					!Distances that get written.
real, allocatable ::  dist_rot(:,:)					!Distances BEFORE rotating.

N = 16
print  *, 'Enter an integer number of ions, N'
read  *, N

allocate( positions(4, N) )
allocate( distances(N, N) )
allocate( posit_rot(4, N) )
allocate(  dist_rot(N, N) )

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

!s = 1
print  *, 'Enter a real density value, rho'
read  *, rho
s = (2.0/rho)**(1.0/3.0)
print  *, 'The value of s is', s

length = s*a								!is the total side length. This is relevant later. 
print *, "Side length is: ", length

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
			counter = counter +1
			z = k*s + s/4
			z_cen = k*s + s/2 + s/4
			x_copy = x
			y_copy = y
			z_copy = z
			posit_rot(1,counter) = counter
			posit_rot(2,counter) = x_copy
			posit_rot(3,counter) = y_copy
			posit_rot(4,counter) = z_copy
			call rotation(x_copy,y_copy,z_copy,theta,phi,psi)
			call boundary_conditions(x_copy, y_copy, z_copy, length)
			write(2,*) x_copy, y_copy, z_copy
			positions(1,counter) = counter
			positions(2,counter) = x_copy
			positions(3,counter) = y_copy
			positions(4,counter) = z_copy
			if (counter .EQ. N) EXIT

!-----------------------Below, we take care of the central ions--------------------------------
			counter = counter + 1
			x_copy = x_cen
			y_copy = y_cen
			z_copy = z_cen
			posit_rot(1,counter) = counter
			posit_rot(2,counter) = x_copy
			posit_rot(3,counter) = y_copy
			posit_rot(4,counter) = z_copy
			call rotation(x_copy,y_copy,z_copy,theta,phi,psi)
			call boundary_conditions(x_copy, y_copy, z_copy, length)
			write(2,*) x_copy, y_copy, z_copy
			positions(1,counter) = counter
			positions(2,counter) = x_copy
			positions(3,counter) = y_copy
			positions(4,counter) = z_copy
			if (counter .EQ. N) EXIT
		end do
	end do	
end do
	close( 2, status='keep')

i=1
j=1
distance = 0
counter = 0
!print *, ""
!print *, "           ", "   i ", "   j ", "  Distance  "
!print *, "--------------------------------------------------------"

!-------------------------------This line checks distances-------------------------
do i=1, N
	do j=i+1, N
		x1		= positions(2, i)
		x2 		= positions(2, j)
		y1 		= positions(3, i)
		y2		= positions(3, j)
		z1 		= positions(4, i)
		z2	 	= positions(4, j)
		call euclid(distance,length, x1, y1, z1, x2, y2, z2)
		distances(i,j)  = distance

		x1		= posit_rot(2, i)
		x2 		= posit_rot(2, j)
		y1 		= posit_rot(3, i)
		y2		= posit_rot(3, j)
		z1 		= posit_rot(4, i)
		z2	 	= posit_rot(4, j)
		call euclid(distance,length, x1, y1, z1, x2, y2, z2)
		dist_rot(i,j)  = distance
		
!		if (abs(distances(i,j) - dist_rot(i,j)) > 0.01*length) then		!Prints differences for debugging.
!			print 2000, i, j, abs(distances(i,j) - dist_rot(i,j))
!	 2000 format ( '   ERROR   ', i5, i5, f10.5)
!		else
!			print 2002, i, j, abs(distances(i,j) - dist_rot(i,j))
!	 2002 format ( '           ', i5, i5, f10.5)
!		endif 
	enddo
enddo

deallocate (positions, distances)
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
	x = s+mod(x,s)
else if (x > s) then
	x = mod(x,s)
else
	x = x
endif

if (y < 0) then
	y = s+mod(y,s)
else if (y > s) then
	y = mod(y,s)
else
	y = y
endif

if (z < 0) then
	z = s+mod(z,s)
else if (z > s) then
	z = mod(z,s)
else
	z = z
endif

end subroutine boundary_conditions

!------------------------------------------------------------------------------------
!This is a testing subroutine.
!These are pulled from wikipedia. 
subroutine euclid(distance,length, x1, y1, z1, x2, y2, z2)
implicit none
real, intent( inout)	:: distance, x1, y1, z1
real, intent( in ) 	:: length, x2, y2, z2

real :: dx, dy, dz

x1 = x2 - floor(x1/length) * length				!Alternate position equation.
dx = x2 - x1
dx = dx - nint(dx / length) * length

y1 = y2 - floor(y1/length) * length
dy = y2 - y1
dy = dy - nint(dy / length) * length

z1 = z2 - floor(z1/length) * length
dz = z2 - z1
dz = dz - nint(dz / length) * length

distance = (dx**2 + dy**2 + dz**2)**(1.0/2.0)

end subroutine euclid


