program bcc_lattice
implicit none

integer, parameter :: dbl=kind(1.0d0)

integer :: a, i, j, k, N, counter								!a = n, the number of ions on a side
																			!N = total num ions
real :: s, length, rho												!crystal side lengths
real :: x, y, z
real :: x_cen, y_cen, z_cen
real :: x_copy, y_copy, z_copy
real :: angx, angy, angz, angle
real, dimension (3,1) :: xyz_prime

real :: x1, x2, y1, y2, z1, z2, distance
real(dbl), allocatable ::   pos_pre_rot(:,:)
real(dbl), allocatable ::  pos_post_rot(:,:)
real, allocatable ::  dist_pre_rot(:,:)
real, allocatable :: dist_post_rot(:,:)
real, allocatable ::   ang_pre_rot(:,:)
real, allocatable ::  ang_post_rot(:,:)
real, dimension (2,2) :: testA
real, dimension (2,1) :: testB
real, dimension (2,1) :: testC

double precision pi
parameter (pi = 3.1415926535897932D0)

double precision rad
parameter (rad = 0.174533)

N = 51200
rho = 0.05	!units are in fermi
s = (2.0/rho)**(1.0/3.0)																	!Test param

allocate(   pos_pre_rot(3,N) )
allocate(  pos_post_rot(3,N) )

a = ceiling((N/2.0)**(1.0/3.0))									!N=2*a^3
counter = 0

print *, "Ions per side (a): ", a
print *, "The distance between each ion(s) is: ", s
print *, "x angle in radians:"
!read *, angx
print *, "y angle in radians:"
!read *, angy
print *, "z angle in radians:"
!read *, angz

angx = 0
angy = 0
angz = rad*4.5

!s        = 1.0

length   = s*a															!length = interion spacing * num ions
distance = 0
angle		= 0

open(unit = 2, file = "bcc.xyz")

!This makes the actual lattice. Re-write version.

z = 0
do i=0, a-1
	x 		= i*s + s/4.0
	x_cen = i*s + s/4.0 + s/2.0
	do j=0, a-1
		y 		= j*s + s/4.0
		y_cen = j*s + s/4.0 + s/2.0
		do k=0, a-1
		z 		= k*s + s/4.0
		z_cen = k*s + s/4.0 + s/2.0
		x_copy = x														!Need to make copies, or we screw
		y_copy = y														!up these values when we need them
		z_copy = z														!for the central ions later.
		if (counter .eq. N) exit
		counter = counter + 1
		pos_pre_rot(1,counter) = x_copy
		pos_pre_rot(2,counter) = y_copy
		pos_pre_rot(3,counter) = z_copy
!		print *, "pre  (m): ", counter, x_copy, y_copy, z_copy
		call rotation (x_copy, y_copy, z_copy, angx, angy, angz)
		call boundary_conditions(x_copy, y_copy, z_copy, length)
		write(2,*) x_copy, y_copy, z_copy
		pos_post_rot(1,counter) = x_copy
		pos_post_rot(2,counter) = y_copy
		pos_post_rot(3,counter) = z_copy
!		print *, "post (m): ", counter, x_copy, y_copy, z_copy
!		print *, ""
		
		if (counter .eq. N) exit
		counter = counter +1
		x_copy = x_cen
		y_copy = y_cen
		z_copy = z_cen
		pos_pre_rot(1,counter) = x_copy
		pos_pre_rot(2,counter) = y_copy
		pos_pre_rot(3,counter) = z_copy
!		print *, "pre  (c): ", counter, x_copy, y_copy, z_copy
		call rotation (x_copy, y_copy, z_copy, angx, angy, angz)
		call boundary_conditions(x_copy, y_copy, z_copy, length)
		write(2,*) x_copy, y_copy, z_copy
		pos_post_rot(1,counter) = x_copy
		pos_post_rot(2,counter) = y_copy
		pos_post_rot(3,counter) = z_copy
!		print *, "post (c): ", counter, x_copy, y_copy, z_copy
!		print *, ""
		if (counter .eq. N) exit
		enddo
	enddo
enddo
	close( 2, status = "keep")

call x8b_format(N, rho, pos_post_rot)

!allocate(  dist_pre_rot(N,N) )
!allocate( dist_post_rot(N,N) )
!allocate(   ang_pre_rot(N,N) )
!allocate(  ang_post_rot(N,N) )

!print *, "          i", "           j", "   distance"
!print *, "----------------------------------------------------"
!do i=1, N
!	do j=i+1, N
!		x1 				= pos_pre_rot(1,i)
!		x2 				= pos_pre_rot(1,j)
!		y1 				= pos_pre_rot(2,i)
!		y2 				= pos_pre_rot(2,j)
!		z1 				= pos_pre_rot(3,i)
!		z2 				= pos_pre_rot(3,j)
!		call euclid(distance, angle, length, x1, y1, z1, x2, y2, z2)
!		dist_pre_rot(i,j)	= distance
!		 ang_pre_rot(i,j) = angle
!
!		x1 				= pos_post_rot(1,i)
!		x2 				= pos_post_rot(1,j)
!		y1 				= pos_post_rot(2,i)
!		y2 				= pos_post_rot(2,j)
!		z1 				= pos_post_rot(3,i)
!		z2 				= pos_post_rot(3,j)
!		call euclid(distance, angle, length, x1, y1, z1, x2, y2, z2)
!		dist_post_rot(i,j)	= distance
!		 ang_post_rot(i,j) 	= angle
!		print *, i, j, dist_pre_rot(i,j)-dist_post_rot(i,j)
!		if (abs(dist_pre_rot(i,j) - dist_post_rot(i,j))>0.1*length .or.&
!			abs(ang_pre_rot(i,j) - ang_post_rot(i,j))>0.00174533) then
!			print *, "ERROR", i, j, abs(dist_pre_rot(i,j) - dist_post_rot(i,j)), &
!				abs(ang_pre_rot(i,j)-ang_post_rot(i,j))
!		end if
!	enddo
!enddo

!deallocate (dist_pre_rot, dist_post_rot(N,N), ang_pre_rot, ang_post_rot )

end program bcc_lattice

!----------------------------------------------------------

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

!----------------------------------------------------------

subroutine euclid(distance, angle, length, x1, y1, z1, x2, y2, z2)
implicit none
real, intent( inout ) :: distance, angle, x1, y1, z1
real, intent( in ) 	 :: length, x2, y2, z2

real 						 :: dx, dy, dz, lena, lenb, dot

x1 = x2 - floor(x1/length) * length								!Alternate position equation.
dx = x2 - x1
dx = dx - nint(dx / length) * length

y1 = y2 - floor(y1/length) * length
dy = y2 - y1
dy = dy - nint(dy / length) * length

z1 = z2 - floor(z1/length) * length
dz = z2 - z1
dz = dz - nint(dz / length) * length

distance = (dx**2 + dy**2 + dz**2)**(1.0/2.0)

dot = x1*x2 + y1*y2 + z1*z2
lena= sqrt(x1**2 + y1**2 + z1**2)
lenb= sqrt(x2**2 + y2**2 + z2**2)

angle = acos( (dot)/( lena*lenb))

end subroutine euclid

!----------------------------------------------------------

subroutine rotation (x, y, z, angx, angy, angz)
implicit none
real, intent ( inout ) :: x, y, z, angx, angy, angz
real, dimension (3,3)  :: rz, rx, ry, rzy, rzyx
real, dimension (3,1)  :: mx
real, dimension (3,1)  :: mx_prime

mx(1,1) 	= x
mx(2,1) 	= y
mx(3,1) 	= z

rz(1,1)	= cos(-angz)												!angz is the appropriate angle:
rz(1,2)	= sin(-angz)												!draw it out and check it. 
rz(2,1)	=-sin(-angz)												!which axis does this actually
rz(2,2)	= cos(-angz)												!describe a rotation around?
rz(3,1)  = 0
rz(3,2)  = 0
rz(1,3)  = 0
rz(2,3)  = 0
rz(3,3)  = 1

ry(1,1) 	= cos(-angy)
ry(1,2)	= 0
ry(1,3)	= sin(-angy)
ry(2,1)	= 0
ry(2,2)	= 1
ry(2,3)	= 0
ry(3,1)	=-sin(-angy)
ry(3,2)	= 0
ry(3,3)	= cos(-angy)

rx(1,1)	= 1
rx(1,2)	= 0
rx(1,3)	= 0
rx(2,1)	= 0
rx(2,2)	= cos(-angx)
rx(2,3)	= sin(-angx)
rx(3,1)	= 0
rx(3,2)	=-sin(-angx)
rx(3,3)	= cos(-angx)

rzy	 	= matmul(rz, ry)
rzyx 		= matmul(rzy, rx)
mx_prime = matmul(rzyx, mx)

!print *, rzyx(1,1), rzyx(1,2), rzyx(1,3)
!print *, rzyx(2,1), rzyx(2,2), rzyx(2,3)
!print *, rzyx(3,1), rzyx(3,2), rzyx(3,3)

x 			= mx_prime(1,1)
y			= mx_prime(2,1)
z			= mx_prime(3,1)


end subroutine rotation

!----------------------------------------------------------

subroutine x8b_format(N, rho, positions)

implicit none

!!!Header
integer, parameter :: dbl=kind(1.0d0)
integer, intent (in) :: N
real, intent (in) :: rho
real(dbl), intent (in) :: positions(3,N)
integer :: i
character(10)	code_name		!name of program which created file
character(8)	code_version		!program version
character(8)	xdate			!date configuration was written to
					!file, from Fortran date_and_time
character(10)	xdaytime		!time of day, from date_and_time
character(5)	xtimezone		!time zone, from date_and_time
character(20)	ion			!simulation type
real(dbl)	time			!simulation time stamp(REQUIRED)
real(dbl)	aspect(3)		!aspect ratio of box edge lengths(REQUIRED)
real(dbl)	ev			!average potential energy per particle
real(dbl)	ek			!average kinetic energy per particle
real(dbl)	px			!pressure
real(dbl)	pp(3,3)			!pressure tensor

!!!Particle data

character(6)	xftype			!file type
logical		xappend			!.true. = append to existing file
					!.false. = write new file
logical		xclose			!.true. = close file after writing
					!.false. = leave file open
logical		xopnd			!.true. unit 14 is already opened

time = 0
aspect(1) = 1
aspect(2) = 1
aspect(3) = 1
ev = 0
ek = 0
px = 0
pp = 0

inquire(14,opened=xopnd)
open(14,FILE="md.00000000000.x8b",STATUS='UNKNOWN',FORM='UNFORMATTED')

call date_and_time(xdate,xdaytime,xtimezone)
write(14) xftype
write(14) code_name,code_version
write(14) xdate,xdaytime,xtimezone
write(14) "ion"
write(14) time, rho, aspect, ev, ek, px, pp, N

write(14) (positions(1,i), positions(2,i), positions(3,i), i=1, N)

 close(14)

end subroutine x8b_format

!------------------------------------------------------------------

!This makes lines. For testing. Inset it into line 63 to impliment.

!do i=1, N/2
!	pos_pre_rot(1,i) = i
!	pos_pre_rot(2,i) = 0
!	pos_pre_rot(3,i) = i
!	x = pos_pre_rot(1,i)
!	y = pos_pre_rot(2,i)
!	z = pos_pre_rot(3,i)
!	call rotation(x,y,z,angx,angy,angz)
!	pos_post_rot(1,i) = x
!	pos_post_rot(2,i) = y
!	pos_post_rot(3,i) = z
!	write(2, *) x, y, z
!	print *, i, x, y, z
	
!	pos_pre_rot(1,i+N/2) = 0
!	pos_pre_rot(2,i+N/2) = i
!	pos_pre_rot(3,i+N/2) = i+N/2
!	x = pos_pre_rot(1,i+N/2)
!	y = pos_pre_rot(2,i+N/2)
!	z = pos_pre_rot(3,i+N/2)
!	call rotation(x,y,z,angx,angy,angz)
!	pos_post_rot(1,i+N/2) = x
!	pos_post_rot(2,i+N/2) = y
!	pos_post_rot(3,i+N/2) = z
!	write(2,*) x, y, z
!	print *, i, x, y, z
!enddo
