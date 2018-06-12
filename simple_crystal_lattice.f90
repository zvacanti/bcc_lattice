program simple_crystal_lattice

implicit none

integer :: N
real :: x
real :: y
real :: z

x = 0.0
y = 0.0
z = 0.0
N = 9

open(unit = 2, file = "simple_crsytal.xyz")

write(2,*) N
write(2,*)
write(2,*) "C", x, y, z
write(2,*) "C", x, y, z + 1
write(2,*) "C", x, y + 1, z
write(2,*) "C", x, y + 1, z + 1
write(2,*) "C", x + 1, y, z
write(2,*) "C", x + 1, y, z + 1
write(2,*) "C", x + 1, y + 1, z
write(2,*) "C", x + 1, y + 1, z + 1
write(2,*) "C", 0.5, 0.5, 0.5

end  program simple_crystal_lattice
