PROGRAM bcc_lattice
  IMPLICIT NONE
  REAL :: a,b,c, Area
  PRINT *, 'Please enter the side lengths.'
  READ *, a, b, c
  PRINT *, 'Triangle''s Area: ', Area(a,b,c)
END PROGRAM bcc_lattice

Function Area(x,y,z)
  IMPLICIT NONE
  REAL :: Area
  REAL, INTENT ( IN ) :: x, y, z
  REAL :: theta, height
  theta = ACOS((x**2 + y**2 - z**2)/(2.0*x*y))
  height = x*SIN(theta); Area = 0.5*y*height
END FUNCTION Area
