# bcc_lattice
This will be a FORTRAN 90 program to generate a bcc_lattice crystal for use in IUMD simulations. It may be folded into a later version of the IUMD code.
The development steps will be:
1. Generate a single bcc lattice
2. Generate a tiled lattice. 
3. Add in density function to dynamically resize the spacing between crystals so that the box size is appropriately dense.
4. Generate only the appropriate number of ions from user input.
5. Generate a tiling angle.
6. Generate a crystal angle.


The output will be a .xyz file type. When read by IUMD, it will assign the particles their velocities via a boltzmann temp. distribution. 
The .xyz file format is as follows:

1 num_of_particles
2 COMMENT LINE
3 <element> <x> <y> <z>
  ...
  
  For iumd, <element>={C,O}, with C being a neutron, O being a proton.
