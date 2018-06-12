# bcc_lattice
This will be a FORTRAN 90 program to generate a bcc_lattice crystal for use in IUMD simulations. It may be folded 
The development steps will be:
1. Generate a single bcc lattice
2. Generate a tiled lattice. 
3. Add in density function to dynamically resize the spacing between crystals so that the box size is appropriately dense.
4. Generate only the appropriate number of ions from user input.
5. Generate a tiling angle.
6. Generate a crystal angle.
