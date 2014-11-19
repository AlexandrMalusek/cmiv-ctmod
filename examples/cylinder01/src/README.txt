To compile and install the code, run:

export CTMODSYS=../../..
make symlinks
make
make install


To visualize the scanner's geometry, run:
root
root [] .x root_init.C
root [] .x x_DrawGeometry.C
