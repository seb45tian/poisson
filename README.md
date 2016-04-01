## poisson

This program will perform an iterative calculation of the 2D poisson equation using Jacobis's method. It outputs two data files potential.dat and field.dat.
Two python scripts are provided for visualising the potential and field (plot.py) and for calculating and plotting an analytical approximation.

Usage
------
A Makefile is provided using the mpicc compiler which outputs a poisson binary.

The python scripts require the matplotlib and numpy library.

    $ make
    $ ./poisson
    $ python plot.py

Author
------
Sebastian Potthoff: <s.potthoff@warwick.ac.uk>

University of Warwick, April 2016
