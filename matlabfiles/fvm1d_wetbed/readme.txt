FVM1D - FOR WET BED PROBLEMS ONLY
by Brett F. Sanders (and pieces of code from Scott F. Bradford)

This is a very simple 1D solver of the shallow-water equations
that uses the Hancock predictor-corrector time-stepping scheme,
the MUSCL method of slope limiting and variable reconstruction
and Roe's approximate Riemann solver to compute fluxes.

The code is kept as simple as possible to emphasize the basic
flow of logic. To account for problems involving a dry bed,
one needs to add a number of "if" statements to avoid
division by zero. This makes the code pretty messy and
therefore these lines have been omitted.

Note that the code can either be run in a first order or
second order accurate mode. The user can select from 
several limiters to see how these impact the solution.

To run this program, copy all the .m files into a directory
and run "fvm1d.m" by either typing "fvm1d" at the matlab
command prompt or pushing the execute button in the matlab
text editor.