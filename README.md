# ABOUT

This project provides a sample program, that provides a significant
numerical computation, written using a number of languages and
approaches for the sake of comparison.  The sample program is a
slightly modified version of a Navier-Stokes fluid dynamics code by
[Griebel et al.][1].

Current implementations include C, Fortran, and Haskell (repa). 

# History 
2012/13 - Accelerate (Haskell) and Ypnos (Haskell) implementations in progress.

# Building and Running
## C

```
cd c/nast2d-dom/
make
./navier
^C
```

^C when you're done running.

## FORTRAN
```
cd fortran
make
./navier
^C
```
^C when you're done running. Then make the video by:

```
./build-vid output
```

# References 
[1]: Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer, Numerical Simulation in Fluid Dynamics, SIAM, 1998. http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html
