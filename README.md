Minim
=====


This is a simple energy minimization program using steepest descents
and a trivially simple forcefield to minimize the energy of
water. This is done in internal coordinates, so there are only 2
parameters to optimize:

- The O-H distance
- The H-O-H angle

`OriginalEmin.py` is J.H. Jensen's Python code from his GIST at https://gist.github.com/jhjensen2/6060440

`Emin.c` is my version, rewritten in C as the basis of something more useful. WHile this version is still in internal coordinates, it uses a structure to hold the coordinates and another to store the energy and the gradients. It also separates the code into Energy/Gradient calculation (as in the original), plus a separate functions for performing the SD minimization. Variable names have been changved to be more descriptive.

Compile the code with:
```
cc -o emin EMin.c -lm
```
