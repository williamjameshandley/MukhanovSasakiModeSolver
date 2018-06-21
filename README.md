# MukhanovSasakiModeSolver

## Set up
An example C++ project.

First check that this works on your system:

```bash
$   cd project
$   make 
$   ./bin/main
```

If the final command produces the output
```bash

Hello World
5
24
M:
1 2
3 4
x^T:
5 6
inverse of M:
  -2    1
 1.5 -0.5
Mx:
17
39
```

Then you know that it's working.

## Inspect the code

* The main driving routine is in `main.cpp`
* In this file, you can see simple usage of [Eigen](http://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html)
* You also usage of the example class defined in `src/`
* An example header/source combination is in `src/example.cpp` and `src/example.hpp`
* This defines a function and a class, which doesn't do anything particularly interesting
