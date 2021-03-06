[![Build Status](https://travis-ci.org/gmedders/ising.svg?branch=master)](https://travis-ci.org/gmedders/ising)

Ising model
===========

 - This project contains a brute-force Monte Carlo approach to obtaining thermodynamic
 properties of the Ising model, as well as the
 [cluster Monte Carlo algorithm by Wolff](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.361).

How to build:
=============

 - git clone https://github.com/gmedders/ising.git
 - cd ising
 - mkdir build
 - cd build
 - cmake ..
 - cmake --build .
 - ctest -VV

To override your default C and C++ compilers, you can use the following instead (for example):
```
cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-8 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-8 ..
```

Acknowledgments
===============
 - The cluster Monte Carlo algorithm was nicely explained in Lecture 8 of the course [Statistical Mechanics: Algorithms and Computations](https://www.coursera.org/course/smac) on Coursera. In particular, a reference implementation of the cluster Monte Carlo algorithm was provided in that course for python. Parts of it formed the basis for this c++ implementation.
