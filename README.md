Ising model
===========

 - Currently implementing the [cluster Monte Carlo algorithm by Wolff](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.361). The acceptance
   probability for cluster expansion is valid for the 2D model with all sites
   filled. Need to determine if it is applicable more generally
   (e.g., with vacancies)

Acknowledgments
===============
 - The cluster Monte Carlo algorithm was nicely explained in Lecture 8 of the course [Statistical Mechanics: Algorithms and Computations](https://www.coursera.org/course/smac) on Coursera. In particular, a reference implementation of the cluster Monte Carlo algorithm was provided in that course for python. Parts of it formed the basis for this c++ implementation. 
