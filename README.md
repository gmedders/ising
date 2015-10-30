Ising model
===========

 - Currently implementing the [cluster Monte Carlo algorithm by Wolff](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.361). The acceptance
   probability for cluster expansion is valid for the 2D model with all sites
   filled. Need to determine if it is applicable more generally
   (e.g., with vacancies)

Acknowledgments
===============
 - The cluster Monte Carlo algorithm was nicely explained in Lecture 8 of the course [Statistical Mechanics: Algorithms and Computations](https://www.coursera.org/course/smac) on Coursera. In particular, a reference implementation of the cluster Monte Carlo algorithm was provided in that course for python. Parts of it formed the basis for this c++ implementation. 
 - The `random_selector` struct that enables an object to be selected at random from a c++ container comes from a [gist](https://gist.github.com/5538174.git), which was provided as part of the response to [this question](http://stackoverflow.com/a/16421677).
