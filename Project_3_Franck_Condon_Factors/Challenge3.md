# Challenge 3

##Difference between methods for FCF calculations

Calculating FCF is equal to calculate transition probability between different modes, which is an overlap integral of wave functions of different states.
We approximate the system is to be a multi variable harmonic oscillator following Franck-Condon principle. 

### a) Hermite polynomials
The photon distribution function in this setup can be expressed by terms of multidimensional Hermite polynomials, and we are able to efficiently calculate the overlap integral by Gaussian integral.(Islampur et al. Journal of Molecular Spectroscopy Volume 194, Issue 2, April 1999, Pages 179-184)

### b) Gaussian Boson Sampling
In this method, photon scattering of specific molecule is directly sampled using photonic quantum computer or its simulator. We create a Gaussian state in Fock basis given a parameter also used in Hermite approach. This yields a population of transition between different states, therefore its spectrum with regard to wave number is proportional to the transition probability = FCF.

### c) Loop Hafnian
This method is a speedup technique of GBS. In this method, we utilize the fact that FCF of a specific setup can be calculated using hafnian of certain matrix that can be constructed using parameters for linear transformations used in GBS.


### Similarity / Differences
- Only b) directly manipulates quantum state. In other 2, the problem is reduced to direct calculatation of matrix elements.
- Doktorov operator are used in all three methods.
