# RBF-PUM
Code sample for the Radial Basis Function Partition of Unity Method.

By Genesis J. Islas

The work found here is a select set of numerical experiments completed for my dissertation in support of my PhD degree in applied mathematics at Arizona State University.

The sample code found here computes the numerical solution of a set of reaction diffusion equations on the surface of the Stanford bunny using the radial basis function partition of unity method.
Details of this method can be found in Chapters 3 and 4 of my dissertation.

##########################################################

Dependencies:

The following libraries are needed:

HQRRP  (https://github.com/flame/hqrrp.git)

This version of this package used in these experiments has been updated by Jason Yalim and Nicholas Chmielweski and can be found in the folder "HQRRP-main".

##########################################################

Instructions:

Run rbf_pum_bunny.m on MATLAB

The code will save MATLAB figures of following:
1. Interpolation Nodes
2. Eigenvalues of the Laplacian Matrix
3. Numerical solution of the final time
