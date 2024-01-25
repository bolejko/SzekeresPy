## Introduction

SzekeresPy is a Python package for cosmological calculations using the Szekeres Cosmological Model. The Szekeres model is an inhomogeneous, anisotropic, cosmological solution of Einstein's equations. In the limit of spherical symmetry it reduces to the Lemaitre-Tolman-Bondi (LTB) class of models, and in the limit of homogeneity, to FLRW models.
The SzekeresPy can therefore be used to investigate LTB models or to perform calculations for FLRW cosmologies, including the ΛCDM model.

The code is still in its early development. Debugging, optimisation, and restructuring of the code is still in progress. Anyone interested in collaboration or using the code is encouraged to contact the author.


## Requirements 

- Fortran compiler, e.g gfortran (the code was tested using GNU Fortran 11.4.0)
- Python 3.5 or higher (the code was tested using Python 3.10.12)
- numpy (the code was tested using numpy 1.21.5)
- files: 
  * szekeres_fortran.f90 
  * SzekeresPy.py 
- matplotlib (for plotting)
- optional: astropy (for comparison and/or importing cosmological parameters) 
- optional: healpy (for sky maps)


## Installation 

To compile the code, navigate to the folder with the above files and run the terminal commands listed either in the file `compile_single` or `compile_openmp` (for openMP). 


## Running the code

Open the file `sample.py` and follow the sample calculations and examples presented in the file:

- Example 1: checking and updating cosmological parameters 
- Example 2: using astropy to set the cosmological parameters
- Example 3: checking the properties of the Szekeres model at a given point (density, expansion rate, shear, Weyl)
- Example 4: plotting a 1d radial profile (for a specific direction) of: density, expansion rate, shear, Weyl 
- Example 5: comparing radial density profile for two different directions
- Example 6: plotting a 2d density profile for a specific slice: x=const, or y=const, or z=const; (comoving domain only)
- Example 7: plotting a path of the light ray, t(R), featuring the apparent horizon (around R ~ 3,000 Mpc - 5,000 Mpc)


## Development plans

- ~~ver 0.1: evolution of a model~~
- ~~ver 0.2: light propagation~~
- ver 0.3: distance-redshift relation
- ver 0.4: ray-tracing and sky maps 
- ver 0.5: volume average with a variable position of the observer and radius of the averaging domain


## Useful references

- https://arxiv.org/abs/astro-ph/0604490
- https://arxiv.org/abs/1412.4976
- https://arxiv.org/abs/1512.07364
- https://arxiv.org/abs/1704.02810
- unpublished notes on light tracing in the Szekeres model (contact the author for a copy)


## Questions and suggestions

In case of any questions, or suggestions for the code, please contact the author.


