## Introduction

SzekeresPy is a Python package for cosmological calculations using the Szekeres Cosmological Model. The code is in its early development form and highly restrictive licence is imposed (educational purposes only). Once the code is more developed the licence will change to some less restrictive one.


## Requirements 

- Fortran compiler, e.g gfortran (the code was tested using GNU Fortran 11.4.0)
- Python 3.5 or higher (the code was tested using Python 3.10.12)
- numpy (the code was tested using numpy 1.21.5)
- files: 
  * szekeres_evolution.f90 
  * SzekeresPy.py 
- matplotlib
- optional: astropy


## Installation 
To compile the code open, navigate to the folder with the above files and run from the file `compile_single` or `compile_openmp` (for openMP)


## Running the code

Open the file `sample.py` and follow the sample calculations


## Development plans

- ~~ver 0.1: evolution of a model~~
- ver 0.2: light propagation 
- ver 0.3: distance-redshift relation
- ver 0.4: volume average position of the observer, radius of the domain
- ver 0.5: ray-tracing and sky maps 
- ver 0.6: change of the type of licence anticipated


## Useful references

- https://arxiv.org/abs/astro-ph/0604490
- https://arxiv.org/abs/1512.07364
- https://arxiv.org/abs/1704.02810
- unpublished notes on light tracing in the Szekeres model (contact the author for a copy)


## Questions and suggestions

In case of any questions, or suggestions for the code, please contact the author


