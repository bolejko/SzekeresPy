#########################################################################################################
# SzekeresPy ver. 0.18 - Python package for cosmological calculations using the Szekeres Cosmological Model
# 
# File: sample.py
# 
# Author: Krzysztof Bolejko
# 
# Intended use: research and education
# 
# Licence: BSD-2-Clause ("FreeBSD License")
# 
# Copyright (c) 2024 Krzysztof Bolejko
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#########################################################################################################


import numpy as np
import matplotlib.pyplot as plt


from   SzekeresPy  import SzekeresModel  as  szekeres_cosmo  # ALTERNATIVELY:  import SzekeresPy as szpy  ->  szekeres_cosmo = szpy.initiate()  
#from astropy.cosmology import Planck18  as  astropy_cosmo


#check the parameters of the background FLRW model
redshift = 0.15
H = szekeres_cosmo.H(redshift)
Om = szekeres_cosmo.Om(redshift)
Ode = szekeres_cosmo.Ode(redshift)
print("FLRW background at redshft {:.2f}, H(z) = {:.2f}, Om(z) = {:.2f}, and Ol(z) = {:.2f}".format(redshift,H,Om,Ode))

#setting the position of the observer, and checking the density, expansion, shear, weyl curvature at their locations 
t = 0 
r = 10.0 
theta = 0.85*np.pi 
phi = np.pi
redshift = 0.0
density, expansion, shear, weyl  = szekeres_cosmo.fluid(t,r,theta,phi,redshift)  
print("Relative to FLRW at redshft {:.2f}, density = {:.2f}, expansion rate = {:.2f}, shear = {:.2f}, Weyl = {:.2f}".format(redshift,density,expansion,shear,weyl))


#creating a 1d radial plot with r_max = r, in the direction of constant theta and phi, at a given redshift
redshift = 0.1
r = 50.0
radius, density,expansion,shear,weyl = szekeres_cosmo.fluid_1d(t,r,theta,phi,redshift)  
plt.figure()
plt.plot(radius,density,radius,expansion,radius,shear,radius,weyl)
plt.grid()
plt.show(block=False)


#creating a 1d radial plot of the path of the light ray
print("work in PROGRESS - so not results in Figure 2 yet")
light_ray =  szekeres_cosmo.null_geodesic(t,r,theta,phi,redshift)  
plt.figure()
plt.plot(light_ray[1,:],light_ray[0,:]) # plot of t(r) 
plt.grid()
plt.show(block=False)






