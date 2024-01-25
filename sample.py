#########################################################################################################
# SzekeresPy ver. 0.21 - Python package for cosmological calculations using the Szekeres Cosmological Model
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


from SzekeresPy import SzekeresModel as szekeres_cosmo  


import numpy as np
import matplotlib.pyplot as plt





#-------------------------------------------------------------
#>>>> EXAMPLE 1: checking and updating cosmological parameters 
#-------------------------------------------------------------
redshift = 0.15
H = szekeres_cosmo.H(redshift)
Om = szekeres_cosmo.Om(redshift)
Ode = szekeres_cosmo.Ode(redshift)
print("FLRW background at redshft {:.2f}, H(z) = {:.2f}, Om(z) = {:.2f}, and Ol(z) = {:.2f}".format(redshift,H,Om,Ode))

#update the parameters (only present day paramaters can be updated)
szekeres_cosmo.update(Ode0 = 0.68, Om0 = 0.32, H0 = 67.5)

#checking if the update was succesful
redshift = 0.0
H = szekeres_cosmo.H(redshift)
Om = szekeres_cosmo.Om(redshift)
Ode = szekeres_cosmo.Ode(redshift)
print("FLRW background at redshft {:.2f}, H(z) = {:.2f}, Om(z) = {:.2f}, and Ol(z) = {:.2f}".format(redshift,H,Om,Ode))



#-------------------------------------------------------------
#>>>> EXAMPLE 2: Using astropy to set the cosmological parameters
#-------------------------------------------------------------
from astropy.cosmology import Planck18  as  astropy_cosmo
import SzekeresPy
szekeres_cosmo = SzekeresPy.initiate(astropy_cosmo) 



#-------------------------------------------------------------
#>>>> EXAMPLE 3: checking the properties of the Szekeres model at a given point
#-------------------------------------------------------------
#setting the position of the observer, and checking the density, expansion, shear, weyl curvature at their locations 
t = 0 
r = 10.0 
theta = 0.25*np.pi 
phi = np.pi
redshift = 0.0
density, expansion, shear, weyl  = szekeres_cosmo.fluid(t,r,theta,phi,redshift)  
print("Relative to FLRW at redshft z= {:.2f}, density = {:.2f}, expansion rate = {:.2f}, shear = {:.2f}, Weyl = {:.2f}".format(redshift,density,expansion,shear,weyl))



#-------------------------------------------------------------
#>>>> EXAMPLE 4: plotting a 1d radial profile with R = r_max:
# density, expansion scalar, shear scalar, Weyl curvature scalar
#-------------------------------------------------------------
redshift = 0.34 
r = 50.0   # this will be (approximetly) the max R
theta = 0.75*np.pi   # direction in which you want to see the 1d profile
phi = np.pi          # direction in which you want to see the 1d profile

radius, density,expansion,shear,weyl = szekeres_cosmo.fluid_1d(t,r,theta,phi,redshift)  
plt.plot(radius,density,   radius,expansion,    radius,shear,    radius,weyl)
plt.grid()
plt.show()



#-------------------------------------------------------------
#>>>> EXAMPLE 5: #creating a 1d radial plot with R = r_max,
# and comparring density in 2 different directions
#-------------------------------------------------------------
redshift = 0.45
r = 50.0   # this will be (approximetly) the max R
theta1 = 0.05*np.pi   # direction in which you want to see the 1d profile
radius1, density1, expansion1,shear1,weyl1 = szekeres_cosmo.fluid_1d(t,r,theta1,phi,redshift)
theta2 = 0.75*np.pi   # direction in which you want to see the 1d profile
radius2, density2, expansion2,shear2,weyl2 = szekeres_cosmo.fluid_1d(t,r,theta2,phi,redshift)
plt.plot(radius1,density1,    radius2,density2)
plt.grid()
plt.show()



#-------------------------------------------------------------
#>>>> EXAMPLE 6: #plotting a 2d map of fluid properties
# in order to create a plot you need to specify a slice
# either x_slice = const, or y_slice = const, or z_slice = const,
#
# The following will result in 3 figures: 2 density and 1 shear
#
# WORK IN PROGRESS, currently you can do the following
# - only comoving domain available at this stage
# - only redshift = 0
#-------------------------------------------------------------
X_plot, Y_plot, Z_plot =  szekeres_cosmo.fluid_2d(x = 5,figure="density")
plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X_plot, Y_plot, Z_plot, cmap='cividis')
plt.show(block =False)


# we can update the dipole paramerter to make the anisotropy more prominent
# WARNING: if the dipole parameter is too large then shell crossing occurs
szekeres_cosmo.update(dipole = 0.6)
X_plot, Y_plot, Z_plot =  szekeres_cosmo.fluid_2d(y = -2,figure="density")
plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X_plot, Y_plot, Z_plot, cmap='cividis')
plt.show(block =False)


szekeres_cosmo.update(dipole = 0.1)
X_plot, Y_plot, Z_plot =  szekeres_cosmo.fluid_2d(z = 8,figure="shear")
plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X_plot, Y_plot, Z_plot, cmap='cividis')
plt.show()



#-------------------------------------------------------------
#>>>> EXAMPLE 7: #light cone
# to get a null cone you need to specify:
#  (i) the position of the observer
#  (ii) the direction
#  (iii) the redshift 
#-------------------------------------------------------------
r =50.0 
theta = 0.25*np.pi   
phi = 1.25*np.pi
RA = 12.0
DEC = -45.0 

observer = r, theta, phi
direction = RA, DEC
redshift = 4.40

#creating a 1d the path of the light ray
light_ray =  szekeres_cosmo.null_geodesic(observer,direction,redshift)  
plt.figure()
plt.plot(light_ray[1,:],light_ray[0,:]) # plot of t(R): t in Gyr, R in Mpc
plt.grid()
plt.show()
    
