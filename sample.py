##########################################################################################################
# SzekeresPy ver. 0.11 - Python package for cosmological calculations using the Szekeres Cosmological Model
# 
# File: sample.py 
# 
# Author: Krzysztof Bolejko
# 
# Licence to use: restrictive licence
# 
# Intended use: educational purposes only
# 
# Copyright (c) 2024 Krzysztof Bolejko
#
# The Author grants a licence, free of charge, to any other person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software with restrictions, including limitation of the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, subject to the following conditions:
#
# The above copyright notice and this permission notice must be included in all
# copies or substantial portions of the Software.
#
# The intended use of this Software is for educational purposes only,
# and it is explicitly prohibited to use, copy, modify, merge, publish, 
# distribute, sublicense, and/or sell copies of the Software for any purposes
# other than educational purposes. 
# 
# Any other use of this Software, or any dealing in this Software other than
# as described in this notice is prohibited. The use of this Software by entities
# other than humans is prohibited.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##########################################################################################################

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
theta = 0.5*np.pi 
phi = np.pi
redshift = 0.0
density, expansion, shear, weyl  = szekeres_cosmo.fluid(t,r,theta,phi,redshift)  
print("Relative to FLRW at redshft {:.2f}, density = {:.2f}, expansion rate = {:.2f}, shear = {:.2f}, Weyl = {:.2f}".format(redshift,density,expansion,shear,weyl))


#creating a 1d radial plot with r_max = r, in the direction of constant theta and phi, at a given redshift
redshift = 0.5
r = 50.0
radius, density,expansion,shear,weyl = szekeres_cosmo.fluid_1d(t,r,theta,phi,redshift)  
plt.plot(radius,density,radius,expansion,radius,shear,radius,weyl)
plt.grid()
plt.show()









