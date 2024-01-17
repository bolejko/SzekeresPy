##########################################################################################################
# SzekeresPy ver. 0.1 - Python package for cosmological calculations using the Szekeres Cosmological Model
# 
# File: SzekeresPy.py 
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

import szekeres_evolution as fortran_evolution


class Szekeres:

    def __init__(self, cospar,szpar):
        self.cospar = cospar
        self.szpar = szpar
        
    def All(self): #returns all values      
        return self.cospar,self.szpar
    def E(self,z):
        Ez = np.sqrt( self.cospar[1]*((1+z)**3) + self.cospar[2] )
        return Ez
    def H(self,z):   #returns H(z)
        Hz = self.cospar[0]*self.E(z)
        return Hz
    def Om(self,z): #returns omega_matter(z)
        Omz =   ((self.cospar[1])*((1+z)**3))/(self.E(z)**2)
        return Omz   
    def Ode(self,z):  #returns omega_lambda(z)
        Odez =   (self.cospar[2])/(self.E(z)**2)
        return Odez  #returns omega_lambda
	
    def fluid(self,t,r,theta,phi,redshift,fluid_flag = None):
        if fluid_flag == None:
            fluid_flag = 'density'
        if fluid_flag == 'density':
            fluid_flag = 0
        if fluid_flag == 'expansion':
            fluid_flag = 1
        if fluid_flag == 'shear':
            fluid_flag = 2
        if fluid_flag == 'weyl':
            fluid_flag = 3
        if fluid_flag == 'ricci':
            fluid_flag = 4
        point = np.zeros(7)
        point = t,r,theta,phi,redshift,1,1
        npoint = len(point)
        ncospar = len(self.cospar)
        ninhomog = len(self.szpar)

        Ngrid = 1   
        nv = Ngrid*np.ones(1)
        aa = np.append(self.cospar,self.szpar)
        bb = np.append(point,nv)
        input_data = np.append(aa,bb)
        rho,tht,shr,wey,ric,arl,prp = fortran_evolution.link_point(input_data)
        return rho,tht,shr,wey
                       
    def fluid_1d(self,t, r, theta, phi, redshift, fluid_flag = None, radius_flag = None, Ngrid = None):

        point = np.zeros(7)
        point = t,r,theta,phi,redshift,1,1
        npoint = len(point)
        ncospar = len(self.cospar)
        ninhomog = len(self.szpar)
        fluid = np.zeros(10) 

        if fluid_flag == None:
            fluid_flag = 'density'
        if radius_flag == None:
            radius_flag = 'areal'
        if fluid_flag == 'density':
            fluid_flag = 0
        if fluid_flag == 'expansion':
            fluid_flag = 1
        if fluid_flag == 'shear':
            fluid_flag = 2
        if fluid_flag == 'weyl':
            fluid_flag = 3
        if fluid_flag == 'ricci':
            fluid_flag = 4
        if radius_flag == 'areal':
            radius_flag  = 6  
        if radius_flag == 'proper':
            radius_flag  = 7
        if Ngrid == None:
            Ngrid = 100      
        if Ngrid < 1 or Ngrid > 1000:
           print("Number of points should be between 1 and 1000, Ngrid = 1000 selected")
           Ngrid = 1000
        Ngrid = 100    
        nv = Ngrid*np.ones(1)
        aa = np.append(self.cospar,self.szpar)
        bb = np.append(point,nv)
        input_data = np.append(aa,bb)
        rho,tht,shr,wey,ric,arl,prp = fortran_evolution.link_multi(input_data)
        radius = prp

        return prp,rho,tht,shr,wey

def initiate(astropy_cosmo=None, inhomog_cosmo=None):

  cospar = np.zeros(15)
  szpar = np.zeros(15)

  background_dict = {
      "H0"   : 70.0,
      "Om0"  :  0.3,
      "Ode0" : 0.7
  }    

  perturbation_dict = {
      "contrast" : -0.0025,    # contrast at the CMB, not today!
      "radius"   :  10.0,
      "slope"    :  0.4,
      "dipole"   :  0.4
  }   


  if astropy_cosmo == None:
      cospar[0] = background_dict["H0"]
      cospar[1] = background_dict["Om0"]
      cospar[2] = background_dict["Ode0"]
  else:
      if "<class 'astropy.cosmology." in str(type(astropy_cosmo)):
          cospar[0] = astropy_cosmo.H0.value
          cospar[1] = astropy_cosmo.Om0
          cospar[2] = astropy_cosmo.Ode0
      else: 
          cospar[0] = background_dict["H0"]
          cospar[1] = background_dict["Om0"]
          cospar[2] = background_dict["Ode0"]


  if inhomog_cosmo == None:
      szpar[0] = perturbation_dict["contrast"]
      szpar[1] = perturbation_dict["radius"]
      szpar[2] = perturbation_dict["slope"]
      szpar[3] = perturbation_dict["dipole"]      
  else:
      szpar[0] = perturbation_dict["contrast"]
      szpar[1] = perturbation_dict["radius"]
      szpar[2] = perturbation_dict["slope"]
      szpar[3] = perturbation_dict["dipole"]  


  sz_cosmo = Szekeres(cospar,szpar)
  

  return sz_cosmo         

SzekeresModel = initiate()  


