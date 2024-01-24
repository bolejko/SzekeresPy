#########################################################################################################
# SzekeresPy ver. 0.20 - Python package for cosmological calculations using the Szekeres Cosmological Model
# 
# File: SzekeresPy.py
# 
# Author: Krzysztof Bolejko
# 
# Intended use: research and education
# 
# Licence to use: BSD-2-Clause ("FreeBSD License")
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


import szekeres_fortran as fortran

import numpy as np
from scipy.interpolate import LinearNDInterpolator


class Szekeres:

    def __init__(self, cospar,szpar,cube_register):
        self.cospar = cospar
        self.szpar = szpar
        self.cube_register = cube_register
        self.box = 4
        self.side = self.box*2+1
        self.grid = (self.side)**3
        self.cube_data = np.zeros((4,self.grid))

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
	

    def update(self,H0 = None, Om0 = None, Ode0 = None):
        if H0 == None:
            print("H0 not updated")
        else:
            self.cospar[0] = H0 
        if Om0 == None:
            print("Om0 not updated")
        else:
            self.cospar[1] = Om0
        if Ode0 == None:
            print("Ode0 not updated")
        else:
            self.cospar[2] = Ode0
        

    def cube_done(self,z):
        if z < 0:
            zi = 0
        elif z>9.99:
            zi = 9999
        else:
            zi = int(np.floor(z*100))    
        self.cube_register[zi] = 1.0 


    def cube_check(self,z):
        if z < 0:
            zi = 0
        elif z>9.99:
            zi = 9999
        else:
            zi = int(np.floor(z*100))
        a = self.cube_register[zi] 
        if a:
            print("CUBE exists and ready to use:",a)
        return a
            

    def fluid(self,t,r,theta,phi,redshift):
        point = np.zeros(7)
        point = t,r,theta,phi,redshift,1,1
        Ngrid = 1   
        input_data = Ngrid*np.ones(1)
        input_data = np.append(input_data,self.cospar)
        input_data = np.append(input_data,self.szpar)
        input_data = np.append(input_data,point)
        rho,tht,shr,wey,ric,arl,prp = fortran.link_point(input_data)
        return rho,tht,shr,wey
                       

    def fluid_1d(self, t, r, theta, phi, redshift,radius_flag = None):
        point = np.zeros(7)
        point = t,r,theta,phi,redshift,1,1
        Ngrid = 100    
        input_data = Ngrid*np.ones(1)
        input_data = np.append(input_data,self.cospar)
        input_data = np.append(input_data,self.szpar)
        input_data = np.append(input_data,point)
        rho,tht,shr,wey,ric,arl,prp = fortran.link_multi(input_data)
        if radius_flag == None:
            radius_flag = 'areal'
        if radius_flag == 'proper':
            radius = prp
        elif radius_flag == 'areal':
            radius = arl
        else:
            radius = arl 
        return radius,rho,tht,shr,wey


    def fluid_2d(self, x=None, figure=None):
        if figure == None:
            figure = "density"
        else:
            figure = "density"
        if x==None:
            x_slice = 0.0
        else:
            x_slice = 0.0
        redshift = 0.0
        point = np.zeros(7)
        point = 0.0,30.0,1,1,redshift,1,1
        Ngrid = 100    
        input_data = Ngrid*np.ones(1)
        input_data = np.append(input_data,self.cospar)
        input_data = np.append(input_data,self.szpar)
        input_data = np.append(input_data,point)

        X_plot = np.outer(np.linspace(-1, 1, self.side), np.ones(self.side))
        Y_plot = X_plot.copy().T
        Z_plot = np.outer(np.linspace(-1, 1, self.side), np.ones(self.side))

        r_range = 4.0*self.szpar[1]
        X_plot = X_plot*r_range
        Y_plot = Y_plot*r_range 

        if not self.cube_check(redshift):
            rho,tht,shr,wey,ric,com,prp = fortran.link_cube(input_data)
            self.cube_data[0] = rho
            self.cube_data[1] = tht
            self.cube_data[2] = shr
            self.cube_data[3] = wey
            self.cube_done(redshift)
        
        ni = 1
        nj = self.side
        nk = self.side

        for i in range(ni):
            for j in range(nj):
                for k in range(nk):

                    pix=self.pixelTO(i,j,k)
                    values = self.cube_data[0][pix]
                    Z_plot[i][j] = values

        
        return X_plot, Y_plot, Z_plot
        

    def pixelFROM(self,i):
        pix = i
        z = int(pix / (self.side*self.side))
        pix = pix - (z * self.side * self.side)
        y = int(pix / self.side)
        x = int(pix % self.side)
        return i,j,k
   
    def pixelTO(self,i,j,k):
        xshift = i   #+self.box
        yshift = j   #+self.box
        zshift = k   #+self.box
        pix = xshift + self.side*(yshift + self.side*zshift)
        return pix
   

    def fluid_3d(self, t, r, theta, phi, redshift,grid_flag = None):
        point = np.zeros(7)
        point = t,r,theta,phi,redshift,1,1
        Ngrid = 100    
        input_data = Ngrid*np.ones(1)
        input_data = np.append(input_data,self.cospar)
        input_data = np.append(input_data,self.szpar)
        input_data = np.append(input_data,point)
        rho,tht,shr,wey,ric,com,prp = fortran.link_cube(input_data)       
        if grid_flag == None:
            grid_flag = 'comoving'
        if grid_flag == 'proper':
            grid = prp
        elif grid_flag == 'comoving':
            grid = com
        else:
            grid = com
        return grid,rho,tht,shr,wey
        
        
        
    def null_geodesic(self,obs,dir,red):
        point = np.zeros(7)
        point = 0.0,obs[0],obs[1],obs[2],red,1,1
        direction = np.zeros(7)
        direction = dir[0],dir[1],1,1,1,1,1
        Ngrid = 200   
        input_data = Ngrid*np.ones(1)
        input_data = np.append(input_data,self.cospar)
        input_data = np.append(input_data,self.szpar)
        input_data = np.append(input_data,point)
        input_data = np.append(input_data,direction)
        light_ray = np.zeros((6,Ngrid)) 
        tempral_position, radial_position, theta_position, phi_position,  testal , redshift_position = fortran.link_null(input_data)
        light_ray[0,:] = tempral_position
        light_ray[1,:] = radial_position
        light_ray[2,:] = theta_position    
        light_ray[3,:] = phi_position
        light_ray[4,:] = testal
        light_ray[5,:] = redshift_position      
        return light_ray
               

               
def initiate(astropy_cosmo=None, inhomog_cosmo=None):
  cospar = np.zeros(15)
  szpar = np.zeros(15)
  cube_register = np.zeros(10000)

  background_dict = {
      "H0"   : 68.8,
      "Om0"  :  0.3,
      "Ode0" : 0.7
  }    

  perturbation_dict = {
      "contrast" : -0.002,    
      "radius"   :  15.0,
      "slope"    :  0.4,
      "dipole"   :  0.4
  }   


  show_logo = False

# FIX needed: assigning values and updating parameters 

  if astropy_cosmo == None:
      cospar[0] = background_dict["H0"]
      cospar[1] = background_dict["Om0"]
      cospar[2] = background_dict["Ode0"]
  else:
      if "<class 'astropy.cosmology." in str(type(astropy_cosmo)):
          cospar[0] = astropy_cosmo.H0.value
          cospar[1] = astropy_cosmo.Om0
          cospar[2] = astropy_cosmo.Ode0
          show_logo = False
          background_dict["H0"] = astropy_cosmo.H0.value
          background_dict["Om0"] = astropy_cosmo.Om0
          background_dict["Ode0"]= astropy_cosmo.Ode0
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



  logo = '''    
   _____            __                       ______  __
  / ___/____  ___  / /_____  ________  _____/ __ \ \/ /
  \__ \/_  / / _ \/ //_/ _ \/ ___/ _ \/ ___/ /_/ /\  / 
 ___/ / / /_/  __/ ,< /  __/ /  /  __(__  ) ____/ / /  
/____/ /___/\___/_/|_|\___/_/   \___/____/_/     /_/   

  '''                                                       

  if show_logo:
    print(logo)
  print("    ")  
  print("INITIATED WITH FOLLOWING PARAMETERS")
  print("FLRW cosmology",background_dict )
  print("Szekeres parameters",perturbation_dict )
  print("    ")  
# FIX needed: optimise for MCMC

  sz_cosmo = Szekeres(cospar,szpar,cube_register)
  
  return sz_cosmo         

SzekeresModel = initiate()  
