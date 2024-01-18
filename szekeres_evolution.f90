!########################################################################################################
! SzekeresPy ver 0.11 - Python package for cosmological calculations using the Szekeres Cosmological Model
! 
! File: szekeres_evolution.f90
! 
! Author: Krzysztof Bolejko
! 
! Licence to use: restrictive licence
! 
! Intended use: educational purposes only
! 
! Copyright (c) 2024 Krzysztof Bolejko
!
! The Author grants a licence, free of charge, to any other person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software with restrictions, including limitation of the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, subject to the following conditions:
!
! The above copyright notice and this permission notice must be included in all
! copies or substantial portions of the Software.
!
! The intended use of this Software is for educational purposes only,
! and it is explicitly prohibited to use, copy, modify, merge, publish, 
! distribute, sublicense, and/or sell copies of the Software for any purposes
! other than educational purposes. 
! 
! Any other use of this Software, or any dealing in this Software other than
! as described in this notice is prohibited. The use of this Software by entities
! other than humans is prohibited.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!########################################################################################################


program szekeres_evolution
implicit none
    integer :: npypac,npyszek,I
    integer, parameter :: npoint = 7
    double precision, allocatable :: pypac(:)
    double precision, allocatable :: pyszek(:)    
    double precision, dimension (npoint) :: szpoint  
    double precision, dimension (100) :: szpac  
    character(len=8), dimension(100) :: szpan

    double precision,  dimension (10) :: fluid
    double precision :: redshift,time
    double precision, parameter :: pi = 4.0d0*datan(1.0d0)

    
    npypac = 15
    allocate(pypac(npypac))
    npyszek = 15
    allocate(pyszek(npyszek))



! overide of model parameters [if you wish to use the code as a standalone fortran program with no dependency on python fornt-end]
   szpac = 0.0d0


   call parameter_values(npypac,pypac, npyszek,pyszek, szpac)

   call parameter_names(szpan)

   szpoint(1) = 0.0 
   szpoint(2) = 1.0
   szpoint(3) = 0.05*pi
   szpoint(4) = pi 
   
   
   redshift = 0.0d0
   call time_evolution(szpac,redshift,time)
   szpoint(5) = redshift 
   szpoint(6) = time
   szpoint(7) = 1.0
   print *, "time_evolution for redshift=",redshift,'is', time*szpac(25)


   open(101,file="fld")
   do I=1,100
     szpoint(2) = 0.5*I
     call fluid_variables(szpac,szpoint, npoint, fluid)
     write(101,*) szpoint(2),fluid
   enddo




end program szekeres_evolution

!------------------------------------

subroutine link_point(input_data,rho,tht,shr,wey,ric,arl,prp)
    implicit none


    double precision, dimension(0:37), intent(in) :: input_data
    
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer :: ngrid
    double precision :: dr
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point
    double precision, dimension (npoint) :: szpoint    
    double precision :: redshift,time
    double precision, dimension (100) :: szpac  
    double precision,  dimension (10) :: fluid
    double precision, intent(out) :: rho,tht,shr,wey,ric,arl,prp

    rho = 0.0d0
    tht = 0.0d0
    shr = 0.0d0
    wey = 0.0d0
    ric = 0.0d0
    ric = 0.0d0
    arl = 0.0d0
    prp = 0.0d0
    pypac(1:15)  =  input_data(0:14)
    pyszek(1:15) =  input_data(15:29)
    point(1:7)  =  input_data(30:36)
    ngrid = int(input_data(37))
    dr = point(2)/(1.0d0*ngrid)   
    redshift = point(5) 
    call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    call time_evolution(szpac,redshift,time)
    szpoint(1) = point(1)
    szpoint(2) = point(2)
    szpoint(3) = point(3)
    szpoint(4) = point(4)
    szpoint(5) = redshift
    szpoint(6) = time
    szpoint(7) = 1.0  
    szpac(100) = 5.0
    call fluid_variables(szpac,szpoint, npoint, fluid)
    rho = fluid(1)
    tht = fluid(2)
    shr = fluid(3)
    wey = fluid(4)
    ric = fluid(5)
    arl = fluid(6)
    prp = fluid(7)
 
    
    
end subroutine link_point


!------------------------------------
subroutine link_multi(input_data,rho,tht,shr,wey,ric,arl,prp)
    implicit none


    double precision, dimension(0:37), intent(in) :: input_data
    
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer :: ngrid, ig
    double precision :: dr


    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point

    double precision, dimension (npoint) :: szpoint    
    double precision :: redshift,time
    double precision, dimension (100) :: szpac  
    double precision,  dimension (10) :: fluid

    double precision,  dimension (100), intent(out) :: rho,tht,shr,wey,ric,arl,prp

    rho = 0.0d0
    tht = 0.0d0
    shr = 0.0d0
    wey = 0.0d0
    ric = 0.0d0
    ric = 0.0d0
    arl = 0.0d0
    prp = 0.0d0
    pypac(1:15)  =  input_data(0:14)
    pyszek(1:15) =  input_data(15:29)
    point(1:7)  =  input_data(30:36)
    ngrid = int(input_data(37))
    dr = point(2)/(1.0d0*ngrid)   
    redshift = point(5) 
    call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    call time_evolution(szpac,redshift,time)
    szpoint = point 
    szpoint(6) = time
    szpoint(7) = 1.0
    


!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(ig,szpac,szpoint,fluid) &
!$OMP SHARED(pypac,pyszek,point,ngrid,dr,redshift,time,input_data,rho,tht,shr,wey,ric,arl,prp)
    
    do ig=1,ngrid
        szpoint(1) = point(1)
        szpoint(2) = ig*dr
        szpoint(3) = point(3)
        szpoint(4) = point(4)
        szpoint(5) = redshift
        szpoint(6) = time
        szpoint(7) = 1.0  
        szpac(100) = 5.0
        call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
        call fluid_variables(szpac,szpoint, npoint, fluid)
        rho(ig) = fluid(1)
        tht(ig) = fluid(2)
        shr(ig) = fluid(3)
        wey(ig) = fluid(4)
        ric(ig) = fluid(5)
        arl(ig) = fluid(6)
        prp(ig) = fluid(7)
    enddo   
!$OMP END PARALLEL DO   
    
    
end subroutine link_multi

!------------------------------------
subroutine fluid_variables(szpac,point,npoint,fluid)
! calulates the fluid variables: 
! density (rho), expansion (tht), shear (shr), Weyl (wey), 3D Ricci (ric)
! based on eq. (3.2), (3.3), (3.4), (3.5), and (3.7) of arxiv:1704.02810
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: npoint
    double precision,  dimension (npoint) :: point
    double precision :: t,r,theta,phi
    double precision :: m,mi,mr,mrr,k,ki,kr,krr,q,qr,qrr,p,pr,prr,s,sr,srr,si,n,epr
    
    double precision :: redshift, zz2, zz3,zz4,ez,gzn,hzn
    double precision :: f12,f13,f16,f19    
    double precision :: rho1,rho2,tht1,tht2,shr1,shr2,wey1,wey2,ric1,ric2
    double precision :: density, expansion, shear, weyl, ricci,areal,proper
    double precision :: aRi,aRi2,aRi3,denumerator
    double precision :: aR, aRt, aRr, aRrr, aRtr, aRtrr, rpR  
    double precision,  dimension (10), intent(out) :: fluid
    
    t = point(1)
    r = point(2)
    theta = point(3)
    phi = point(4)
    
    call szekeres_specifics(szpac,r)
    
    m    = szpac(61) 
    mr   = szpac(62) 
    mrr  = szpac(63) 
    k    = szpac(64)
    kr   = szpac(65) 
    krr  = szpac(66) 
    q    = szpac(67) 
    qr   = szpac(68) 
    qrr  = szpac(69) 
    p    = szpac(70) 
    pr   = szpac(71) 
    prr  = szpac(72) 
    s    = szpac(73) 
    sr   = szpac(74) 
    srr  = szpac(75) 

    redshift = point(5)
    zz2 = (1.0 + redshift)**2
    zz3 = (1.0 + redshift)**3
    zz4 = (1.0 + redshift)**4
    ez = dsqrt(szpac(7)*zz4 + szpac(3)*zz3 + szpac(9)*zz2 + szpac(2))
    gzn = 1.0d0/zz3
    hzn = 1.0d0/ez

    call areal_radius(szpac,point,npoint,aR,aRt, aRr, aRrr, aRtr, aRtrr)
    call radial_proper_distance(szpac,point,npoint,aR,aRr,rpR)
   
    areal = aR
    proper = rpR

    aRi = 1.0d0/aR
    aRi2 = aRi*aRi
    aRi3 = aRi*aRi*aRi
    f12 = 0.5d0
    f13 = 1.0d0/3.0d0
    f16 = f12*f13
    f19 = f13*f13
    
    
    ki = 1.0d0/k
    mi = 1.0d0/m
    si = -1.0d0/s
    
    n = pr*cos(phi) + qr*sin(phi)
    epr = si*(sr*cos(theta) + n*sin(theta))
    denumerator = 1.0d0/(aRr - aR*epr)

    rho1 = 2.0d0*(mr-3.0d0*m*epr) 
    rho2 = aRi2*denumerator
    density = rho1*rho2*szpac(23)*gzn 
    
    tht1 = aRtr + 2*aRt*aRr*aRi - 3.0d0*aRt*epr
    tht2 = denumerator*f13
    expansion = tht1*tht2*szpac(21)*hzn

    shr1 = -(aRtr - aRt*aRr*aRi)
    shr2 = f13*denumerator
    shear = shr1*shr2*szpac(21)*hzn
    
    wey1 = m*(3*aRr - aR*mr*mi) 
    wey2 = f13*aRi3*denumerator
    weyl = wey1*wey2*szpac(23)*gzn 

    ric1 = 2.0*k*aRi2
    ric2 = 1.0d0 + (aR*kr*ki - 2.0*aR*epr)*denumerator
    ricci = ric1*ric2


    fluid(1) = density
    fluid(2) = expansion
    fluid(3) = shear
    fluid(4) = weyl
    fluid(6) = ricci
    fluid(7) = areal  * 1d-3   ! changing units from Kpc to Mpc
    fluid(8) = proper * 1d-3   ! changing units from Kpc to Mpc

	
end subroutine fluid_variables
!--------------------------------------------------

subroutine areal_radius(szpac,point,npoint,aR,aRt, aRr, aRrr, aRtr, aRtrr)
! calulates areal distace R and its derivatives  
! areal distance at time t: aR
! time derivatice areal distance at time t: aRt
! first radial derivatice areal distance at time t: aRr
! time derivative of the first radial derivatice areal distance at time t: aRtr
! second radial derivatice areal distance at time t: aRrr
! time derivative of the second radial derivatice areal distance at time t: aRtrr
! based on  the unpublished notes: "sphesze04.pdf"
! start with initial conditions (40) via evolution of (37)-(39)

    implicit none
        
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: npoint
    double precision,  dimension (npoint) :: point
    integer :: I,Ni
    double precision :: dti,dt,time
    double precision :: aR,aRt, aRr, aRrr, aRtr, aRtrr
    double precision :: q,qr,qrr,p,pr,prr,s,sr,srr
    double precision :: m,mr,mrr,k,kr,krr,l
    double precision :: kr1,kr2,kr3,kr4
    double precision :: krr1,krr2,krr3,krr4
    double precision :: krrr1,krrr2,krrr3,krrr4
    double precision :: r,rp,rt
    double precision :: rr,rrp,rtr
    double precision :: rrr,rrrp,rtrr


    if(point(7)<0.5) then   
       time = point(1)
    else
       time = point(6)
    endif
    
    dti = 105.0
    Ni = int( abs(time/dti) )
    if(Ni.le.3) Ni = 3
    l = szpac(26)
    dt = time/(Ni*1.0d0)
	
	
    r = point(2)
    rr = 1.0d0
    rrr = 0.0d0

    m    = szpac(61) 
    mr   = szpac(62) 
    mrr  = szpac(63) 
    k    = szpac(64)
    kr   = szpac(65) 
    krr  = szpac(66) 
    q    = szpac(67) 
    qr   = szpac(68) 
    qrr  = szpac(69) 
    p    = szpac(70) 
    pr   = szpac(71) 
    prr  = szpac(72) 
    s    = szpac(73) 
    sr   = szpac(74) 
    srr  = szpac(75) 
    

    do I=1,Ni
      rp = r
      rrp = rr
      rrrp = rrr
      rt = dsqrt(2.0d0*(m/rp) - k + (l/3.0d0)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l/3d0)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l/3.0d0)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr1   = dt*rt
      krr1  = dt*rtr
      krrr1 = dt*rtrr
      rp   = r   + kr1*5d-1
      rrp  = rr  + krr1*5d-1
      rrrp = rrr + krrr1*5d-1
      rt = dsqrt(2.0d0*(m/rp) - k + (l/3.0d0)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l/3d0)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l/3.0d0)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr2   = dt*rt
      krr2  = dt*rtr
      krrr2 = dt*rtrr
      rp   = r   + kr2*5d-1
      rrp  = rr  + krr2*5d-1
      rrrp = rrr + krrr2*5d-1
      rt = dsqrt(2.0d0*(m/rp) - k + (l/3.0d0)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l/3d0)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l/3.0d0)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr3   = dt*rt
      krr3  = dt*rtr
      krrr3 = dt*rtrr
      rp   = r   + kr3
      rrp  = rr  + krr3
      rrrp = rrr + krrr3
      rt = dsqrt(2.0d0*(m/rp) - k + (l/3.0d0)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l/3d0)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l/3.0d0)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr4   = dt*rt
      krr4  = dt*rtr
      krrr4 = dt*rtrr
      rp   = r   + (kr1  +2.0d0*kr2  +2.0d0*kr3  +kr4)/6.0d0
      rrp  = rr  + (krr1 +2.0d0*krr2 +2.0d0*krr3 +krr4)/6.0d0
      rrrp = rrr + (krrr1+2.0d0*krrr2+2.0d0*krrr3+krrr4)/6.0d0
      r   = rp
      rr  = rrp
      rrr = rrrp
    enddo

    rt = dsqrt(2.0d0*(m/rp) - k + (l/3.0d0)*rp*rp)
    rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l/3d0)*rp*rrp)/rt
    rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
    rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
    rtrr = rtrr + (l/3.0d0)*(rrp**2 + rp*rrrp) - rtr*rtr
    rtrr = rtrr/rt
    aR    = r
    aRt   = rt
    aRr   = rr
    aRrr  = rrr
    aRtr  = rtr
    aRtrr = rtrr

end subroutine areal_radius
!--------------------------------------------------

subroutine radial_proper_distance(szpac,point,npoint,aR,aRr,rpR)
    implicit none
    double precision,  intent(inout), dimension (100) :: szpac 
    integer, intent(in) :: npoint
    double precision,  dimension (npoint) :: point
    double precision, intent(in) :: aR,aRr
    double precision, intent(out)  :: rpR
    integer :: I,Ni
    double precision :: r, theta, phi, n
    double precision :: m,mr,mrr,k,ki1,kr,krr,q,qr,qrr,p,pr,prr,s,sr,srr,si
    double precision :: rout,rint,ri,dr,kr1,kr2,kr3,kr4,f16,f

    r = point(2)
    theta = point(3)
    phi = point(4)

    f16 = 1.0d0/6.0d0
    Ni = 124
    dr = r/(1.0d0*Ni)
    ri = dr
    rout = 0.0d0

    do I =1,Ni
      rint = ri 
      call szekeres_specifics(szpac,rint)
      m    = szpac(61) 
      mr   = szpac(62) 
      mrr  = szpac(63) 
      k    = szpac(64)
      kr   = szpac(65) 
      krr  = szpac(66) 
      q    = szpac(67) 
      qr   = szpac(68) 
      qrr  = szpac(69) 
      p    = szpac(70) 
      pr   = szpac(71) 
      prr  = szpac(72) 
      s    = szpac(73) 
      sr   = szpac(74) 
      srr  = szpac(75) 
      n = pr*cos(phi) + qr*sin(phi)
      si = 1.0d0/s
      ki1 = 1.0d0/(1.0d0-k)
      f = (abs((sr*cos(theta) + n*sin(theta))*aR*si + aRr))*sqrt(ki1)
      kr1 = dr*f
      rint = ri + 0.5d0*dr
      call szekeres_specifics(szpac,rint)
      m    = szpac(61) 
      mr   = szpac(62) 
      mrr  = szpac(63) 
      k    = szpac(64)
      kr   = szpac(65) 
      krr  = szpac(66) 
      q    = szpac(67) 
      qr   = szpac(68) 
      qrr  = szpac(69) 
      p    = szpac(70) 
      pr   = szpac(71) 
      prr  = szpac(72) 
      s    = szpac(73) 
      sr   = szpac(74) 
      srr  = szpac(75) 
      n = pr*cos(phi) + qr*sin(phi)
      si = 1.0d0/s
      ki1 = 1.0d0/(1.0d0-k)
      f = (abs((sr*cos(theta) + n*sin(theta))*aR*si + aRr))*sqrt(ki1)
      kr2 = dr*f
      kr3 = dr*f
      rint = ri + dr
      call szekeres_specifics(szpac,rint)
      m    = szpac(61) 
      mr   = szpac(62) 
      mrr  = szpac(63) 
      k    = szpac(64)
      kr   = szpac(65) 
      krr  = szpac(66) 
      q    = szpac(67) 
      qr   = szpac(68) 
      qrr  = szpac(69) 
      p    = szpac(70) 
      pr   = szpac(71) 
      prr  = szpac(72) 
      s    = szpac(73) 
      sr   = szpac(74) 
      srr  = szpac(75) 
      n = pr*cos(phi) + qr*sin(phi)
      si = 1.0d0/s
      ki1 = 1.0d0/(1.0d0-k)
      f = (abs((sr*cos(theta) + n*sin(theta))*aR*si + aRr))*sqrt(ki1)
      kr4 = dr*f
      rout = rout + f16*(kr1+2.0d0*(kr2+kr3)+kr4)
      ri = rint
    enddo

rpR = rout

end subroutine radial_proper_distance
!--------------------------------------------------
subroutine look_back_time(szpac,redshift,time)

! calculates the lookback time, BUT....
! the convention here is that its is negative
! output is in Myr

    implicit none
    integer :: I,Ni
    double precision,  dimension (100), intent(in) :: szpac  
    double precision, intent(in) :: redshift
    double precision, intent(out) :: time
    double precision :: omega_matter,omega_lambda,omega_curvature,omega_radiation
    double precision :: Ho, Hzz,dz,zz,zi,z
    double precision :: t_rk1,t_rk2,t_rk3,t_rk4,ti,t,f16
    Ho = szpac(11)
    omega_matter =  szpac(3)
    omega_lambda =  szpac(2)
    omega_radiation = szpac(7)
    omega_curvature = szpac(9)
    if(redshift > 10.0) then
       print *, "--------------------- warning ---------------------------------"
       print *, "! one should not use the function lookback for high redshifts !"
       print *, "! one should use a different function, or accept the errors   !"
       print *, "--------------------- warning ---------------------------------"
    endif
    Ni = 1024
    dz = redshift/(Ni*1d0)
    zi = 0d0
    f16 = 1.0d0/6.0d0
    ti = 0d0
    t = ti
    do I=1,Ni
      z = zi 
      zz = 1.0d0 + z
      Hzz=zz*Ho*dsqrt(omega_radiation*(zz**4)+omega_matter*(zz**3)+omega_curvature*(zz**2)+omega_lambda)
      t_rk1 = dz/Hzz
      z = zi + 5.0d-1*dz
      zz = 1.0d0 + z
      Hzz=zz*Ho*dsqrt(omega_radiation*(zz**4)+omega_matter*(zz**3)+omega_curvature*(zz**2)+omega_lambda)
      t_rk2 = dz/Hzz
      t_rk3 = dz/Hzz
      z = zi + dz
      zz = 1.0d0 + z
      Hzz=zz*Ho*dsqrt(omega_radiation*(zz**4)+omega_matter*(zz**3)+omega_curvature*(zz**2)+omega_lambda)
      t_rk4 = dz/Hzz
      t = ti + f16*(t_rk1+2.0d0*(t_rk2+t_rk3)+t_rk4)
      ti = t
      zi = z
    enddo
    time = -t*szpac(25)
end subroutine look_back_time
!--------------------------------------------------
subroutine time_evolution(szpac,redshift,time)
! time_evolution time is evaluated from the initial redshift
! it provide the amount of time needed to evolve the strcuture 
! from the initial instant (eg. CMB) to the a given redshift 
! the output is in Kpc, i.e. length unit: c * time
    implicit none
    double precision,  dimension (100), intent(in) :: szpac  
    double precision, intent(in) :: redshift
    double precision, intent(out) :: time
    double precision :: lookback
    call look_back_time(szpac,redshift,lookback)
    time = szpac(10) + lookback*szpac(24)
end subroutine time_evolution
!--------------------------------------------------
subroutine age_from_initial(szpac,z_initial,time)
    implicit none
    integer :: I,Ni
    double precision,  dimension (100), intent(in) :: szpac  
    double precision, intent(in) :: z_initial
    double precision, intent(out) :: time
    double precision :: omega_matter,omega_lambda,omega_curvature,omega_radiation
    double precision :: Ho, Haa,da,ai,a,ia1,ia2,ia3,ia4
    double precision :: t_rk1,t_rk2,t_rk3,t_rk4,ti,t,f16 
    Ho = szpac(11)
    omega_matter =  szpac(3)
    omega_lambda =  szpac(2)
    omega_radiation = szpac(7)
    omega_curvature = szpac(9)
    Ni = 4096
    f16 = 1.0d0/6.0d0
    ti = 0d0
    t = ti
    ai = 1.0d0/(z_initial + 1.0d0)
    da = (1.0d0-ai) / (Ni*1.0d0)
    do I=1,Ni
        a = ai 
        ia1 = 1.0d0/a
        ia2 = ia1*ia1
        ia3 = ia2*ia1
        ia4 = ia2*ia2
        Haa=a*Ho*dsqrt(omega_radiation*ia4+omega_matter*ia3+omega_curvature*ia2+omega_lambda)
        t_rk1 = da/Haa
        a = ai + 5.0d-1*da
        ia1 = 1.0d0/a
        ia2 = ia1*ia1
        ia3 = ia2*ia1
        ia4 = ia2*ia2
        Haa=a*Ho*dsqrt(omega_radiation*ia4+omega_matter*ia3+omega_curvature*ia2+omega_lambda)
        t_rk2 = da/Haa
        t_rk3 = da/Haa
        a = ai + da
        ia1 = 1.0d0/a
        ia2 = ia1*ia1
        ia3 = ia2*ia1
        ia4 = ia2*ia2
        Haa=a*Ho*dsqrt(omega_radiation*ia4+omega_matter*ia3+omega_curvature*ia2+omega_lambda)
        t_rk4 = da/Haa
        t = ti + f16*(t_rk1+2.0d0*(t_rk2+t_rk3)+t_rk4)
        ti = t
        ai = a
    enddo
    time = t 
end subroutine age_from_initial
!--------------------------------------------------
subroutine redlcdm(szpac,time,redshift)
    implicit none
    double precision,  dimension (100), intent(in) :: szpac  
    double precision, intent(in) :: time
    double precision, intent(out) :: redshift
    double precision :: ctt,lambda,rhoz,rhoi,zo
    ctt = time*szpac(24)
    lambda =  szpac(26)
    rhoz = lambda/((sinh(ctt*(dsqrt(75d-2*lambda))))**2)
    rhoi = szpac(23)
    zo = ((rhoz*rhoi)**(1.0d0/3.0d0)) - 1.0d0
    redshift = zo
end subroutine redlcdm
!--------------------------------------------------
subroutine time_lcdm(szpac,redshift,time)
    implicit none
    double precision,  dimension (100), intent(in) :: szpac  
    double precision, intent(in) :: redshift
    double precision, intent(out) :: time
    double precision :: zi,lambda,rhoz,x,arsh,tz 
    zi = redshift
    lambda =  szpac(26)
    rhoz = szpac(22)*( (1.0d0 + zi )**3 )
    x = dsqrt(lambda/rhoz)
    arsh = dlog(x + dsqrt(x*x + 1.0d0))
    tz = (dsqrt((4.0d0)/(3.0d0*lambda)))*arsh 
    time = tz*szpac(25)      
end subroutine time_lcdm
!------------------------------------

subroutine szekeres_specifics(szpac,r)
    implicit none
    double precision, intent(inout), dimension (100) :: szpac    
    double precision, intent(in) :: r
    double precision :: m,mr,mrr  
    double precision :: k,kr,krr
    double precision :: s,sr,srr    
    double precision :: p,pr,prr    
    double precision :: q,qr,qrr        
    double precision :: dta0,dta1,dta2,d0el,d1el,d2el
    double precision :: r0,dlr,Om,Ak,Am,amp,ho,alpha

    Om = szpac(3)
    amp = szpac(51)
    ho = szpac(11)
    r0 = szpac(52)
    dlr = szpac(53)*r0
    
    
    Am = (1.0d0/6.0d0)*szpac(31)
    Ak = (14.0d0/3.0d0)*Am
 

    dta0 = dtanh((r-r0)/(2.0d0*dlr))
    dta1 = (1.0d0 - dta0**2)*(1.0d0/(2.0d0*dlr))
    dta2 = -dta0*dta1/dlr
    d0el = amp*0.5d0*(1.0d0 - dta0)
    d1el = -amp*0.5d0*dta1
    d2el = -amp*0.5d0*dta2

    m = Am*(1 + d0el)*r*r*r
    mr = Am*d1el*r*r*r + 3.0d0*(m/r)
    mrr = Am*d2el*(r**3)+3d0*Am*d1el*r*r+3d0*(mr/r)-3.0d0*(m/(r*r))

    k = Ak*d0el*r*r
    kr = Ak*d1el*r*r + 2.0*(k/r)
    krr = Ak*d2el*r*r + 2.0*Ak*d1el*r + 2.0*kr/r - 2.0*(k/(r*r))
    
    alpha = szpac(54) 

    q = 1.0d0
    qr = 0.0d0
    qrr = 0.0d0

    p = 1.0d0
    pr = 0.0d0
    prr = 0.0d0

    if (alpha < 0.01 .or. alpha > 0.99) then
       s = 1.0d0
       sr = 0.0d0
       srr = 0.0d0
    else
       s = r**alpha
       sr = alpha*(r**(alpha-1.0))
       srr = alpha*(alpha-1.0)*(r**(alpha-2.0))
    endif


    szpac(60) = r
    szpac(61) = m
    szpac(62) = mr
    szpac(63) = mrr
    szpac(64) = k
    szpac(65) = kr
    szpac(66) = krr
    szpac(67) = q
    szpac(68) = qr
    szpac(69) = qrr
    szpac(70) = p
    szpac(71) = pr
    szpac(72) = prr
    szpac(73) = s
    szpac(74) = sr
    szpac(75) = srr




end subroutine szekeres_specifics
!--------------------------------------------------

subroutine parameter_names(szpan)
    implicit none
    character(len=8), dimension(100) :: szpan

    szpan = ''
    szpan(1) = "H0"
    szpan(2) = "Ol0"
    szpan(3) = "Om0"

end subroutine parameter_names

!------------------------------------
subroutine parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    implicit none
    integer :: npypac, npyszek
    double precision,  dimension (npypac) :: pypac
    double precision,  dimension (npyszek) :: pyszek        
    double precision,  dimension (100) :: szpac    


    double precision :: H0,Ho,little_h,age,z_initial
    double precision :: contrast, radius, slope, dipole    
    double precision :: omega_matter,omega_baryon,omega_cold
    double precision :: omega_lambda,omega_photon,omega_radiation,omega_curvature
    double precision :: w_matter,w_lambda,w_baryon,w_cold,w_photon,w_radiation
    double precision :: mass_unit,length_unit,time_unit,pi,G0,c0,T0,Neff,sbc,gcons,light_speed,kap,kapc2,gkr,gcr,lambda
    

! szpac:: contains parameters by convention:
! szpac(1-50)   -> parameters relate to homogeneous background
! szpac(51-90)  -> parameters relate to inhomogeneous setup 
! szpac(91-100) -> reserve for special flags 
    
  

! from python:    
    if(szpac(100) > 1.0) then
       H0 = pypac(1)
       omega_matter = pypac(2)
       omega_lambda = pypac(3)    
       
       contrast = pyszek(1)
       radius = pyszek(2)
       slope = pyszek(3)
       dipole = pyszek(4)
       
    else
       H0 = 70.0d0
       omega_matter = 0.3
       omega_lambda = 1.0d0-omega_matter
       
       contrast = -0.0025
       radius = 10.0
       slope  = 0.4
       dipole = 0.4
       
    endif



    
! units: time in 10^6 years, length in Kpc, mass in 10^15 M_{\odot}
    mass_unit=1.989d45
    length_unit=3.085678d19
    time_unit=31557600*1d6
! and other constants
    little_h = H0*0.01d0
    pi = 4.0d0*datan(1.0d0)
    G0 = 6.67430d-11
    c0 = 299792458.0d0
    T0 = 2.728d0                 *    0.0
    Neff = 3.046d0
    sbc = 5.670374419d-8
    gcons= G0*((mass_unit*(time_unit**2))/(length_unit**3))
    light_speed=c0*(time_unit/length_unit)
    kap=8.0d0*pi*gcons*(1d0/(light_speed*light_speed*light_speed*light_speed))
    kapc2=8d0*pi*gcons*(1d0/(light_speed*light_speed))
    Ho=(time_unit/(length_unit))*H0
    gkr=3.0d0*(((Ho*Ho))/(8.0d0*pi*gcons))
    lambda=3.0d0*omega_lambda*(((Ho*Ho))/(light_speed*light_speed))
    gkr=kapc2*gkr*omega_matter
    gcr = gkr/kapc2    
    w_photon = 8.0d0*pi*G0*4*sbc*(T0**4)/(c0*c0*c0*3*(100.0/length_unit)**2)
    w_matter = omega_matter/(little_h*little_h)
    w_lambda = omega_lambda/(little_h*little_h)    
    omega_photon = w_photon/(little_h*little_h)
    omega_radiation = omega_photon* (1.0d0 + (7.0d0/8.0d0)*Neff*((4.0d0/11.0d0)**(4.0d0/3.0d0)))  
    w_radiation = w_photon* (1.0d0 + (7.0d0/8.0d0)*Neff*((4.0d0/11.0d0)**(4.0d0/3.0d0)))  
    omega_curvature = 1.0d0 - omega_matter - omega_lambda - omega_radiation
    z_initial = 1089.80  
    
    omega_baryon = 0.0
    omega_cold = 0.0 
    w_baryon = omega_baryon/(little_h*little_h)
    w_cold = omega_cold/(little_h*little_h)  

    szpac(1) = H0
    szpac(2) = omega_lambda
    szpac(3) = omega_matter
    szpac(4) = omega_baryon 
    szpac(5) = omega_cold   
    szpac(7) = omega_radiation
    szpac(9) = omega_curvature
    
    szpac(11)  = Ho/light_speed
    szpac(12)  = w_lambda
    szpac(13)  = w_matter
    szpac(14)  = w_baryon    
    szpac(15)  = w_cold
    szpac(16)  = w_photon 
    szpac(17)  = w_radiation
    
    
    szpac(21)  = light_speed/Ho 
    szpac(22)  = gkr
    szpac(23)  = 1.0d0/gkr    
    szpac(24)  = light_speed
    szpac(25)  = 1.0d0/light_speed
    szpac(26)  = lambda

    szpac(30) = z_initial
    szpac(31) = gkr*( (1.0d0 + z_initial)**3)

    szpac(51) = contrast
    szpac(52) = radius    
    szpac(53) = slope
    szpac(54) = dipole

    call age_from_initial(szpac,z_initial, age)
    szpac(10) = age
    
end subroutine parameter_values



