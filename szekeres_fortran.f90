!########################################################################################################
! SzekeresPy ver. 0.18 - Python package for cosmological calculations using the Szekeres Cosmological Model
! 
! File: szekeres_fortran.f90
! 
! Author: Krzysztof Bolejko
! 
! Intended use: research and education
! 
! Licence to use: BSD-2-Clause ("FreeBSD License")
! 
! Copyright (c) 2024 Krzysztof Bolejko
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 
! THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!########################################################################################################


program szekeres_fortran
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




end program szekeres_fortran

!------------------------------------

subroutine link_point(input_data,rho,tht,shr,wey,ric,arl,prp)
    implicit none


    double precision, dimension(0:37), intent(in) :: input_data
    
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer, parameter :: ngrid = 100
    integer :: ngrid000
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point
    double precision, dimension (npoint) :: szpoint    
    double precision :: redshift,time
    double precision, dimension (100) :: szpac  
    double precision,  dimension (10) :: fluid
    double precision, intent(out) :: rho,tht,shr,wey,ric,arl,prp

    szpac(100) = 2


    rho = 0.0d0
    tht = 0.0d0
    shr = 0.0d0
    wey = 0.0d0
    ric = 0.0d0
    ric = 0.0d0
    arl = 0.0d0
    prp = 0.0d0

    ngrid000 = int(input_data(0))
    pypac(1:15)  =  input_data(1:15)
    pyszek(1:15) =  input_data(16:30)
    point(1:7)  =  input_data(31:37)
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

    integer, parameter :: ngrid = 100

    integer :: ngrid000, ig
    double precision :: dr


    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point

    double precision, dimension (npoint) :: szpoint    
    double precision :: redshift,time
    double precision, dimension (100) :: szpac  
    double precision,  dimension (10) :: fluid

    double precision, dimension (0:(ngrid-1)), intent(out) :: rho,tht,shr,wey,ric,arl,prp


    szpac(100) = 2

    rho = 0.0d0
    tht = 0.0d0
    shr = 0.0d0
    wey = 0.0d0
    ric = 0.0d0
    ric = 0.0d0
    arl = 0.0d0
    prp = 0.0d0

    ngrid000 = int(input_data(0))
    pypac(1:15)  =  input_data(1:15)
    pyszek(1:15) =  input_data(16:30)
    point(1:7)  =  input_data(31:37)


    dr = point(2)/(1.0d0*ngrid)   
    redshift = point(5) 
    call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    call time_evolution(szpac,redshift,time)
    szpoint = point 
    szpoint(6) = time
    szpoint(7) = 1.0
    


!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(ig,szpac,szpoint,fluid) &
!$OMP SHARED(pypac,pyszek,point,dr,redshift,time,input_data,rho,tht,shr,wey,ric,arl,prp)
    
    do ig=0,ngrid-1
        szpoint(1) = point(1)
        szpoint(2) = (ig+1)*dr
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
subroutine link_null(input_data, temporal, radial, thetal, phial,redshiftal)
    implicit none
    double precision, dimension(0:44), intent(in) :: input_data
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer, parameter :: ngrid = 100
    integer :: ngrid000, ig
    double precision :: redshift,ngz
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point,direction
    double precision, dimension (100) :: szpac  
    double precision, dimension (0:(ngrid-1)), intent(out) :: temporal, radial, thetal, phial,redshiftal

    szpac(100) = 2
    temporal = 0.0d0
    radial = 0.0d0
    thetal = 0.0d0
    phial = 0.0d0
   
    ngrid000 = int(input_data(0))

    pypac(1:15)  =  input_data(1:15)
    pyszek(1:15) =  input_data(16:30)
    point(1:7)  =  input_data(31:37)
    direction(1:7)  =  input_data(38:44)
    redshift = point(5) 
    ngz = redshift/(ngrid*1.0d0 - 1.0d0)
    do ig = 0,ngrid-1
        redshiftal(ig) = ig*ngz
    enddo
    
    call parameter_values(npypac,pypac,npyszek,pyszek,szpac)
    call light_propagation(szpac,ngrid,npoint,point,direction,temporal,radial,thetal,phial,redshiftal)
    
end subroutine link_null

!------------------------------------
subroutine light_propagation(szpac,ngrid,npoint,point,direction,temporal,radial,thetal,phial,redshiftal)
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: ngrid, npoint
    double precision, dimension (npoint) :: point, direction    
    double precision, dimension (0:(ngrid-1)), intent(in) :: redshiftal
    double precision, dimension (0:(ngrid-1)), intent(out) :: temporal, radial, thetal, phial
    double precision, dimension(4) :: PV,PVi,PVii,NV,NVi,NVii,AV
    double precision, dimension(7,4) :: PRK,NRK
    integer :: Ui,I
    double precision :: ds,dss


    temporal = 0.0
    radial = 0.0 
    thetal = 0.0
    phial = 0.0
      

    Ui = 2500000
	dss = 0.0001d0
	ds = dss


    call initial_conditions(szpac,npoint,point,direction,PV,NV,AV)

    PVi=PV
	NVi=NV
	PVii=PV
	NVii=NV
	PRK = 0.0d0
	NRK = 0.0d0





end subroutine light_propagation
!------------------------------------
subroutine initial_conditions(szpac,npoint,point,direction,PV,NV,AV)
    implicit none
    double precision,  intent(inout), dimension (100) :: szpac 
    integer, intent(in) :: npoint    
    double precision, dimension (npoint) :: point, direction
    double precision, intent(out), dimension(4) :: PV,NV,AV
    double precision, dimension (npoint) :: DRV 
    double precision, dimension(4,4)   :: GIJ
    double precision, dimension(4,4,4) :: GAM
    double precision :: r

    PV(1) = point(1)
    PV(2) = point(2)
    PV(3) = point(3)
    PV(4) = point(4)



	call christoffels(szpac,PV,DRV,GIJ,GAM)
	call null_initial(szpac,npoint,point,direction,DRV,GIJ,NV)

    




end subroutine initial_conditions
!------------------------------------
subroutine null_initial(szpac,npoint,point,direction,DRV,GIJ,NV)
    implicit none
    double precision,  intent(inout), dimension (100) :: szpac 
    integer, intent(in) :: npoint    
    double precision, dimension (npoint) :: point, direction 
    double precision, dimension(npoint), intent(in)   :: DRV
    double precision, dimension(4,4), intent(in)   :: GIJ
    double precision, dimension(4), intent(out) :: NV
    integer :: I,J
    double precision :: k,kr,qr,pr,s,sr,si,s3,s4,c3,c4,s3i,ri,mkis
    double precision :: RA,DEC,th,ph,wt,bN,bN4,W,Wi

    RA = direction(1)
    DEC = direction(2)
    th = point(3)
    ph = point(4)
    k    = szpac(64)
    kr   = szpac(65) 
    qr   = szpac(68) 
    pr   = szpac(71) 
    s    = szpac(73) 
    sr   = szpac(74) 
    si = 1.0/s
    s3= sin(th)
    c3= cos(th)
    s4= sin(ph)
    c4= cos(ph)
    s3i = 1.0/s3
    ri = 1.0d0/DRV(1)
    mkis = 1.0/sqrt(1.0d0 - k)
    bN = pr*c4+qr*s4  ! big N
    bN4 = -pr*s4 + qr*c4 ! derivative of big N with respect to phi
    W = mkis*(DRV(2)+ DRV(1)*si*(sr*c3+bN*s3) )  !big delta
    Wi = 1.0d0/W

    NV(1) = 1.0d0
    NV(2) = dcos(RA)*dcos(DEC)*Wi
    NV(3) = dcos(RA)*dcos(DEC)*(sr*s3+bN*(1d0-c3)*Wi*si) - ri*dsin(DEC)
    NV(4) = -NV(1)*bN4*(1d0-c3)+ri*s3i*dsin(RA)*dcos(DEC) 

	NV = -1.0*NV
	wt = 0d0
	do J=2,4
		do I=2,4
		wt = wt + (GIJ(I,J)*NV(I)*NV(J))
		enddo
	enddo
	NV = NV/dsqrt(dabs(wt))
	NV(1) = -1.0d0


end subroutine null_initial
    






!------------------------------------
subroutine christoffels(szpac,position_vector,derivatives,metric,christoffel)
    implicit none
    integer, parameter :: npoint = 7
    double precision,  intent(inout), dimension (100) :: szpac 
    double precision, dimension(4), intent(in)        :: position_vector
    double precision, dimension(npoint)                  :: point = 0.0d0  
    double precision, dimension(npoint), intent(out)   :: derivatives 
    double precision, dimension(4,4), intent(out)   :: metric 
    double precision, dimension(4,4,4), intent(out) :: christoffel 
    double precision :: t,r,th,ph,s3,c3,s4,c4,mc3,mk
    double precision :: PSr,QSr,SSr,PSrr,QSrr,SSrr
    double precision :: p2,ppt,ppr,pprt,pprr
    double precision :: fp,f1,fp1,f3,fr3,f4,fp4,f5
    double precision :: grr,gr3,gr4,rr3,rp3,frrr
    double precision :: aR,aRt, aRr, aRrr, aRtr, aRtrr
    double precision :: ro0,rr,rt,rtr,rrr
    double precision :: m,mr,mrr,k,kr,krr,q,qr,qrr,p,pr,prr,s,sr,srr

    t  = position_vector(1)
    r  = position_vector(2)
    th = position_vector(3)
    ph = position_vector(4)
    point(1:4) = position_vector(1:4)
    call areal_radius(szpac,point,npoint,aR,aRt, aRr, aRrr, aRtr, aRtrr)
    ro0 = aR
    rr  = aRr
    rt  = aRt
    rtr = aRtr
    rrr = aRrr
    derivatives(1) = aR
    derivatives(2) = aRr
    derivatives(3) = aRrr
    derivatives(4) = aRt
    derivatives(5) = aRtr
    derivatives(6) = aRtrr
    derivatives(7) = r

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
    s3= sin(th)
    c3= cos(th)
    s4= sin(ph)
    c4= cos(ph)
    mc3= 1.0d0-c3
    mk= 1.0d0-k
    kr= kr
    PSr= pr/s
    QSr= qr/s
    SSr= sr/s
    PSrr= prr/s
    QSrr= qrr/s
    SSrr= srr/s
    p2= ro0*ro0
    ppt= rt/ro0
    ppr= rr/ro0
    pprt= rtr/ro0
    pprr= rrr/ro0
    f1= PSr*c4+QSr*s4
    fp1= PSrr*c4+QSrr*s4
    f3= QSr*c4-PSr*s4
    fr3= QSrr*c4-PSrr*s4-f3*SSr
    f4= f1*s3+SSr*c3
    fp4= fp1*s3+SSrr*c3
    f5= SSr**2*s3**2+2.0d0*SSr*f1*s3*mc3+mc3**2*(PSr**2+QSr**2)
    fp= ppr+f4
    grr= f5+fp**2/mk
    gr3= f1*mc3+SSr*s3
    gr4= mc3*s3*f3
    rr3= (fp1-f1*SSr)*mc3 + (SSrr-SSr**2)*s3
    rp3= f1*c3-SSr*s3
    frrr= 5d-1*kr/mk 
    frrr=frrr+(pprr-ppr**2+fp4-SSr*f4-gr3*rp3+mc3*f3**2-f5*mk)/fp
    metric=0.0d0
    christoffel=0.0d0
    metric(1,1)= 1.0d0
    metric(2,2)= -p2*grr
    metric(2,3)=  p2*gr3
    metric(3,2)= metric(2,3)
    metric(2,4)= -p2*gr4
    metric(4,2)= metric(2,4)
    metric(3,3)= -p2
    metric(4,4)= metric(3,3)*s3**2
    christoffel(1,2,2)= p2*(ppt*f5+(pprt+ppt*f4)*fp/mk)
    christoffel(1,2,3)= -ppt*metric(2,3)
    christoffel(1,3,2)= -ppt*metric(2,3)
    christoffel(1,3,3)= -ppt*metric(3,3)
    christoffel(1,2,4)= -ppt*metric(2,4)
    christoffel(1,4,2)= -ppt*metric(2,4)
    christoffel(1,4,4)= -ppt*metric(4,4)
    christoffel(2,1,2)= (pprt+ppt*f4)/fp
    christoffel(2,2,1)= (pprt+ppt*f4)/fp
    christoffel(2,2,2)=  ppr+ frrr
    christoffel(2,2,3)= (gr3*mk+rp3)/fp
    christoffel(2,3,2)= (gr3*mk+rp3)/fp
    christoffel(2,2,4)= (1.0d0-mc3*mk)*s3*f3/fp
    christoffel(2,4,2)= (1.0d0-mc3*mk)*s3*f3/fp
    christoffel(2,3,3)= -mk/fp
    christoffel(2,4,4)= -mk*s3**2/fp
    christoffel(3,1,2)= (pprt-ppr*ppt)*gr3/fp
    christoffel(3,2,1)= (pprt-ppr*ppt)*gr3/fp
    christoffel(3,1,3)= ppt
    christoffel(3,3,1)= ppt
    christoffel(3,2,2)= gr3*frrr -fp*(gr3+rp3/mk) -rr3 -mc3*s3*f3**2
    christoffel(3,2,3)= gr3*christoffel(2,2,3) +ppr
    christoffel(3,3,2)= gr3*christoffel(2,2,3) +ppr
    christoffel(3,2,4)= gr3*christoffel(2,2,4) -s3**2*f3
    christoffel(3,4,2)= gr3*christoffel(2,2,4) -s3**2*f3
    christoffel(3,3,3)= -mk*gr3/fp
    christoffel(3,4,4)= gr3*christoffel(2,4,4) -s3*c3
    christoffel(4,1,2)= f3*mc3/s3*(ppt-christoffel(2,1,2))
    christoffel(4,2,1)= f3*mc3/s3*(ppt-christoffel(2,1,2))
    christoffel(4,1,4)= ppt
    christoffel(4,4,1)= ppt
    christoffel(4,2,2)= (mc3*(f3*(ppr-frrr-SSr) +fr3)-fp/mk*f3)/s3
    christoffel(4,2,3)= f3-mc3*f3/s3*christoffel(2,2,3)
    christoffel(4,3,2)= f3-mc3*f3/s3*christoffel(2,2,3)
    christoffel(4,2,4)= ppr-mc3*f3/s3*christoffel(2,2,4)
    christoffel(4,4,2)= ppr-mc3*f3/s3*christoffel(2,2,4)
    christoffel(4,3,3)= mc3*mk*f3/fp/s3
    christoffel(4,3,4)= c3/s3
    christoffel(4,4,3)= c3/s3
    christoffel(4,4,4)= -mc3*s3*f3*christoffel(2,3,3)

end subroutine christoffels
!-------------------------------------------






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
    fluid(5) = ricci
    fluid(6) = areal  * 1d-3   ! changing units from Kpc to Mpc
    fluid(7) = proper * 1d-3   ! changing units from Kpc to Mpc

	
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
! the convention here:
! it is negative
! it is in Kpc !


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
    time = -t
end subroutine look_back_time
!--------------------------------------------------
subroutine time_evolution(szpac,redshift,time)
! time_evolution time is evaluated from the initial redshift
! it provide the amount of time needed to evolve the strcuture 
! from the initial instant (eg. CMB) to the a given redshift 
! the output is in Kpc!
    implicit none
    double precision,  dimension (100), intent(in) :: szpac  
    double precision, intent(in) :: redshift
    double precision, intent(out) :: time
    double precision :: lookback
    call look_back_time(szpac,redshift,lookback)
    time = szpac(10) + lookback
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
subroutine redshift_lcdm(szpac,time,redshift)
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
end subroutine redshift_lcdm
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



