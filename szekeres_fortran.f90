!########################################################################################################
! SzekeresPy ver. 0.4 - Python package for cosmological calculations using the Szekeres Cosmological Model
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
  !  double precision, dimension (npoint) :: szdirection  
    double precision, dimension (100) :: szpac  
    character(len=8), dimension(100) :: szpan
    double precision,  dimension (10) :: fluid
    double precision :: redshift,time
    double precision, parameter :: pi = 4.0d0*datan(1.0d0)
    integer, parameter :: ngrid = 100
   ! double precision, dimension (0:(ngrid-1)) :: temporal, radial, thetal, phial,extral,redshiftal
    logical :: print_names = .True.

   npypac = 15
   allocate(pypac(npypac))
   npyszek = 15
   allocate(pyszek(npyszek))
   szpac = 0.0d0

   call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
   call parameter_names(print_names,szpac,szpan)

   szpoint(1) = 0.0 
   szpoint(2) = 30.0
   szpoint(3) = 0.05*pi
   szpoint(4) = pi 
   
   redshift = 0.5d0
   
   call time_evolution(szpac,redshift,time)
   szpoint(5) = redshift 
   szpoint(6) = time
   szpoint(7) = 1.0
   print *, "time_evolution for redshift=",redshift,'is', time*szpac(25)


   ! testing the evolution of the fluid
   open(101,file="szekeres_fluid_output.txt")
   write(101,*) "# areal,proper,density,expansion,shear,weyl,ricci-3d" 
   do I=1,100
     szpoint(2) = 0.5*I
     call fluid_variables(szpac,npoint,szpoint,fluid)
     write(101,*) fluid(6:7),fluid(1:5)
   enddo

    ! testing null geodesics 
   

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
    ngrid000 = int(input_data(0))
    pypac(1:15)  =  input_data(1:15)
    pyszek(1:15) =  input_data(16:30)
    point(1:7)  =  input_data(31:37)
    redshift = point(5) 
    szpac(100) = 2
    call parameter_values(npypac,pypac,npyszek,pyszek,szpac)
    call age_from_initial(szpac)
    call time_evolution(szpac,redshift,time)
    point(6) = time
    szpac(100) = 5.0
    call fluid_variables(szpac,npoint,point,fluid)
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
    double precision :: redshift,time,age
    double precision, dimension (100) :: szpac  
    double precision,  dimension (10) :: fluid
    double precision, dimension (0:(ngrid-1)), intent(out) :: rho,tht,shr,wey,ric,arl,prp

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
    szpac(100) = 2
    call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    call age_from_initial(szpac)
    call time_evolution(szpac,redshift,time)
    age = szpac(10)
    
    
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(ig,szpac,szpoint,fluid) &
!$OMP SHARED(pypac,pyszek,point,dr,redshift,time,age,rho,tht,shr,wey,ric,arl,prp)
    do ig=0,ngrid-1
        szpoint(1) = point(1)
        szpoint(2) = (ig+1)*dr
        szpoint(3) = point(3)
        szpoint(4) = point(4)
        szpoint(5) = redshift
        szpoint(6) = time
        szpoint(7) = 1.0
        szpac(10)  = age  
        szpac(100) = 5.0
        call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
        call fluid_variables(szpac,npoint,szpoint,fluid)
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
!------------------------------------
subroutine link_cube(input_data,rho,tht,shr,wey,ric,com,prp)
    implicit none
    double precision, dimension(0:37), intent(in) :: input_data
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer, parameter :: ngrid = 100
    integer, parameter :: nbox = 32
    integer, parameter :: nsize = nbox*2 +1
    integer, parameter :: npix = nsize*nsize*nsize 
    integer :: ngrid000, igx,igy,igz,ib
    double precision :: dr,x,y,z,r
    integer :: xshift,yshift,zshift
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point
    double precision, dimension (npoint) :: szpoint    
    double precision :: redshift,time,age
    double precision, dimension (100) :: szpac  
    double precision,  dimension (10) :: fluid
    double precision, dimension (0:npix-1), intent(out) :: rho,tht,shr,wey,ric
    double precision, dimension (0:npix-1,3), intent(out) :: com,prp

    rho = 0.0d0
    tht = 0.0d0
    shr = 0.0d0
    wey = 0.0d0
    ric = 0.0d0
    ric = 0.0d0
    com = 0.0d0
    prp = 0.0d0
    ngrid000 = int(input_data(0))
    pypac(1:15)  =  input_data(1:15)
    pyszek(1:15) =  input_data(16:30)
    point(1:7)  =  input_data(31:37)
    dr = point(2)/(1.0d0*nbox)   
    redshift = point(5) 
    szpac(100) = 2
    call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    call age_from_initial(szpac)
    call time_evolution(szpac,redshift,time)
    age = szpac(10)
    ib = 0

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(igx,igy,igz,ib,x,y,z,xshift,yshift,zshift,r,szpac,szpoint,fluid) &
!$OMP SHARED(pypac,pyszek,point,dr,redshift,time,age,rho,tht,shr,wey,ric,com,prp)
    do igx=-nbox,nbox
        do igy=-nbox,nbox
            do igz=-nbox,nbox
    
        x=dr*igx
        y=dr*igy
        z=dr*igz
        r = sqrt(x*x + y*y + z*z)
        if(r < 0.01) r = 0.01
    
        szpoint(1) = point(1)
        szpoint(2) = r
        szpoint(3) = acos(z/r)
        szpoint(4) = atan2(y,x) 
        szpoint(5) = redshift
        szpoint(6) = time
        szpoint(7) = 1.0
        szpac(10)  = age  
        szpac(100) = 5.0
        call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
        call fluid_variables(szpac,npoint,szpoint,fluid)

        xshift = igx+nbox
        yshift = igy+nbox
        zshift = igz+nbox
        ib = xshift + nsize*(yshift + nsize*zshift)
        rho(ib) = fluid(1)
        tht(ib) = fluid(2)
        shr(ib) = fluid(3)
        wey(ib) = fluid(4)
        ric(ib) = fluid(5)       
       
        com(ib,1) = x
        com(ib,2) = y
        com(ib,3) = z       

! FIX needed: adjusment of the dispacment of the spheres needed
! at this stage, no displacement implementged,
! ie. R is taken to be the areal distance
       prp(ib,1) = fluid(6)*sin(szpoint(3))*cos(szpoint(4))
       prp(ib,2) = fluid(6)*sin(szpoint(3))*sin(szpoint(4))
       prp(ib,3) = fluid(6)*cos(szpoint(3))

!    write(303,*) com(ib,:),rho(ib)
!    write(304,*) prp(ib,:),rho(ib)
   
   

            enddo
        enddo
    enddo   
!$OMP END PARALLEL DO   
    

    
end subroutine link_cube
!------------------------------------
subroutine link_temperature(INTERFACE_FIX_REQUIRED,ND,input_data,temperature,rmax,dmax)
    implicit none
    integer :: INTERFACE_FIX_REQUIRED
    integer, intent(in) :: ND
    double precision, dimension(0:ND-1), intent(in) :: input_data

    double precision, dimension (0:((ND-45)/2)-1), intent(out) :: temperature
    double precision, intent(out) :: rmax,dmax
    
    double precision, dimension (0:((ND-45)/2)-1) :: RA,DEC
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer :: ig,ngrid,N1,N2,imax
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point,direction
    double precision, dimension(npoint)  :: szpoint,szdirection    
    double precision, dimension(0:5) :: tempi
    double precision :: age,tav
    double precision, dimension (100) :: szpac  
    INTERFACE_FIX_REQUIRED = 1 
    szpac(100) = 2
    ngrid = int(input_data(0))

    if(ngrid.ne.((ND-45)/2)) then
        print *, "Number of data: Python->Fortran error [Error 45-T]"
        stop
    endif


    N1 = 1;  N2 = (ND-45)/2
    RA = input_data(N1:N2)
    N1 = N2+1;  N2 = ND-45
    DEC = input_data(N1:N2)
    N1 = N2+1;  N2 = N1+15-1
    pypac(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+15-1
    pyszek(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    point(1:7) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    direction(1:7)  =  input_data(N1:N2)


    call parameter_values(npypac,pypac,npyszek,pyszek,szpac)
    call age_from_initial(szpac)
    age = szpac(10)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(ig,szpac,szpoint,szdirection,tempi) &
!$OMP SHARED(ngrid,RA,DEC,pypac,pyszek,point,age,temperature)
    do ig=0,ngrid-1

        szpoint = point
        szdirection(1) = RA(ig)
        szdirection(2) = DEC(ig)
        szpac(100) = 2
        szpac(10)  = age 

        call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
        call cmb_temperature(szpac,npoint,szpoint,szdirection,tempi)
        temperature(ig) = tempi(0)

    enddo   
!$OMP END PARALLEL DO  

    tav = sum(temperature)/(1d0*(ngrid))
    temperature = 2.725d0*((temperature/tav)- 1.0d0)

    imax= maxloc(temperature,1)
    rmax = RA(imax-1)
    dmax = DEC(imax-1)


   ! call invrotmap(RA,DEC,Lmax,Bmax,lt,bt)



end subroutine link_temperature


!------------------------------------
subroutine link_density(INTERFACE_FIX_REQUIRED,ND,input_data,density,expansion,shear,weyl)
    implicit none
    integer :: INTERFACE_FIX_REQUIRED
    integer, intent(in) :: ND
    double precision, dimension(0:ND-1), intent(in) :: input_data
    double precision, dimension (0:((ND-45)/2)-1), intent(out) :: density,expansion,shear,weyl
    double precision, dimension (0:((ND-45)/2)-1) :: RA,DEC

    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer :: ig,ngrid,N1,N2
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point,direction
    double precision, dimension(npoint)  :: szpoint,szdirection    
    double precision, dimension(0:5) :: tempi
    double precision :: age,dav,zlim

    double precision, dimension (100) :: szpac  
    INTERFACE_FIX_REQUIRED = 1 
    szpac(100) = 2
    ngrid = int(input_data(0))

    if(ngrid.ne.((ND-45)/2)) then
        print *, "Number of data: Python->Fortran error [Error 45-DN]"
        stop
    endif


    N1 = 1;  N2 = (ND-45)/2
    RA = input_data(N1:N2)
    N1 = N2+1;  N2 = ND-45
    DEC = input_data(N1:N2)
    N1 = N2+1;  N2 = N1+15-1
    pypac(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+15-1
    pyszek(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    point(1:7) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    direction(1:7)  =  input_data(N1:N2)


    call parameter_values(npypac,pypac,npyszek,pyszek,szpac)
    call age_from_initial(szpac)
    age = szpac(10)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(ig,szpac,szpoint,szdirection,zlim,tempi) &
!$OMP SHARED(ngrid,RA,DEC,pypac,pyszek,point,age,density,expansion,shear,weyl)
    do ig=0,ngrid-1

        szpoint = point
        szdirection(1) = RA(ig)
        szdirection(2) = DEC(ig)
        szpac(100) = 2
        szpac(10)  = age 
        zlim = point(5)
        call parameter_values(npypac,pypac, npyszek,pyszek, szpac)
        call rho_map(zlim,szpac,npoint,szpoint,szdirection,tempi)

        density(ig) = tempi(1)
        expansion(ig) = tempi(2)
        shear(ig) = tempi(3)
        weyl(ig) = tempi(4)
    enddo   
!$OMP END PARALLEL DO  

    
    dav = sum(density)/(1d0*(ngrid))
    density = ((density/dav)- 1.0d0)




end subroutine link_density


!------------------------------------
subroutine link_null(INTERFACE_FIX_REQUIRED,ND,input_data,temporal,radial,thetal,phial,extral)
    implicit none
    integer :: INTERFACE_FIX_REQUIRED
    integer, intent(in) :: ND
    double precision, dimension(0:ND-1), intent(in) :: input_data
    double precision, dimension (0:ND-45-1), intent(out) :: temporal, radial, thetal, phial,extral
    double precision, dimension(0:ND-45-1)  :: redshiftal
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer :: ngrid,N1,N2
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point,direction

    double precision, dimension (100) :: szpac  
    INTERFACE_FIX_REQUIRED = 1 
    szpac(100) = 2
    ngrid = int(input_data(0))
    if(ngrid.ne.(ND-45)) then
        print *, "Number of data: Python->Fortran error [Error 45-N]"
        stop
    endif
    N1 = 1;  N2 = ND-45
    redshiftal = input_data(N1:N2)
    N1 = N2+1;  N2 = N1+15-1
    pypac(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+15-1
    pyszek(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    point(1:7) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    direction(1:7)  =  input_data(N1:N2)

    call parameter_values(npypac,pypac,npyszek,pyszek,szpac)
    call age_from_initial(szpac)
    call null_geodesic(szpac,npoint,point,direction,ngrid,redshiftal,temporal,radial,thetal,phial,extral)
    
end subroutine link_null

!------------------------------------
subroutine link_distance(INTERFACE_FIX_REQUIRED,ND,input_data,distance)
    implicit none
    integer :: INTERFACE_FIX_REQUIRED
    integer, intent(in) :: ND
    double precision, dimension(0:ND-1), intent(in) :: input_data
    double precision, dimension (0:ND-45-1), intent(out) :: distance
    double precision, dimension(0:ND-45-1)  :: redshiftal
    integer, parameter :: npypac = 15
    integer, parameter :: npyszek = 15
    integer, parameter :: npoint = 7
    integer :: ngrid,N1,N2
    double precision, dimension(npypac)  :: pypac 
    double precision, dimension(npyszek) :: pyszek   
    double precision, dimension(npoint)  :: point,direction

    double precision, dimension (100) :: szpac  
    INTERFACE_FIX_REQUIRED = 1 
    szpac(100) = 2

    ngrid = int(input_data(0))
    if(ngrid.ne.(ND-45)) then
        print *, "Number of data: Python->Fortran error [Error 45-D]"
        stop
    endif
    N1 = 1;  N2 = ND-45
    redshiftal = input_data(N1:N2)
    N1 = N2+1;  N2 = N1+15-1
    pypac(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+15-1
    pyszek(1:15) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    point(1:7) =  input_data(N1:N2)
    N1 = N2+1; N2 = N1+7-1
    direction(1:7)  =  input_data(N1:N2)

    call parameter_values(npypac,pypac,npyszek,pyszek,szpac)
    call age_from_initial(szpac)
    call angular_distance(szpac,npoint,point,direction,ngrid,redshiftal,distance)
    
end subroutine link_distance

!------------------------------------


subroutine rho_map(zlim,szpac,npoint,point,direction,tempi)
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: npoint
    double precision, dimension(npoint) :: point, direction    
    double precision, dimension(0:5), intent(out) :: tempi
    double precision, dimension(0:3)   :: PV,PVi,PVii,NV,NVi,NVii,AV
    double precision, dimension(4,0:3) :: PRK,NRK
    integer :: Ui,I,J,iz
    double precision,dimension(0:10) :: yp,yp1,yp2,yp3
    double precision,dimension(4:5) :: D,Di,Dii
    double precision,dimension(6:10) :: F,Fi,Fii
    double precision :: DA,ADA
    double precision :: xp,xp1,xp2,xp3
    double precision :: ds,dss,rmin
    double precision :: null_test,f16
    double precision :: t,zlim,z

    f16 = 1.0d0/6.0d0
    iz = 0
    Ui = 100000
    dss = 100.d0
    ds = dss
    rmin = 12.0


    ! FIX needed: adaptive step

    call initial_conditions(szpac,npoint,point,direction,PV,NV)

    tempi = 0.0

    F = 0.0
    F(6) = szpac(87)*szpac(23)
    F(7) = szpac(88)*szpac(21)
    F(8) = szpac(88)*szpac(21)
    F(9) = szpac(88)*szpac(23)

    Fi = F; Fii = F

    D(4) = szpac(81)*1d-3
    Di = D; Dii = D
    PVi=PV; PVii=PV
    NVi=NV; NVii=NV
    PRK = 0.0d0
    NRK = 0.0d0

    do I=1,Ui

        if(PV(1).le.rmin) ds = 0.001 + dss*(abs(PV(1)/rmin))
        if(PV(1).ge.rmin) ds = dss

! FIX needed:  safegaurd for orgin-crossing

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(1,J) = NV(J)*ds
        NRK(1,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(1,J)
        NV(J) = NVi(J) + 0.5*NRK(1,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(2,J) = NV(J)*ds
        NRK(2,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(2,J)
        NV(J) = NVi(J) + 0.5*NRK(2,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(3,J) = NV(J)*ds
        NRK(3,J) = AV(J)*ds
        PV(J) = PVi(J) + PRK(3,J)
        NV(J) = NVi(J) + NRK(3,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(4,J) = NV(J)*ds
        NRK(4,J) = AV(J)*ds
        PV(J) = PVi(J) + f16*(PRK(1,J) + 2.0d0*(PRK(2,J) + PRK(3,J)) + PRK(4,J))
        NV(J) = NVi(J) + f16*(NRK(1,J) + 2.0d0*(NRK(2,J) + NRK(3,J)) + NRK(4,J))
        enddo


! FIX needed: adjust the step to the grid
    call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
    D(4) = szpac(81)*1d-3

! FIX needed, the following is for testing only, should be replaced with proper integration
    F(6) =  F(6)  + ds*szpac(87)*szpac(23)
    F(7) =  F(7)  + ds*szpac(88)*szpac(21)
    F(8) =  F(8)  + ds*szpac(88)*szpac(21)
    F(9) =  F(9)  + ds*szpac(88)*szpac(23)
    F(10) = F(10) + ds

    t = PV(0)
    z = abs(NV(0)) - 1.0
        if (z > zlim .and. I>1) then
            yp1(0:3) = NVii; yp1(4:5) = Dii(4:5); yp1(6:10) = Fii(6:10)
            yp2(0:3) = NVi;  yp2(4:5) = Di(4:5);  yp2(6:10) = Fi(6:10)
            yp3(0:3) = NV;   yp3(4:5) = D(4:5);;  yp3(6:10) = F(6:10)

            yp = 0.0d0
            xp1 = dabs(NVii(0)) - 1.0d0
            xp2 = dabs(NVi(0)) - 1.0d0
            xp3 = dabs(NV(0)) - 1.0d0
            xp = zlim
            yp = yp + yp1*((xp - xp2)/(xp1-xp2))*((xp - xp3)/(xp1-xp3))
            yp = yp + yp2*((xp - xp1)/(xp2-xp1))*((xp - xp3)/(xp2-xp3))
            yp = yp + yp3*((xp - xp1)/(xp3-xp1))*((xp - xp2)/(xp3-xp2))
        
                 
            tempi(1) = F(6)/F(10)
            tempi(2) = F(7)/F(10)
            tempi(3) = F(8)/F(10)
            tempi(4) = F(9)/F(10)
            exit
        endif

        PVii = PVi
        PVi  = PV
        NVii = NVi
        NVi  = NV
        Dii = Di
        Di = D

        Fii = Fi
        Fi = F

    enddo

end subroutine rho_map
!----------------------------------------------


subroutine cmb_temperature(szpac,npoint,point,direction,tempi)
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: npoint
    double precision, dimension(npoint) :: point, direction    
    double precision, dimension(0:5), intent(out) :: tempi
    double precision, dimension(0:3)   :: PV,PVi,PVii,NV,NVi,NVii,AV
    double precision, dimension(4,0:3) :: PRK,NRK
    integer :: Ui,I,J,iz
    double precision,dimension(0:10) :: yp,yp1,yp2,yp3
    double precision,dimension(4:5) :: D,Di,Dii
    double precision,dimension(6:10) :: F,Fi,Fii
    double precision :: DA,ADA
    double precision :: xp,xp1,xp2,xp3
    double precision :: ds,dss,rmin
    double precision :: null_test,f16
    double precision :: tlim,t

    f16 = 1.0d0/6.0d0
    iz = 0
    Ui = 100000
    dss = 100.d0
    ds = dss
    rmin = 12.0
    tlim = -2.0d0**18

    ! FIX needed: adaptive step

    call initial_conditions(szpac,npoint,point,direction,PV,NV)

    tempi = 0.0

    F = 0.0
    F(6) = szpac(87)*szpac(23)
    Fi = F; Fii = F

    D(4) = szpac(81)*1d-3
    Di = D; Dii = D
    PVi=PV; PVii=PV
    NVi=NV; NVii=NV
    PRK = 0.0d0
    NRK = 0.0d0

    do I=1,Ui

        if(PV(1).le.rmin) ds = 0.001 + dss*(abs(PV(1)/rmin))
        if(PV(1).ge.rmin) ds = dss

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(1,J) = NV(J)*ds
        NRK(1,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(1,J)
        NV(J) = NVi(J) + 0.5*NRK(1,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(2,J) = NV(J)*ds
        NRK(2,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(2,J)
        NV(J) = NVi(J) + 0.5*NRK(2,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(3,J) = NV(J)*ds
        NRK(3,J) = AV(J)*ds
        PV(J) = PVi(J) + PRK(3,J)
        NV(J) = NVi(J) + NRK(3,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(4,J) = NV(J)*ds
        NRK(4,J) = AV(J)*ds
        PV(J) = PVi(J) + f16*(PRK(1,J) + 2.0d0*(PRK(2,J) + PRK(3,J)) + PRK(4,J))
        NV(J) = NVi(J) + f16*(NRK(1,J) + 2.0d0*(NRK(2,J) + NRK(3,J)) + NRK(4,J))
        enddo


! FIX needed: adjust the step to the grid
    call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
    D(4) = szpac(81)*1d-3

    F(6) =  F(6)  + ds*szpac(87)*szpac(23)
    F(10) = F(10) + ds

    t = PV(0)
        if (t < tlim .and. I>1) then
            yp1(0:3) = NVii; yp1(4:5) = Dii(4:5); yp1(6:10) = Fii(6:10)
            yp2(0:3) = NVi;  yp2(4:5) = Di(4:5);  yp2(6:10) = Fi(6:10)
            yp3(0:3) = NV;   yp3(4:5) = D(4:5);;  yp3(6:10) = F(6:10)

            yp = 0.0d0
            xp1 = PVii(0)
            xp2 = PVi(0)
            xp3 = PV(0)
            xp = tlim 
            yp = yp + yp1*((xp - xp2)/(xp1-xp2))*((xp - xp3)/(xp1-xp3))
            yp = yp + yp2*((xp - xp1)/(xp2-xp1))*((xp - xp3)/(xp2-xp3))
            yp = yp + yp3*((xp - xp1)/(xp3-xp1))*((xp - xp2)/(xp3-xp2))
        
            tempi(0) = 1.0d0/abs(yp(0)) 
            
            tempi(1) = F(6)/F(10)

            exit
        endif

        PVii = PVi
        PVi  = PV
        NVii = NVi
        NVi  = NV
        Dii = Di
        Di = D

        Fii = Fi
        Fi = F

    enddo

end subroutine cmb_temperature
!----------------------------------------------

subroutine angular_distance(szpac,npoint,point,direction,ngrid,redshiftal,distance)
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: ngrid, npoint
    double precision, dimension(npoint) :: point, direction    
    double precision, dimension(0:(ngrid-1)), intent(in) :: redshiftal
    double precision, dimension(0:(ngrid-1)), intent(out) :: distance
    double precision, dimension(0:3)   :: PV,PVi,PVii,NV,NVi,NVii,AV
    double precision, dimension(4,0:3) :: PRK,NRK
    double precision, dimension(4) :: DRK,VRK
    integer :: Ui,I,J,iz
    double precision,dimension(0:5) :: yp,yp1,yp2,yp3
    double precision,dimension(4:5) :: D,Di,Dii
    double precision :: DA,DAi,VDA,VDAi,ADA
    double precision :: xp,xp1,xp2,xp3
    double precision :: ds,dss,z,rmin
    double precision :: null_test,f16

 
    f16 = 1.0d0/6.0d0
    iz = 0
    Ui = 100000
    dss = 150.d0
    ds = dss
    rmin = 5.0
! FIX needed: adaptive step

    call initial_conditions(szpac,npoint,point,direction,PV,NV)

    DA=0d0
    VDA=1d0
    ADA=0d0


    D(4) = szpac(81)*1d-3; D(5) = 0.0d0

    DAi = DA; VDAi = VDA 
    Di = D; Dii = D
    PVi=PV; PVii=PV
    NVi=NV; NVii=NV
    PRK = 0.0d0; NRK = 0.0d0
    DRK = 0.0d0; VRK = 0.0d0; 

    do I=1,Ui

        if(PV(1).le.rmin) ds = 0.001 + dss*(PV(1)/rmin)
        if(PV(1).ge.rmin) ds = dss

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(1,J) = NV(J)*ds
        NRK(1,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(1,J)
        NV(J) = NVi(J) + 0.5*NRK(1,J)
        enddo
        VRK(1) = ADA*ds
        DRK(1) = VDA*ds
        VDA = VDAi+0.5d0*VRK(1)
        DA  = DAi +0.5d0*DRK(1)


        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(2,J) = NV(J)*ds
        NRK(2,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(2,J)
        NV(J) = NVi(J) + 0.5*NRK(2,J)
        enddo
        VRK(2) = ADA*ds
        DRK(2) = VDA*ds
        VDA = VDAi+0.5d0*VRK(2)
        DA  = DAi +0.5d0*DRK(2)

  

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(3,J) = NV(J)*ds
        NRK(3,J) = AV(J)*ds
        PV(J) = PVi(J) + PRK(3,J)
        NV(J) = NVi(J) + NRK(3,J)
        enddo
        VRK(3) = ADA*ds
        DRK(3) = VDA*ds
        VDA = VDAi + VRK(3)
        DA  = DAi  + DRK(3)

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(4,J) = NV(J)*ds
        NRK(4,J) = AV(J)*ds
        PV(J) = PVi(J) + f16*(PRK(1,J) + 2.0d0*(PRK(2,J) + PRK(3,J)) + PRK(4,J))
        NV(J) = NVi(J) + f16*(NRK(1,J) + 2.0d0*(NRK(2,J) + NRK(3,J)) + NRK(4,J))
        enddo
        VRK(4) = ADA*ds
        DRK(4) = VDA*ds
        VDA = VDAi + f16*(VRK(1)+2.0d0*(VRK(2) + VRK(3))+VRK(4))
        DA  = DAi  + f16*(DRK(1)+2.0d0*(DRK(2) + DRK(3))+DRK(4))

   
! FIX needed: adjust the step to the grid
    call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
    D(4) = szpac(81)*1d-3
    D(5) = DA

    z = abs(NV(0)) - 1.0
        if (z > redshiftal(iz) .and. I>1) then
            yp1(0:3) = PVii; yp1(4:5) = Dii(4:5)
            yp2(0:3) = PVi; yp2(4:5) = Di(4:5)
            yp3(0:3) = PV; yp3(4:5) = D(4:5)

            yp = 0.0d0
            xp1 = dabs(NVii(0)) - 1.0d0
            xp2 = dabs(NVi(0)) - 1.0d0
            xp3 = dabs(NV(0)) - 1.0d0
            xp = redshiftal(iz)
            yp = yp + yp1*((xp - xp2)/(xp1-xp2))*((xp - xp3)/(xp1-xp3))
            yp = yp + yp2*((xp - xp1)/(xp2-xp1))*((xp - xp3)/(xp2-xp3))
            yp = yp + yp3*((xp - xp1)/(xp3-xp1))*((xp - xp2)/(xp3-xp2))
            distance(iz) = yp(5) * 1d-3
            iz=iz+1
            if(iz == ngrid) exit
        endif

        PVii = PVi
        PVi  = PV
        NVii = NVi
        NVi  = NV
        Dii = Di
        Di = D

        VDAi    = VDA
        DAi     = DA

  
    enddo

end subroutine angular_distance
!------------------------------------


subroutine null_geodesic(szpac,npoint,point,direction,ngrid,redshiftal,temporal,radial,thetal,phial,extral)
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: ngrid, npoint
    double precision, dimension(npoint) :: point, direction    
    double precision, dimension(0:(ngrid-1)), intent(in) :: redshiftal
    double precision, dimension(0:(ngrid-1)), intent(out) :: temporal, radial, thetal, phial, extral
    double precision, dimension(0:3)   :: PV,PVi,PVii,NV,NVi,NVii,AV
    double precision, dimension(4,0:3) :: PRK,NRK
    integer :: Ui,I,J,iz
    double precision,dimension(0:5) :: yp,yp1,yp2,yp3
    double precision,dimension(4:5) :: D,Di,Dii
    double precision :: DA,ADA
    double precision :: xp,xp1,xp2,xp3
    double precision :: ds,dss,z,rmin
    double precision :: null_test,f16


    f16 = 1.0d0/6.0d0
    iz = 0
    Ui = 100000
    dss = 100.d0
    ds = dss
    rmin = 5.0

    ! FIX needed: adaptive step

    call initial_conditions(szpac,npoint,point,direction,PV,NV)
    D(4) = szpac(81)*1d-3
    Di = D; Dii = D
    PVi=PV; PVii=PV
    NVi=NV; NVii=NV
    PRK = 0.0d0
    NRK = 0.0d0

    do I=1,Ui

        if(PV(1).le.rmin) ds = 0.001 + dss*(PV(1)/rmin)
        if(PV(1).ge.rmin) ds = dss

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(1,J) = NV(J)*ds
        NRK(1,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(1,J)
        NV(J) = NVi(J) + 0.5*NRK(1,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(2,J) = NV(J)*ds
        NRK(2,J) = AV(J)*ds
        PV(J) = PVi(J) + 0.5*PRK(2,J)
        NV(J) = NVi(J) + 0.5*NRK(2,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(3,J) = NV(J)*ds
        NRK(3,J) = AV(J)*ds
        PV(J) = PVi(J) + PRK(3,J)
        NV(J) = NVi(J) + NRK(3,J)
        enddo

        call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
        do J=0,3
        PRK(4,J) = NV(J)*ds
        NRK(4,J) = AV(J)*ds
        PV(J) = PVi(J) + f16*(PRK(1,J) + 2.0d0*(PRK(2,J) + PRK(3,J)) + PRK(4,J))
        NV(J) = NVi(J) + f16*(NRK(1,J) + 2.0d0*(NRK(2,J) + NRK(3,J)) + NRK(4,J))
        enddo


! FIX needed: adjust the step to the grid
    call light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
    D(4) = szpac(81)*1d-3

    z = abs(NV(0)) - 1.0
        if (z > redshiftal(iz) .and. I>1) then
            yp1(0:3) = PVii; yp1(4:5) = Dii(4:5)
            yp2(0:3) = PVi; yp2(4:5) = Di(4:5)
            yp3(0:3) = PV; yp3(4:5) = D(4:5)

            yp = 0.0d0
            xp1 = dabs(NVii(0)) - 1.0d0
            xp2 = dabs(NVi(0)) - 1.0d0
            xp3 = dabs(NV(0)) - 1.0d0
            xp = redshiftal(iz)
            yp = yp + yp1*((xp - xp2)/(xp1-xp2))*((xp - xp3)/(xp1-xp3))
            yp = yp + yp2*((xp - xp1)/(xp2-xp1))*((xp - xp3)/(xp2-xp3))
            yp = yp + yp3*((xp - xp1)/(xp3-xp1))*((xp - xp2)/(xp3-xp2))
            temporal(iz) = yp(0)*szpac(25)*1d-3
            radial(iz)   = yp(4)  ! alternativel use: yp(1)  
            phial(iz)    = yp(2)
            thetal(iz)   = yp(3)
            extral(iz) = null_test
            iz=iz+1
            if(iz == ngrid) exit
        endif

! FIX needed: further tests/diagnostics

        PVii = PVi
        PVi  = PV
        NVii = NVi
        NVi  = NV
        Dii = Di
        Di = D

    enddo

end subroutine null_geodesic
!------------------------------------
subroutine light_propagation(szpac,PV,NV,AV,DA,ADA,null_test)
    implicit none
    double precision,  intent(inout), dimension (100) :: szpac 
    double precision, dimension(0:3), intent(in) :: PV,NV
    double precision, dimension(0:3), intent(out) :: AV 
    double precision, dimension(0:3,0:3)   :: GIJ
    double precision, dimension(0:3,0:3,0:3) :: GAM
    integer :: G,I,J
    double precision :: DA,ADA
    double precision, intent(out) :: null_test

    call christoffels2(szpac,PV,GIJ,GAM)

    do G=0,3
        AV(G) = 0.0
        do I=0,3
            do J=0,3
                AV(G) = AV(G) - GAM(G,I,J)*NV(I)*NV(J)
            enddo
        enddo
    enddo


    null_test = 0.0
    do I=0,3
        do J=0,3
            null_test = null_test + GIJ(I,J)*NV(I)*NV(J)
        enddo
    enddo   

  
    ADA = -5d-1*szpac(87)*(NV(0)*NV(0))*DA
  

end subroutine light_propagation
!------------------------------------
subroutine initial_conditions(szpac,npoint,point,direction,PV,NV)
    implicit none
    double precision,  intent(inout), dimension (100) :: szpac 
    integer, intent(in) :: npoint    
    double precision, dimension (npoint) :: point, direction
    double precision, intent(out), dimension(0:3) :: PV,NV
    double precision, dimension(0:3,0:3)   :: GIJ
    double precision, dimension(0:3,0:3,0:3) :: GAM

    PV(0) = 0.0d0
    PV(1) = point(2)
    PV(2) = point(3)
    PV(3) = point(4)

    call christoffels2(szpac,PV,GIJ,GAM)
    call null_initial(szpac,npoint,point,direction,NV)

end subroutine initial_conditions
!------------------------------------
subroutine null_initial(szpac,npoint,point,direction,NV)
    implicit none
    double precision,  intent(in), dimension (100) :: szpac 
    integer, intent(in) :: npoint    
    double precision, dimension (npoint), intent(in) :: point, direction 
    double precision, dimension(0:3,0:3) :: E
    double precision, dimension(0:3) :: LV
    double precision, dimension(0:3), intent(out) :: NV
    integer :: I,J
    double precision :: RA,DEC,th,ph,epsilon_ang,epsilon_dir
    double precision :: aR,aRi,aRr
    double precision :: sth,cth,sph,cph,sthi,n,nd,mcth
    double precision :: mk,mki,mkis,sna,snc
    double precision :: k,qr,pr,s,sr,si,W,Wi,N1,N2,N3

    epsilon_dir = 1d-8; epsilon_ang = 1d-10
    717 continue
    RA = direction(1)*szpac(33) + epsilon_ang
    DEC = direction(2)*szpac(33) + epsilon_ang
    th = point(3)
    ph = point(4)
    aRi   = szpac(80)
    aR    = szpac(81)
    aRr   = szpac(82)
    k    = szpac(64)
    qr   = szpac(68) 
    pr   = szpac(71) 
    s    = szpac(73) 
    sr   = szpac(74) 
    si = 1.0/s
    mk = 1.0 - k
    mki = 1.0/mk
    mkis = sqrt(mki)
    sth = sin(th)
    cth = cos(th)
    sph = sin(ph)
    cph = cos(ph)
    sthi = 1.0d0/sth
    mcth = 1.0d0 - cth
    n   =  pr*cph + qr*sph
    nd  = -pr*sph + qr*cph

    sna = n*si*sth+sr*si*cth 
    snc = n*si*mcth+sr*si*sth
    W = mkis*(aRr+aR*sna)
    Wi = 1.0d0/W

! Null vector in the Lorentz's frame:
    LV(0) = -1.0d0
    LV(1) = sin(DEC)
    LV(2) = cos(DEC)*cos(RA)
    LV(3) = cos(DEC)*sin(RA)

! vierbein:  e_A^alpha = E(A,alpha)
    E(0,0) = 1.0
    E(0,1) = 0.0
    E(0,2) = 0.0
    E(0,3) = 0.0    
    E(1,0) = 0.0
    E(1,1) = Wi
    E(1,2) = (sr*sth*si+n*mcth*si)*Wi
    E(1,3) = - nd*mcth*si*Wi*sthi
    E(2,0) = 0.0
    E(2,1) = 0.0
    E(2,2) = aRi
    E(2,3) = 0.0  
    E(3,0) = 0.0
    E(3,1) = 0.0
    E(3,2) = 0.0
    E(3,3) = aRi*sthi

! null vector in the Szekeres coordinates
    do I=0,3
        NV(I) = 0.0
        do J=0,3
            NV(I) = NV(I) + LV(J)*E(J,I)
        enddo
    enddo

    N1 = NV(1)
    N2 = NV(2)*aR
    N3 = NV(3)*aR*sth

    if (N1.le.epsilon_dir.and.(abs(N2).le.epsilon_dir.and.abs(N3).le.epsilon_dir)) then
        epsilon_ang = epsilon_ang + 1d-3
        if(epsilon_ang>2d-2) then
            print *, "Initial condition problem, terminting [Error 717]"
            stop
        endif
    goto 717
    endif


end subroutine null_initial
    

!------------------------------------

subroutine christoffels2(szpac,PV,GIJ,GAM)
    implicit none

    double precision,  intent(inout), dimension (100)     :: szpac 
    double precision, dimension(0:3), intent(in)          :: PV
    double precision, dimension(0:3,0:3), intent(out)     :: GIJ
    double precision, dimension(0:3,0:3,0:3), intent(out) :: GAM
    double precision, dimension(0:3,0:3,0:3) :: GDER
    double precision, dimension(0:3,0:3) :: GINV
    double precision :: th,ph
    double precision :: aR,aRt, aRr, aRrr, aRtr,aRtrr,aRi,aR2,rs1,rs2
    double precision :: sth,cth,sph,cph,sthi,sth2,n,nd,nr,ndr,mcth,mcth2
    double precision :: si,si2,si3,mk,mki,mki2,mksi,mns2,mnss,mns,nsc,nd2
    double precision :: m,mr,mrr,k,kr,krr,q,qr,qrr,p,pr,prr,s,sr,srr,epr
    double precision :: GW,GWi, g11_part_1, g11_part_2,denumerator,f13
    integer :: I,J,G,H

    call areal_evolution(szpac,PV)
    f13 = 1.0d0/3.0d0
    aRi   = szpac(80)
    aR    = szpac(81)
    aRr   = szpac(82)
    aRrr  = szpac(83)
    aRt   = szpac(84)
    aRtr  = szpac(85)
    aRtrr = szpac(86)
    aR2 = aR*aR

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
    si = 1.0/s
    si2 = si*si
    si3 = si2*si
    mk = 1.0 - k
    mki = 1.0/mk
    mksi = mk*si
    mki2 = mki*mki
    rs1 = aR*si
    rs2 = rs1*rs1
    th = PV(2)
    ph = PV(3)
    sth = sin(th)
    cth = cos(th)
    sph = sin(ph)
    cph = cos(ph)
    sthi = 1.0d0/sth
    sth2 = sth*sth
    mcth = 1.0d0 - cth
    mcth2 = mcth*mcth
    n   =  pr*cph + qr*sph
    nr  =  prr*cph + qrr*sph
    nd  = -pr*sph + qr*cph
    ndr = -prr*sph + qrr*cph
    nd2 = nd*nd
    nsc = n*sth+sr*cth
    mns = mcth*n + sr*sth
    mnss = mns*s  
    mns2 = mns*mns
    epr = nsc*si

    denumerator = 1.0d0/(aRr - aR*epr)
    szpac(87) = 2.0d0*(mr-3.0d0*m*epr)*aRi*aRi*denumerator
    szpac(88) = (aRtr + 2*aRt*aRr*aRi - 3.0d0*aRt*epr)*denumerator*f13
    szpac(89) = -(aRtr - aRt*aRr*aRi)*f13*denumerator
    szpac(90) = m*(3*aRr - aR*mr*m)*f13*aRi*aRi*aRi*denumerator

! FIX required: above should be provided by a separate function/subroutine 

    g11_part_1 = -mki*(  (aRr + rs1*(sr*cth + n*sth) )**2    )
    g11_part_2 =  -rs2*( (sr*sth + n*mcth)**2 + (nd*mcth)**2  )

    GIJ = 0.0d0
    GIJ(0,0) = 1.0d0
    GIJ(1,1) = g11_part_1 + g11_part_2
    GIJ(1,2) = rs2*(s*sr*sth + s*n*mcth )
    GIJ(2,1) = GIJ(1,2)
    GIJ(1,3) = -rs2*s*nd*sth*mcth
    GIJ(3,1) = GIJ(1,3)
    GIJ(2,2) = -aR*aR
    GIJ(3,3) = GIJ(2,2)*sth*sth
    
    GW = GIJ(1,1)*GIJ(2,2)*GIJ(3,3) - GIJ(1,2)*GIJ(2,1)*GIJ(3,3) - GIJ(1,3)*GIJ(2,2)*GIJ(3,1)
    GWi = 1.0d0/GW
    GINV = 0.0d0
    GINV(0,0) = 1.0d0
    GINV(0,1) = 0.0d0
    GINV(0,2) = 0.0d0
    GINV(0,3) = 0.0d0
    GINV(1,0) = 0.0d0
    GINV(1,1) = GIJ(2,2)*GIJ(3,3)*GWi
    GINV(1,2) = -GIJ(1,2)*GIJ(3,3)*GWi
    GINV(1,3) = -GIJ(1,3)*GIJ(2,2)*GWi
    GINV(1,0) = 0.0d0
    GINV(2,1) = -GIJ(2,1)*GIJ(3,3)*GWi
    GINV(2,2) = (GIJ(1,1)*GIJ(3,3) - GIJ(1,3)*GIJ(3,1))*GWi
    GINV(2,3) = GIJ(1,3)*GIJ(2,1)*GWi
    GINV(3,0) = 0.0d0
    GINV(3,1) = -GIJ(2,2)*GIJ(3,1)*GWi
    GINV(3,2) =  GIJ(1,2)*GIJ(3,1)*GWi
    GINV(3,3) = (GIJ(1,1)*GIJ(2,2) - GIJ(1,2)*GIJ(2,1))*GWi

    GDER(0,0,0) = 0.0d0 
    GDER(0,0,1) = 0.0d0 
    GDER(0,0,2) = 0.0d0 
    GDER(0,0,3) = 0.0d0 
    GDER(0,1,0) = 0.0d0 
    GDER(0,1,1) = 0.0d0 
    GDER(0,1,2) = 0.0d0 
    GDER(0,1,3) = 0.0d0 
    GDER(0,2,0) = 0.0d0 
    GDER(0,2,1) = 0.0d0 
    GDER(0,2,2) = 0.0d0 
    GDER(0,2,3) = 0.0d0 
    GDER(0,3,0) = 0.0d0 
    GDER(0,3,1) = 0.0d0 
    GDER(0,3,2) = 0.0d0 
    GDER(0,3,3) = 0.0d0 
    GDER(1,0,0) = 0.0d0 
    GDER(1,0,1) = 0.0d0 
    GDER(1,0,2) = 0.0d0 
    GDER(1,0,3) = 0.0d0 
    GDER(1,1,0) = -2*mcth2*nd2*aR*aRt*si2- 2*mns2*aR*aRt*si2- (nsc*aR*si+ aRr)*(2*nsc*aRt*si+ 2*aRtr)*mki 
    GDER(1,1,1) = 2*mcth2*nd2*aR2*sr*si3- 2*mcth2*nd2*aR*aRr*si2- mcth2*nd*2.0*ndr*aR2*si2+ 2*mns2*aR2*sr*si3
    GDER(1,1,1)=GDER(1,1,1)-2*mns2*aR*aRr*si2-mns*(2*mcth*nr + 2*sth*srr)*aR2*si2
    GDER(1,1,1)=GDER(1,1,1)-(nsc*aR*si+ aRr)*(-2*nsc*aR*sr*si2+2*nsc*aRr*si+2*(nr*sth+cth*srr)*aR*si+2*aRrr)*mki
    GDER(1,1,1)=GDER(1,1,1)-(nsc*aR*si+ aRr)**2*kr*mki2 
    GDER(1,1,2) = -2*mcth*nd2*aR2*sth*si2-mns*(2*n*sth+2*sr*cth)*aR2*si2-2*(n*cth-sr*sth)*(nsc*aR*si+aRr)*aR*mksi
    GDER(1,1,3) = mcth2*nd*2.0*n*aR2*si2- 2*mcth*mns*nd*aR2*si2- 2*nd*(nsc*aR*si+ aRr)*aR*sth*mksi
    GDER(1,2,0) = 2*mnss*aR*aRt*si2
    GDER(1,2,1) = -2*mnss*aR2*sr*si3+ 2*mnss*aR*aRr*si2+ (mcth*n*sr + mcth*nr*s + s*sth*srr + sr*sth*sr)*aR2*si2
    GDER(1,2,2) = (n*s*sth + s*sr*cth)*aR2*si2
    GDER(1,2,3) = mcth*nd*aR2*si
    GDER(1,3,0) = -2*mcth*nd*aR*sth*aRt*si
    GDER(1,3,1) = mcth*nd*aR2*sth*sr*si2- 2*mcth*nd*aR*sth*aRr*si- mcth*ndr*aR2*sth*si
    GDER(1,3,2) = -mcth*nd*aR2*cth*si- nd*aR2*sth2*si
    GDER(1,3,3) = mcth*n*aR2*sth*si
    GDER(2,0,0) = 0.0d0 
    GDER(2,0,1) = 0.0d0 
    GDER(2,0,2) = 0.0d0 
    GDER(2,0,3) = 0.0d0 
    GDER(2,1,0) = 2*mnss*aR*aRt*si2
    GDER(2,1,1) = -2*mnss*aR2*sr*si3+ 2*mnss*aR*aRr*si2+ (mcth*n*sr + mcth*nr*s + s*sth*srr + sr*sth*sr)*aR2*si2
    GDER(2,1,2) = (n*s*sth + s*sr*cth)*aR2*si2
    GDER(2,1,3) = mcth*nd*aR2*si
    GDER(2,2,0) = -2*aR*aRt 
    GDER(2,2,1) = -2*aR*aRr 
    GDER(2,2,2) = 0.0d0 
    GDER(2,2,3) = 0.0d0 
    GDER(2,3,0) = 0.0d0 
    GDER(2,3,1) = 0.0d0 
    GDER(2,3,2) = 0.0d0 
    GDER(2,3,3) = 0.0d0 
    GDER(3,0,0) = 0.0d0 
    GDER(3,0,1) = 0.0d0 
    GDER(3,0,2) = 0.0d0 
    GDER(3,0,3) = 0.0d0 
    GDER(3,1,0) = -2*mcth*nd*aR*sth*aRt*si
    GDER(3,1,1) = mcth*nd*aR2*sth*sr*si2- 2*mcth*nd*aR*sth*aRr*si- mcth*ndr*aR2*sth*si
    GDER(3,1,2) = -mcth*nd*aR2*cth*si- nd*aR2*sth2*si
    GDER(3,1,3) = mcth*n*aR2*sth*si
    GDER(3,2,0) = 0.0d0 
    GDER(3,2,1) = 0.0d0 
    GDER(3,2,2) = 0.0d0 
    GDER(3,2,3) = 0.0d0 
    GDER(3,3,0) = -2*aR*sth2*aRt 
    GDER(3,3,1) = -2*aR*sth2*aRr 
    GDER(3,3,2) = -2*aR2*sth*cth 
    GDER(3,3,3) = 0.0d0 

    do I=0,3
        do J=0,3
            do G=0,3

    GAM(I,J,G) = 0.0
    do H=0,3
    GAM(I,J,G)=GAM(I,J,G)+0.5d0*GINV(I,H)*(GDER(H,G,J)+GDER(J,H,G)-GDER(J,G,H))
    enddo

            enddo
        enddo
    enddo


end subroutine christoffels2
!-------------------------------------------
subroutine fluid_variables(szpac,npoint,point,fluid)
! calulates the fluid variables: 
! density (rho), expansion (tht), shear (shr), Weyl (wey), 3D Ricci (ric)
! based on eq. (3.2), (3.3), (3.4), (3.5), and (3.7) of arxiv:1704.02810
    implicit none
    double precision,  dimension (100) :: szpac  
    integer, intent(in) :: npoint
    double precision,  dimension (npoint) :: point
    double precision,  dimension (10), intent(out) :: fluid
    double precision :: t,r,theta,phi
    double precision :: m,mi,mr,mrr,k,ki,kr,krr,q,qr,qrr,p,pr,prr,s,sr,srr,si,n,epr
    double precision :: redshift, zz2, zz3,zz4,ez,gzn,hzn
    double precision :: f12,f13,f16,f19    
    double precision :: rho1,rho2,tht1,tht2,shr1,shr2,wey1,wey2,ric1,ric2
    double precision :: density, expansion, shear, weyl, ricci,areal,proper
    double precision :: aRi,aRi2,aRi3,denumerator
    double precision :: aR, aRt, aRr, aRrr, aRtr, aRtrr, rpR  

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

    call areal_radius(szpac,npoint,point,aR,aRt, aRr, aRrr, aRtr, aRtrr)
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

subroutine areal_radius(szpac,npoint,point,aR,aRt, aRr, aRrr, aRtr, aRtrr)
! calulates areal distance R and its derivatives  
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
    double precision :: aR,aRt, aRr, aRrr, aRtr, aRtrr, aRi
    double precision :: q,qr,qrr,p,pr,prr,s,sr,srr
    double precision :: m,mr,mrr,k,kr,krr,l
    double precision :: kr1,kr2,kr3,kr4
    double precision :: krr1,krr2,krr3,krr4
    double precision :: krrr1,krrr2,krrr3,krrr4
    double precision :: r,rp,rt
    double precision :: rr,rrp,rtr
    double precision :: rrr,rrrp,rtrr

! FIX needed: time to point(6) only

    if(point(7)<0.5) then   
       time = point(1)
    else
       time = point(6)
    endif
    

! FIX needed: adaptive step
    dti = 2048.0
    Ni = int( abs(time/dti) )
    if(Ni.le.3) Ni = 3

    Ni = 128
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
    
! FIX needed: rt,rtr,rtrr -> into array
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

! FIX needed: adaptive step interpolation 

    rt = sqrt(2.0d0*(m/rp) - k + (l/3.0d0)*rp*rp)
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
    aRi   = 1.0d0/r

    szpac(80) = aRi 
    szpac(81) = aR
    szpac(82) = aRr
    szpac(83) = aRrr
    szpac(84) = aRt
    szpac(85) = aRtr
    szpac(86) = aRtrr




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
    Ni = 128
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
subroutine areal_evolution(szpac,PV)
    implicit none
        
    double precision,  dimension (100) :: szpac  
    double precision, dimension(0:3), intent(in) :: PV

    integer :: I,Ni
    double precision :: dti,dt,time
    double precision :: aR,aRt, aRr, aRrr, aRtr, aRtrr, aRi
    double precision :: q,qr,qrr,p,pr,prr,s,sr,srr
    double precision :: m,mr,mrr,k,kr,krr,l
    double precision :: kr1,kr2,kr3,kr4,f16,f13
    double precision :: krr1,krr2,krr3,krr4
    double precision :: krrr1,krrr2,krrr3,krrr4
    double precision :: r,rp,rt
    double precision :: rr,rrp,rtr
    double precision :: rrr,rrrp,rtrr


    time = szpac(10)  + PV(0)

! FIX needed: adaptive step

    dti = 1024.0
    Ni = int( abs(time/dti) )
    if(Ni.le.3) Ni = 3
    l = szpac(26)

    f13 = 1.0d0/3.0d0
    f16 = 0.5d0*f13

    Ni = 256

    dt = time/(Ni*1.0d0)
	

    r = PV(1)
    rr = 1.0d0
    rrr = 0.0d0
    
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
    
! FIX needed: rt,rtr,rtrr -> into array (for adaptive step)

    do I=1,Ni
      rp = r
      rrp = rr
      rrrp = rrr
      rt = dsqrt(2.0d0*(m/rp) - k + (l*f13)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l*f13)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l*f13)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr1   = dt*rt
      krr1  = dt*rtr
      krrr1 = dt*rtrr
      rp   = r   + kr1*5d-1
      rrp  = rr  + krr1*5d-1
      rrrp = rrr + krrr1*5d-1
      rt = dsqrt(2.0d0*(m/rp) - k + (l*f13)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l*f13)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l*f13)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr2   = dt*rt
      krr2  = dt*rtr
      krrr2 = dt*rtrr
      rp   = r   + kr2*5d-1
      rrp  = rr  + krr2*5d-1
      rrrp = rrr + krrr2*5d-1
      rt = dsqrt(2.0d0*(m/rp) - k + (l*f13)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l*f13)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l*f13)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr3   = dt*rt
      krr3  = dt*rtr
      krrr3 = dt*rtrr
      rp   = r   + kr3
      rrp  = rr  + krr3
      rrrp = rrr + krrr3
      rt = dsqrt(2.0d0*(m/rp) - k + (l*f13)*rp*rp)
      rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l*f13)*rp*rrp)/rt
      rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
      rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
      rtrr = rtrr + (l*f13)*(rrp**2 + rp*rrrp) - rtr*rtr
      rtrr = rtrr/rt
      kr4   = dt*rt
      krr4  = dt*rtr
      krrr4 = dt*rtrr
      rp   = r   + (kr1  +2.0d0*kr2  +2.0d0*kr3  +kr4)*f16 
      rrp  = rr  + (krr1 +2.0d0*krr2 +2.0d0*krr3 +krr4)*f16
      rrrp = rrr + (krrr1+2.0d0*krrr2+2.0d0*krrr3+krrr4)*f16
      r   = rp
      rr  = rrp
      rrr = rrrp
    enddo

! FIX needed: adaptive step interpolation 

    rt = sqrt(2.0d0*(m/rp) - k + (l*f13)*rp*rp)
    rtr = ((mr/rp) - (m/(rp**2))*rrp - 0.5d0*kr +  (l*f13)*rp*rrp)/rt
    rtrr = (mrr/rp)  - 2.0d0*(mr/(rp**2))*rrp - (m/(rp**2))*rrrp
    rtrr = rtrr + 2.0d0*(m/(rp**3))*(rrp**2) - 0.5d0*krr 
    rtrr = rtrr + (l*f13)*(rrp**2 + rp*rrrp) - rtr*rtr
    rtrr = rtrr/rt
    aR    = r
    aRt   = rt
    aRr   = rr
    aRrr  = rrr
    aRtr  = rtr
    aRtrr = rtrr
    aRi   = 1.0d0/r

    szpac(80) = aRi 
    szpac(81) = aR
    szpac(82) = aRr
    szpac(83) = aRrr
    szpac(84) = aRt
    szpac(85) = aRtr
    szpac(86) = aRtrr




end subroutine areal_evolution
!--------------------------------------------------

subroutine look_back_time(szpac,redshift,time)
! calculates the lookback time, BUT....
! the convention here:
! it is negative, it is in Kpc !
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
subroutine age_from_initial(szpac)
    implicit none
    integer :: I,Ni
    double precision,  dimension (100) :: szpac  
    double precision :: omega_matter,omega_lambda,omega_curvature,omega_radiation
    double precision :: Ho, Haa,da,ai,a,ia1,ia2,ia3,ia4
    double precision :: z_initial,t_rk1,t_rk2,t_rk3,t_rk4,ti,t,f16
    Ho = szpac(11)
    omega_matter =  szpac(3)
    omega_lambda =  szpac(2)
    omega_radiation = szpac(7)
    omega_curvature = szpac(9)
    z_initial = szpac(30) 
    Ni = 2048
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
    szpac(10) = t 
    
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
    double precision :: r0,dlr,Ak,Am,amp,alpha

    amp = szpac(51)
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
!--------------------------------------------------
subroutine invrotmap(ra,dec,Lmaxi,Bmaxi,l,b)

    implicit none
    double precision ra,dec,pi,pi180,Lmaxi,Bmaxi
    double precision Lmax,Bmax,lcmb,bcmb,x,y,z,DV(4),DA(4)
    double precision ex,ey,ez
    double precision l,b,rai,deci

    pi = 4d0*datan(1d0)
    pi180 = pi/180.0d0
    rai = (ra)*pi180
    deci = dec*pi180

    DV(2) = dsin(deci)
    DV(3) = dcos(rai)*dcos(deci)
    DV(4) = dsin(rai)*dcos(deci)
	
    Lmax = Lmaxi*pi180
    Bmax = Bmaxi*pi180
    
    lcmb  = 276.4d0*pi180
    bcmb = 29.3d0*pi180
        x = DV(2)
        y = DV(3)
        z = DV(4)

! Rotate OZ(-Lmax) then R-OY(+Bmax) then R-OY(-bcmb) then R-OZ(+lcmb):
! ANTI-Rotate: R-OZ(-lcmb), then R-OY(+bcmb), then R-OY(-Bmax), then OZ(+Lmax)
    call rotxyz(3,-lcmb,DV,DA)
        DV=DA
    call rotxyz(2,+bcmb,DV,DA)
        DV=DA
    call rotxyz(2,-Bmax,DV,DA)
        DV=DA
    call rotxyz(3,+Lmax,DV,DA)

    ez = DA(4)
    ey = DA(3)
    ex = DA(2)

    if(ez.ge.0d0) b = dasin(dabs(ez))
    if(ez.le.0d0) b = -dasin(dabs(ez)) 
    if(ey.ge.0d0 .and. ex.ge.0d0) l = datan(dabs(ey/ex))
    if(ey.ge.0d0 .and. ex.le.0d0) l = pi - datan(dabs(ey/ex))
    if(ey.le.0d0 .and. ex.le.0d0) l = pi + datan(dabs(ey/ex))
    if(ey.le.0d0 .and. ex.ge.0d0) l = 2*pi-(datan(dabs(ey/ex)))
    if(ex.eq.0d0 .and. ey.gt.0d0) l = 0.5d0*pi
    if(ex.eq.0d0 .and. ey.lt.0d0) l = 0.5d0*pi
    if(ex.eq.0d0 .and. ey.eq.0d0) l = 0.0d0*pi
    l=l/pi180
    b=b/pi180

end subroutine invrotmap
!--------------------------------------------
subroutine rotxyz(n,theta,ri,ro)
    implicit none
    integer n
    double precision ri(4),ro(4)
    double precision xo,yo,zo,xi,yi,zi,theta
        xi = ri(2)
        yi = ri(3)
        zi = ri(4)
        if(n.eq.1) then
            xo = xi
            yo = dcos(theta)*yi - dsin(theta)*zi
            zo = dsin(theta)*yi + dcos(theta)*zi
        endif
        if(n.eq.2) then
            xo = dcos(theta)*xi + dsin(theta)*zi
            yo = yi
            zo = -dsin(theta)*xi + dcos(theta)*zi
        endif
        if(n.eq.3) then
            xo = dcos(theta)*xi - dsin(theta)*yi
            yo = dsin(theta)*xi + dcos(theta)*yi
            zo = zi
        endif
            ro(2) = xo
            ro(3) = yo
            ro(4) = zo
end subroutine rotxyz
!--------------------------------------------

subroutine parameter_names(print_names,szpac,szpan)
    implicit none
    integer :: I
    logical :: print_names 
    double precision, dimension (100) :: szpac  
    character(len=8), dimension(100) :: szpan

    szpan = ''
    szpan(1) = "H0"
    szpan(2) = "Ol0"
    szpan(3) = "Om0"
    szpan(4) = "Omb"
    szpan(5) = "Omc"
    szpan(7) = "Om_rad"
    szpan(9) = "Om_cur"
    szpan(11) = "Ho/c"
    szpan(21) = "c/Ho"
    szpan(22) = "rho_0"
    szpan(24) = "c"
    szpan(30) = "z_init"
    szpan(51) = "contrast"
    szpan(52) = "radius"
    szpan(53) = "slope" 
    szpan(54) = "dipole"

    if ( print_names ) then
        do I=1,100
            print *, szpan(I),szpac(I)
        enddo
    endif

end subroutine parameter_names

!------------------------------------
subroutine parameter_values(npypac,pypac, npyszek,pyszek, szpac)
    implicit none
    integer :: npypac, npyszek
    double precision,  dimension (npypac) :: pypac
    double precision,  dimension (npyszek) :: pyszek        
    double precision,  dimension (100) :: szpac    


    double precision :: H0,Ho,little_h,z_initial
    double precision :: contrast, radius, slope, dipole    
    double precision :: omega_matter,omega_baryon,omega_cold
    double precision :: omega_lambda,omega_photon,omega_radiation,omega_curvature
    double precision :: w_matter,w_lambda,w_baryon,w_cold,w_photon,w_radiation
    double precision :: mass_unit,length_unit,time_unit,pi,G0,c0,T0,Neff,sbc,gcons,light_speed,kap,kapc2,gkr,gcr,lambda
    

! szpac:: contains parameters by convention:
! szpac(1-50)   -> parameters relate to homogeneous background
! szpac(51-90)  -> parameters relate to inhomogeneous setup 
! szpac(91-100) -> reserve for special flags 
    

  ! FIX needed: ensure this is just done one, not later
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
       contrast = -0.0015
       radius = 10.0
       slope  = 0.4
       dipole = 0.3
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
    z_initial = 1000.0  
    
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
    szpac(32) = pi 
    szpac(33) = pi/180.0d0
  

    szpac(51) = contrast
    szpac(52) = radius    
    szpac(53) = slope
    szpac(54) = dipole


    
end subroutine parameter_values



