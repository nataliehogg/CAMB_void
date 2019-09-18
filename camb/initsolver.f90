!#############################################################
!#This module computes the main equations for the interactive#
!#void model.                                                #
!#Current version updated on 08/06/2018                      #
!# NH modifications to improve speed 24/06/19                #
!#############################################################

module initsolver
use precision
use constants
use ModelParams

implicit none

!global variable for solver
real                               :: initial_z                   !starting scale factor for ODE
real                               :: final_z                     !final scale factor for ODE
integer                            :: nsteps    = 10000           !number of integration steps
real(dl), dimension(:),allocatable :: z_ode, solmat, solvoid      !8piG/3 * rho_m and rho_v
real(dl), dimension(:),allocatable :: ddsolmat, ddsolvoid         !same but derivatives obtained from spline

real                               :: coupling                    !value of the coupling parameter as taken from CAMB
integer                            :: model                       !choice of the interaction model we want to use
integer, parameter                 :: smooth_void=2 !possible options for q(z) binned reconstruction

logical                            :: debugging = .false.         !if T prints some files to check solver

contains

subroutine getcoupling(CP,z,rhov,Q)
!this subroutine just returns Q at any redshift
Type(CAMBparams) CP
real(dl), intent(in)  :: z
real, intent(in)      :: rhov
real(dl), intent(out) :: Q
integer               :: i

! old version here
!         if (z.gt.CP%zbins(CP%numvoidbins)) then
!            Q = -CP%qbins(CP%numvoidbins)*rhov
!         else
!            Q = CP%qbins(1)
!            do i=1,CP%numvoidbins-1
!               if (i.eq.1) then
!                  Q = Q + (CP%qbins(i+1)-CP%qbins(i))/2 * (1+tanh( CP%smoothfactor*(z-CP%zbins(i))/((CP%zbins(i))/2)  ) )
!               else
!                  Q = Q + (CP%qbins(i+1)-CP%qbins(i))/2 * (1+tanh( CP%smoothfactor*(z-CP%zbins(i))/((CP%zbins(i)-CP%zbins(i-1))/2)  ) )
!               end if
!            end do
!            Q = -Q*rhov
!         end if

if (z.gt.CP%zbins(CP%numvoidbins)) then
   Q = -CP%qbins(CP%numvoidbins)*rhov
else
   Q = CP%qbins(1)
   i=1
   Q = Q + (CP%qbins(i+1)-CP%qbins(i))/2 * (1+tanh( CP%smoothfactor*(z-CP%zbins(i))/((CP%zbins(i))/2)  ) )
   do i=2, CP%numvoidbins-1
     Q = Q + (CP%qbins(i+1)-CP%qbins(i))/2 * (1+tanh( CP%smoothfactor*(z-CP%zbins(i))/((CP%zbins(i)-CP%zbins(i-1))/2)  ) ) 
   end do
   Q = -Q*rhov
end if  

end subroutine getcoupling

subroutine getrhos(a,rho_m,rho_v) !NH removed declaration of error
!this subroutine gets the output of the differential equation
!interpolates it at the requested redshift
!outputs 8piG * rho_i

real(dl), intent(in)      :: a
real(dl), intent(out)     :: rho_m
real(dl), intent(out)     :: rho_v
real(dl)                  :: z
real(dl), parameter       :: azero = 0.0000001 !10**(-8.) !NH mod
! integer, optional :: error !Zero if OK

if (a.gt.azero) then
   z =  (1./a) - 1. !NHmod
else
   z = (1./azero) - 1. !NHmod
end if

if ((z.ge.initial_z).and.(z.le.final_z)) then
   call extrasplint(z_ode,solmat,ddsolmat,nsteps,z,rho_m)
   call extrasplint(z_ode,solvoid,ddsolvoid,nsteps,z,rho_v)
else
   rho_m = solmat(nsteps) * ( (1+z)/(1+z_ode(nsteps)) )**3.
   rho_v = solvoid(nsteps)
end if

! if ((rho_m.le.0._dl).or.(rho_v.le.0._dl)) then
!    write(*,*) 'negative densities in initsolver' !NH mod
!    return
! end if

end subroutine getrhos

subroutine deinterface(CP)
      Type(CAMBparams) CP
      integer, parameter      :: n = 1
      real, dimension(0:n)    :: x                              !dependent variables: rho_m, rho_v
      real                    :: h                              !step size
      real(dl)                :: rhoc_init, rhov_init
      integer                 :: i,j,k
      !debugging stuff
      real(dl)                :: debug_a, debug_c, debug_v, first_a_debug
      real(dl)                :: debug_q

      !initializing global ODE solver parameters from CAMB
      initial_z = 0._dl
      final_z   = CP%endred
      nsteps    = CP%numstepsODE

      !allocating arrays
      if (allocated(z_ode) .eqv. .false.) allocate(z_ode(nsteps+1), solmat(nsteps+1), solvoid(nsteps+1))
      if (allocated(ddsolmat) .eqv. .false.) allocate(ddsolmat(nsteps+1), ddsolvoid(nsteps+1))
      !if (allocated(GP_z) .eqv. .false.) allocate(GP_z(nsteps+1)) !NH mod
      !if (allocated(GP_q) .eqv. .false.) allocate(GP_q(nsteps+1))
      !if (allocated(dd_GP_q) .eqv. .false.) allocate(dd_GP_q(nsteps+1))

      !setting initial conditions for rho_c and rho_v at z=0
      rhoc_init = 3*(1000*CP%H0/c)**2.*CP%omegac               !8 pi G * rho_c^0
      rhov_init = 3*(1000*CP%H0/c)**2.*CP%omegav               !8 pi G * rho_V^0
      x = (/rhoc_init, rhov_init/)                             !initial conditions
      h = (log(1/(1+final_z)) - log(1/(1+initial_z)))/nsteps   !this is the step size for the ode

      !if (debugging) write(*,*) 'Started solving ODE'

      call rk4sys(CP,n,h,x)

      !if (debugging) write(*,*) 'solution done'

      !getting everything ready to interpolate
      call spline(z_ode,solmat,nsteps,1d30,1d30,ddsolmat)
      call spline(z_ode,solvoid,nsteps,1d30,1d30,ddsolvoid)

end subroutine deinterface

!differential equations utilities

!
! Numerical Mathematics and Computing, Fifth Edition
! Ward Cheney & David Kincaid
! Brooks/Cole Publ. Co.
! (c) 2003
!
! Section 11.1
!
! File: rk4sys.f90
!
! Runge-Kutta method of order 4 for a system of ode's (rk4sys,xpsys)


subroutine xpsys(CP,n,k,h,x,f)
      Type(CAMBparams) CP
      real, dimension (0:n) ::  x, f
      integer n,k
      real :: h      !stepsize
      real :: Hubble !Hubble parameter
      real(dl) :: Q      !interaction term
      real(dl) :: redshift, scalefac

      !Gets redshift for this step and the coupling
      scalefac = log(1/(1+initial_z)) + k*h
      redshift = -1+1/exp(scalefac)
      call getcoupling(CP,redshift,x(1),Q)

      !These are the actual derivatives
      !x'(1) = f(1)
      !x'(2) = f(2)

      f(0) = -(Q+3*x(0))!/(1+redshift)
      f(1) = Q!/(1+redshift)

end subroutine xpsys

subroutine rk4sys(CP,n,h,x)
      Type(CAMBparams) CP
      real ::  x(0:n)
      real, allocatable :: y(:), f(:,:)
      integer :: i, k, n
      real :: h


      z_ode(1) = initial_z
      solmat(1) = x(0)
      solvoid(1) = x(1)

      allocate (y(0:n), f(0:n,4))
out:  do k = 1,nsteps
        call xpsys(CP,n,k,h,x,f(0,1))
in1:    do i = 0,n
          y(i) = x(i) + 0.5*h*f(i,1)
        end do in1
        call xpsys(CP,n,k,h,y,f(0,2))
in2:    do i = 0,n
          y(i) = x(i) + 0.5*h*f(i,2)
        end do in2
        call xpsys(CP,n,k,h,y,f(0,3))
in3:    do i = 0,n
          y(i) = x(i) + h*f(i,3)
        end do in3
        call xpsys(CP,n,k,h,y,f(0,4))
in4:    do i = 0,n
          x(i) = x(i) + (h/6.0)* (f(i,1) + 2.0*(f(i,2) + f(i,3)) + f(i,4))
        end do in4
!        print *, k, x
        !storing functions at each step
        z_ode(k+1) = -1+1/exp(log(1/(1+initial_z)) + k*h)
        solmat(k+1)  = x(0)
        solvoid(k+1) = x(1)
      end do out
end subroutine rk4sys

    !-------------------------------------------------------------------
    SUBROUTINE extrasplint(xa,ya,y2a,n,x,y)
    INTEGER n
    real(dl)x,y,xa(n),y2a(n),ya(n)
    INTEGER k,khi,klo
    real(dl)a,b,h
    klo=1
    khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
            khi=k
        else
            klo=k
        endif
        goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+&
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
    END SUBROUTINE extrasplint
    !--------------------------------------------------------------------

end module initsolver

