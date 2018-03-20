!#############################################################
!#This module computes the main equations for the interactive#
!#void model.                                                #
!#Current version updated on 15/03/2018                      #
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
real(dl), dimension(:),allocatable :: GP_z, GP_q, dd_GP_q         !output arrays of GP reconstruction
real                               :: coupling                    !value of the coupling parameter as taken from CAMB
integer                            :: model                       !choice of the interaction model we want to use
integer, parameter                 :: theta_void=1, smooth_void=2 !possible options for q(z) binned reconstruction
integer, parameter                 :: GP_void=3, baseline_void=4  !possible options for q(z) gaussian process reconstruction

logical                            :: debugging = .false.         !if T prints some files to check solver

contains

subroutine getcoupling(CP,z,rhov,Q)
!this subroutine just returns Q at any redshift
Type(CAMBparams) CP
real(dl), intent(in)  :: z
real, intent(in)      :: rhov
real(dl), intent(out) :: Q
real                  :: multitheta !double theta function for binning
integer               :: i

      if (CP%void_model.eq.theta_void) then
         !Working with binned qV. No smoothing.

         if (z.gt.CP%zbins(CP%numvoidbins)) then
            Q = -CP%qbins(CP%numvoidbins)*rhov
         else
            Q = CP%qbins(1)
            do i=1,CP%numvoidbins-1
               multitheta = (sign(1d0,(z-CP%zbins(i)))+1)/2 - (sign(1d0,(z-CP%zbins(i+1)))+1)/2
               Q = Q + (CP%qbins(i+1)-CP%qbins(1))*multitheta
            end do
            Q = -Q*rhov
         end if

      else if (CP%void_model.eq.smooth_void) then
         !Working with binned qV smoothed with tanh.
         !Binned function used is based on Eq. 6 of 1703.01271


         if (z.gt.CP%zbins(CP%numvoidbins)) then
            Q = -CP%qbins(CP%numvoidbins)*rhov
         else
            Q = CP%qbins(1)
            do i=1,CP%numvoidbins-1
               if (i.eq.1) then
                  Q = Q + (CP%qbins(i+1)-CP%qbins(i))/2 * (1+tanh( CP%smoothfactor*(z-CP%zbins(i))/((CP%zbins(i))/2)  ) )
               else
                  Q = Q + (CP%qbins(i+1)-CP%qbins(i))/2 * (1+tanh( CP%smoothfactor*(z-CP%zbins(i))/((CP%zbins(i)-CP%zbins(i-1))/2)  ) )
               end if
            end do
            Q = -Q*rhov
         end if

      else if ((CP%void_model.eq.GP_void).or.(CP%void_model.eq.baseline_void)) then
         !interpolates what is obtained by GP in deinterface
         if (z.lt.CP%endred) then
            call extrasplint(GP_z,GP_q,dd_GP_q,nsteps,z,Q)
         else
            Q=0._dl
         end if
         Q = -Q*rhov
      else
         write(*,*) 'wait for it'
      end if

      !SP: debugging
      if (z>final_z) Q = 0._dl

end subroutine getcoupling

subroutine getrhos(a,rho_m,rho_v,error)
!this subroutine gets the output of the differential equation
!interpolates it at the requested redshift
!outputs 8piG * rho_i

real(dl), intent(in)      :: a
real(dl), intent(out)     :: rho_m
real(dl), intent(out)     :: rho_v
real(dl)                  :: z
real(dl), parameter       :: azero = 10**(-8.)
integer, optional :: error !Zero if OK

if (a.gt.azero) then
   z = -1.+1./a
else
   z = -1+1./azero
end if

if ((z.ge.initial_z).and.(z.le.final_z)) then
   call extrasplint(z_ode,solmat,ddsolmat,nsteps,z,rho_m)
   call extrasplint(z_ode,solvoid,ddsolvoid,nsteps,z,rho_v)
else
   rho_m = solmat(nsteps) * ( (1+z)/(1+z_ode(nsteps)) )**3.
   rho_v = solvoid(nsteps)
end if

if ((rho_m.le.0._dl).or.(rho_v.le.0._dl)) then
   global_error_flag         = 1
   global_error_message      = 'INITSOLVER: negative densities'
   if (present(error)) error = global_error_flag
   return
end if

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

      !Interface with GP python script
      character(LEN= 1000)                :: redbin
      character(LEN= 1000)                :: qbin
      character(LEN= 1000)                :: steps_de
      character(LEN= 1000)                :: z_ini
      character(LEN= 1000)                :: z_end
      character(LEN= 1000)                :: lencorr
      character(LEN= 20)                  :: feature_file="tmp_GPqz_000000.dat"
      character(LEN=10000)                :: command_plus_arguments
      real(dl), dimension(CP%numvoidbins) :: gpreds
      integer :: status
      integer :: getpid
      integer :: system


      !initializing global ODE solver parameters from CAMB
      initial_z = 0._dl
      final_z   = CP%endred
      nsteps    = CP%numstepsODE

      !allocating arrays
      if (allocated(z_ode) .eqv. .false.) allocate(z_ode(nsteps+1), solmat(nsteps+1), solvoid(nsteps+1))
      if (allocated(ddsolmat) .eqv. .false.) allocate(ddsolmat(nsteps+1), ddsolvoid(nsteps+1))
      if (allocated(GP_z) .eqv. .false.) allocate(GP_z(nsteps+1))
      if (allocated(GP_q) .eqv. .false.) allocate(GP_q(nsteps+1))
      if (allocated(dd_GP_q) .eqv. .false.) allocate(dd_GP_q(nsteps+1))

      if (debugging) then
         if ((CP%void_model.eq.theta_void).or.(CP%void_model.eq.smooth_void)) then
            write(*,*) 'num_bins=',CP%numvoidbins
            do k=1,CP%numvoidbins
               write(*,*) 'redshift',k,'=',CP%zbins(k)
               write(*,*) 'coupling',k,'=',CP%qbins(k)
            end do
         end if
      end if

      !setting initial conditions for rho_c and rho_v at z=0
      rhoc_init = 3*(1000*CP%H0/c)**2.*CP%omegac               !8 pi G * rho_c^0
      rhov_init = 3*(1000*CP%H0/c)**2.*CP%omegav               !8 pi G * rho_V^0
      x = (/rhoc_init, rhov_init/)                             !initial conditions
      h = (final_z - initial_z)/nsteps                         !step size for runge-kutta



      !Gaussian process interface
      if ((CP%void_model.eq.GP_void).or.(CP%void_model.eq.baseline_void)) then

         !Setting GP redshift to median redshift of each bin
         gpreds(1) = CP%zbins(1)/2
         do i=2,CP%numvoidbins
            gpreds(i) = (CP%zbins(i)+CP%zbins(i-1))/2.
         end do

         !Creating command line 
  
         !Generate tmp file name based on PID
         write (feature_file(11:16), "(Z6.6)"), getpid()
         !1. Prepare command and launch it!
         write(z_ini, "(E15.7)"      ) initial_z
         write(z_end, "(E15.7)"      ) final_z
         write(steps_de, "(I10)"     ) nsteps
         write(redbin, "(10E15.7)"   ) (gpreds(k),k=1,CP%numvoidbins)
         write(qbin, "(10f15.7)"     ) (CP%qbins(k),k=1,CP%numvoidbins) !python parser struggles with scientific notation negatives: using floats here
         write(lencorr, "(10E15.7)"  ) CP%corrlen

      
         if (CP%void_model.eq.GP_void) then
            if (debugging) write(*,*) 'WORKING WITH GP'
            !here needs the call to script with no baseline
            write(0,*) 'GP with no baseline not implemented yet'
            stop


         else if (CP%void_model.eq.baseline_void) then

            if (debugging) write(*,*) 'WORKING WITH GP (with baseline)'

            command_plus_arguments = "python camb/GP.py --inired "//trim(adjustl(z_ini))//" --endred "//trim(adjustl(z_end))//" --ODEsteps "//trim(adjustl(steps_de))// & 
            & " --redshifts "//trim(adjustl(redbin))// " --couplings "//trim(adjustl(qbin))// " --l "//trim(adjustl(lencorr))//" --outfile " // feature_file

            !calling script!!!
            if (debugging) then 
               write(*,*) 'Calling Gaussian process script with command line:'
               write(*,*) trim(adjustl(command_plus_arguments))
            end if
            status = system(trim(adjustl(command_plus_arguments)))
            if (status/=0) then
               print *, "Failure in GP reconstruction of q(z) -- see above."
               call abort
            end if

         end if

         !Reading temporary file generated by GP script--------------
         open(unit=17, file=feature_file, action='read')
         do i=1,nsteps
            read(17, "(E15.8, 1X, E15.8)", iostat=status) GP_z(i), GP_q(i)
            if (status>0) then
               print *, "Error reading the tmp q(z) file."
               call abort
            end if
         end do
         close(17, status='delete')
         !-----------------------------------------------------------

         !Setting interpolation for GP arrays------------------------
         call spline(GP_z,GP_q,nsteps,1d30,1d30,dd_GP_q)
         !-----------------------------------------------------------
      end if

      if (debugging) then
         write(*,*) '---------------------------------------------'
         write(*,*) 'STARTING DIFFERENTIAL EQUATION FOR COUPLED DE'
         write(*,*) 'initial settings:'
         write(*,*) 'model           =',CP%void_model
         write(*,*) 'initial redshift=',initial_z
         write(*,*) 'final redshift  =',final_z
         write(*,*) 'rho_c initial   =',rhoc_init
         write(*,*) 'rho_v initial   =',rhov_init
         write(*,*) 'points for ODE   =',CP%numstepsODE
         write(*,*) '---------------------------------------------'
      end if

      if (debugging) write(*,*) 'Started solving ODE'

      call rk4sys(CP,n,h,x)

      if (debugging) write(*,*) 'solution done'

      !getting everything ready to interpolate
      call spline(z_ode,solmat,nsteps,1d30,1d30,ddsolmat)
      call spline(z_ode,solvoid,nsteps,1d30,1d30,ddsolvoid)


      if (debugging) then
         first_a_debug = 1.e-4
         if (CP%void_model.eq.theta_void) then
            open(42, file='solutions_thetabin.dat')
            open(666,file='binned_coupling_thetabin.dat')
         else if (CP%void_model.eq.smooth_void) then
            open(42, file='solutions_smoothbin.dat')
            open(666,file='binned_coupling_smoothbin.dat')
         else if (CP%void_model.eq.GP_void) then
            open(42, file='solutions_GP.dat')
            open(666,file='binned_coupling_GP.dat')
         else if (CP%void_model.eq.baseline_void) then
            open(42, file='solutions_GPbaseline.dat')
            open(666,file='binned_coupling_GPbaseline.dat')
         end if
         do k=1,nsteps
            debug_a = first_a_debug+k*(1.-first_a_debug)/nsteps
            call getrhos(debug_a,debug_c, debug_v)
            write(42,*) debug_a, debug_c/(debug_c+debug_v), debug_v/(debug_c+debug_v)
            !write(42,*) debug_a, debug_c, debug_v
            call getcoupling(CP,-1+1/debug_a,real(debug_v),debug_q)
            write(666,*) -1+1/debug_a, -debug_q/debug_v
         end do
         close(42)
         close(666)
      end if

      if (debugging) then
        write(*,'("Om_b h^2             = ",f9.6)') CP%omegab*(CP%H0/100)**2
        write(*,'("Om_c h^2             = ",f9.6)') CP%omegac*(CP%H0/100)**2
        write(*,'("Om_nu h^2            = ",f9.6)') CP%omegan*(CP%H0/100)**2
        write(*,'("Om_Lambda            = ",f9.6)') CP%omegav
        write(*,'("Om_K                 = ",f9.6)') CP%omegak
        write(*,'("Om_m (1-Om_K-Om_L)   = ",f9.6)') 1-CP%omegak-CP%omegav
        write(*,'("100 theta (CosmoMC)  = ",f9.6)') 100*CosmomcTheta()
      end if

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
      real(dl) :: redshift

      !Gets redshift for this step and the coupling
      redshift = initial_z + k*h
      call getcoupling(CP,redshift,x(1),Q)

      !These are the actual derivatives
      !x'(1) = f(1)
      !x'(2) = f(2)

      f(0) = (Q+3*x(0))/(1+redshift)
      f(1) = -Q/(1+redshift)

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
        z_ode(k+1) = initial_z + k*h
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


!subroutine getfunc(zbin,fbin,z0,f0,fder0,z,outfunc,outder)

!real(dl), dimension(nbin), intent(in) :: fbin
!real(dl), dimension(nbin), intent(in) :: zbin
!real(dl), dimension(nbin-1) :: fderbin
!real(dl), dimension(nbin-1) :: zderbin
!real(dl), intent(in) :: z
!real(dl), intent(in) :: z0
!real(dl), intent(in) :: f0
!real(dl), intent(in) :: fder0
!real(dl), intent(out) :: outfunc, outder

!integer i,j,k


!if (z.lt.zbin(1)) then

!   outfunc = fbin(1)
!   outder = 0.d0!fderbin(1)

!else if (z.gt.zbin(nbin)) then

!   outfunc = fbin(nbin)
!   outder = 0.d0!fderbin(nbin-1)

!else

!   outfunc = (fbin(1)+f0)/2 + (fbin(1)-f0)/2 * tanh(((z-zbin(1))/(zbin(2)-zbin(1)))*h)

!   do i=1,nbin-1
!      outfunc = outfunc + (fbin(i+1)-fbin(i))/2 * (1+tanh(((z-zbin(i+1))/(zbin(i+1)-zbin(i)))*h))
!   end do


!end if

!end subroutine getfun

end module initsolver
