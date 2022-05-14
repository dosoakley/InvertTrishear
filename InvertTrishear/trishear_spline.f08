!   This file is part of the program InvertTrishear
!    Copyright (C) 2015-2021  David Oakley
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along
!    with this program; if not, see <http://www.gnu.org/licenses/> or write 
!    to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
!    Boston, MA  02110-1301, USA..
!
!-----------------------------------------------------------------------------
!   Invert Trishear implements methods described in:
!    Oakley D.O.S. and Fisher, D.M, 2015, Inverse trishear modeling of bedding dip 
!    data using Markov chain Monte Carlo methods, Journal of Structural Geology, 
!    v. 80, p. 157-172.
!
!   David Oakley: david.o.oakley@uis.no
!   
!   If using this program in any academic or other publication, please acknowledge use of
!   the program. (This is not required by the GPL license, but is requested.)
!
!   If you make any changes to this program, please make note of that fact in this header.
!   (Prominent notice of changes is required by the terms of the GPL license).
!

module trishear_spline
!Spline fault.
!Clamped spline defined by (x,y) coordinates of control points and fault dip at the two ends.
!For simplicity, it is required that the fault tip that is fit for be the final fault tip.
!The splines are of the form yi = ai+bi*(x-xi)+ci*(x-xi)^2+di*(x-xi)^3
!This splits the spline into a number of straight segments and then calls functions from trishear_multibend.
!For now, the interpolation is into a fixed number of segments along a constant x interval, but this is not ideal.
!Something with a variable number of segments at a fixed change in angle would be best, but even a fixed number of segments along length 
!or a fixed number of segments evenly spaced by angle would be better.
!The fault tip is not a knot. Instead, the fault goes straight above the top knot and the fault tip is the end point of that straight segment.

logical :: nparam_vars_fixed = .false. !Tells whether we have yet fixed the nparam_start_terr, nparam_start_growth, and nparam_start_restored_fit variables. We should only have to do this once.
integer :: nknots !Number of knots (control points) for the spline. This includes the fault tip.

contains

subroutine trishear_spline_options
use options, only: options_id
use trishear_multibend, only: trishear_multibend_options,nangles
use parameters, only: nparams
implicit none
print*,'Enter number of spline knots (control points).'
print*,'This includes the bottom end point, but not the fault tip.'
do
    read(options_id,*) nknots
    if (nknots>=2) then
        exit
    else
        print*,'Error: You must have at least 2 knots'
    end if
end do
!Otherwise, this uses the same options as trishear_multibend.
print*,'Note: number of fault segments is the number of segments the spline will be divided into.'
print*,'This includes the horizontal detachment and the upper straight segment.'
call trishear_multibend_options
nparams = nparams-2*nangles+2+2*nknots !Correct the number of parameters from what trishear_multibend_options will have assigned to it.
!(This nparams removes the 2*nangles-1 multibend parameters and adds 2 ramp angles (top and bottom), x for all knots, and y for all but top knot.
end subroutine trishear_spline_options

subroutine run_model_spline(output,params)
use constants
use parameters, only: nparams,nparam_start_growth,nparam_start_terr,nparam_start_restored_fit
use trishear_multibend,only: nangles,run_model_multibend
use options, only: errs_id,ResultType
use data_uncertainties, only: nparam_start_sigma_bed,nparam_start_sigma_dip,nparam_start_sigma_terr,nparam_start_sigma_fault,&
    nparam_start_sigma_bed_restored,nparam_lc
use math, only: matinv
!use options, only: TipToSolve
implicit none
double precision,intent(out) :: output !the output result (RMS or probability)
double precision,dimension(nparams),intent(inout) :: params !values of the parameters
double precision :: tipx, tipy !Fault tip position (may be initial or final; doesn't really matter).
double precision :: ramp_angle_max,ramp_angle_max_deg !Maximum ramp angle, above which the fault becomes straight.
double precision :: ramp_angle_min,ramp_angle_min_deg !Minimum ramp angle, below which the fault becomes straight.
double precision,dimension(:),allocatable :: x_knots,y_knots !x and y coordinates of spline knots. These are listed from the fault tip down, same as ramp angles and y for bends in trishear_multi_bend
integer :: ramp_dir !1 = fault dips left, -1 = fault dips right.
double precision,dimension(:),allocatable :: a,b,c,d !Spline coefficients.
double precision :: alpha,beta !Spline slopes at the two ends (tan of ramp_angle_max and ramp_angle_min)
double precision,dimension(:),allocatable :: h !x_i+1-x_i for knots i = 1,...,n-1
double precision,dimension(:,:),allocatable :: A_mat !The 2D matrix for the matrix equation A*c = nu
double precision,dimension(:),allocatable :: nu !The 1D matrix nu for the matrix equation A*c = nu
double precision,dimension(nangles-1) :: bend_x,bend_y !x and y coordinates of bends in the fault.
double precision,dimension(nangles) :: ramp_angle !Fault ramp angles of all segments of the fault.
double precision :: x_step !x increment for interpolation.
!double precision :: angle_increment !Increment by which the ramp angle changes between segments.
integer :: i,n !Counters for do loops
double precision,dimension(:),allocatable :: params_for_multibend !The 5 elliptic fault parameters are replaced by nangles ramp angles and nangles-1 fault bends, for a net increase in the number of parameters of 2*nangles-6 (or 2*nangles-5 for a circle).
double precision :: output_from_multibend !The output value returned by trishear_multibend.
logical :: goodrun !tells if the run is good
double precision :: NaN = 0 !For creating NaN by dividing by zero
NaN = NaN/NaN
!Read in the relevant parameters.
tipx = params(1)
tipy = params(2)
ramp_angle_max_deg = params(4) !Maximum ramp angle
ramp_angle_max = ramp_angle_max_deg*pi/180.
ramp_angle_min_deg = params(5) !Minimum ramp angle. If you want a horizontal detachment, make this zero.
ramp_angle_min = ramp_angle_min_deg*pi/180.
!Determine whether the fault dips left or right
if (ramp_angle_max<=pi/2) then
    ramp_dir = 1 !Dips to the left
else
    ramp_dir = -1 !Dips to the right.
end if
!Get the knots.
allocate(x_knots(nknots),y_knots(nknots))
!x_knots(1) = tipx
!y_knots(1) = tipy
goodrun = .true.
do n = 1,nknots
    if (n==1) then !First knot has only an x coordinate b/c from tip, tipy, ramp_angle_max, and one coordinate, you can find the other one.
        x_knots(n) = params(6)
        y_knots(n) = tan(ramp_angle_max)*(x_knots(n)-tipx)+tipy
        if (x_knots(n)*ramp_dir > tipx*ramp_dir .or. y_knots(n)>tipy) then  !Each knot should be farther down the fault than the last one.
            goodrun = .false.
        end if
    else
        x_knots(n) = params(5+n)
        y_knots(n) = params(5+(nknots-1)+n) !Parameters file will list all the knot x coordinates and then all the knot y coordinates.
        if (x_knots(n)*ramp_dir > x_knots(n-1)*ramp_dir .or. y_knots(n)>y_knots(n-1)) then  !Each knot should be farther down the fault than the last one.
            goodrun = .false. !I could relax these requirements a bit by sorting the knots and only rejecting it if they can't be sorted in a way that is monotonic in both x and y.
        end if
    end if
end do
if (goodrun .eqv. .true.) then
    !Find the coefficients for the spline equations.
    allocate(a(nknots),b(nknots-1),c(nknots),d(nknots-1),h(nknots-1)) !There are only nknots-1 coefficients, but we need the extra entries in a and c for some calculations.
    a = y_knots
    h = x_knots(2:nknots)-x_knots(1:nknots-1)
    allocate(A_mat(nknots,nknots),nu(nknots))
    alpha = tan(ramp_angle_max)
    beta = tan(ramp_angle_min)
    A_mat = 0
    A_mat(1,1) = 2*h(1)
    A_mat(1,2) = h(1)
    nu(1) = 3*((1/h(1))*(y_knots(2)-y_knots(1))-alpha)
    A_mat(nknots,nknots) = 2*h(nknots-1)
    A_mat(nknots,nknots-1) = h(nknots-1)
    nu(nknots) = 3*(beta-(1/h(nknots-1))*(y_knots(nknots)-y_knots(nknots-1)))
    do n = 2,nknots-1
        A_mat(n,n-1) = h(n-1)
        A_mat(n,n) = 2*(h(n-1)+h(n))
        A_mat(n,n+1) = h(n)
        nu(n) = 3*((1/h(n))*(y_knots(n+1)-y_knots(n))-(1/h(n-1))*(y_knots(n)-y_knots(n-1)))
    end do
    c = matmul(matinv(nknots,A_mat),nu)
    b = (1/h)*(a(2:nknots)-a(1:nknots-1))-(h/3)*(2*c(1:nknots-1)+c(2:nknots))
    d = (1/(3*h))*(c(2:nknots)-c(1:nknots-1))
    !Resample the spline.
    !This is just a simple resample along x.
    x_step = (maxval(x_knots)-minval(x_knots))/(nangles-2)
    ramp_angle(1) = ramp_angle_max_deg
    ramp_angle(nangles) = ramp_angle_min_deg
    do n = 1,(nangles-1)
        bend_x(n) = x_knots(1)-ramp_dir*(n-1)*x_step
        do i = 1,nknots-1
            if (ramp_dir*bend_x(n)<=ramp_dir*x_knots(i) .and. ramp_dir*bend_x(n)>=ramp_dir*x_knots(i+1)) then
                bend_y(n) = a(i)+b(i)*(bend_x(n)-x_knots(i))+c(i)*(bend_x(n)-x_knots(i))**2+d(i)*(bend_x(n)-x_knots(i))**3
            end if
        end do
    end do
    !ramp_angle(1) = atan2(tipy-bend_y(1),tipx-bend_x(1))*180/pi
    do n = 2,(nangles-1) !Find the fault segment ramp angles.
        ramp_angle(n) = atan2(bend_y(n-1)-bend_y(n),bend_x(n-1)-bend_x(n))*180/pi
    end do
    !Prepare the parameters for sending to run_model_multibend
    allocate(params_for_multibend(nparams+2*nangles-2-2*nknots))
    params_for_multibend(1:3) = params(1:3)
    params_for_multibend(4:(3+nangles)) = ramp_angle
    params_for_multibend((4+nangles):(2+2*nangles)) = bend_y
    params_for_multibend((3+2*nangles):(nparams+2*nangles-2-2*nknots)) = params(5+2*nknots:nparams)
    if (nparam_vars_fixed .eqv. .false.) then !Fitting for restored dips or bed positions as parameters
        nparam_vars_fixed = .true.
        nparam_start_growth = nparam_start_growth+2*nangles-2-2*nknots
        nparam_start_terr = nparam_start_terr+2*nangles-2-2*nknots
        nparam_start_restored_fit = nparam_start_restored_fit+2*nangles-2-2*nknots
        nparam_start_sigma_bed = nparam_start_sigma_bed+2*nangles-2-2*nknots
        nparam_start_sigma_dip = nparam_start_sigma_dip+2*nangles-2-2*nknots
        nparam_start_sigma_terr = nparam_start_sigma_terr+2*nangles-2-2*nknots
        nparam_start_sigma_fault = nparam_start_sigma_fault+2*nangles-2-2*nknots
        nparam_start_sigma_bed_restored = nparam_start_sigma_bed_restored+2*nangles-2-2*nknots
        nparam_lc = nparam_lc+2*nangles-3-2*nknots
    end if
end if
!Check for some things not allowed, and if okay, proceed to run_model_multibend.
if (ramp_dir*ramp_angle_min>ramp_dir*ramp_angle_max) then
    goodrun = .false.
else if ((ramp_angle_min_deg == 0 .or. ramp_angle_min_deg == 180) .and. tipy<y_knots(nknots)) then
    goodrun = .false.
!else
!    goodrun = .true.
end if
if (goodrun .eqv. .true.) then
    call run_model_multibend(output_from_multibend,params_for_multibend)
    output = output_from_multibend
else
    if (ResultType == 1 .or. ResultType == 3) then !RMS or chisq
        output = NaN !This gives NaN.
    else !Probability
        output = -huge(0.) !Because these models have 0 probability of being correct.
    end if
end if
end subroutine run_model_spline

end module trishear_spline