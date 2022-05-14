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

module trishear_elliptic_approx
!This model approximates an elliptic listric fault, using a finite number of straight segments to approximate the ellipse.
!The approximate fault, composed of straight segments, is then passed to trishear_multi_bend.
!Currently, the ellipse is required to have its axes parallel to the horizontal and vertical axes of the world Cartesian coordinate system,
!and it is required to shallow into a horizontal detachment.
!The parameters to be fit for are: detachment depth, horizontal semiaxis (termed a, whether truly the major semiaxis or not), 
!vertical semiaxis (termed b, whether truly the minor semiaxis or not), and the maximum ramp angle.
!Note that the maximum angle is what it would be if this were a true ellipse. Since we approximate it with straight segments, the ramp angle
!of the steepest segment within the ellipse will in fact be slightly less, but the ramp angle continuing up from this point will be at the 
!maximum angle.
!Above the maximum angle, the fault becomes straight, as a single segment.
!The segments are created by dividing the ellipse into segments with an equal change in ramp angle (dip) from one segment to the next.
!It seems that the change in angle from the detachment to the elliptical part or the elliptical part to the straight part are half what they are for the other steps.
!This seems somewhat inevitable, since at that point, the ramp dip is equal to the dip of the tangent to the ellipse, while elsewhere the dip between two points on the ellipse is equal to the average of the tangents at those two points.
!I think the only way to have all dips be equal would be to approximate the ellipse with a series of tangents to it rather than with a series of straight lines between points on it.

logical :: nparam_vars_fixed = .false. !Tells whether we have yet fixed the nparam_start_terr, nparam_start_growth, and nparam_start_restored_fit variables. We should only have to do this once.
integer :: fault_shape !Tells the shape of the fault to be approximated by segments. 1 = circle, 2 = ellipse.
integer :: circle !This is 1 if it's a circle and 0 otherwise. For subtracting 1 from the number of paremters for a circle rather than an ellipse.
integer :: ellipse_parameterization !How to parameterize the ellipse. (1) y_base, a, b; (2) y_base, y_str, e
integer :: decoly_at_tipy !Tells whether decollement depth and tipy must be the same (1) or not (0)

contains

subroutine trishear_elliptic_options
use options, only: options_id,TipToSolve
use trishear_multibend, only: trishear_multibend_options,nangles,bend_at_tipy
use parameters, only: nparams
implicit none
character (len = 25) :: str
print*,'Shape: '
print*,'(1) Circle'
print*,'(2) Ellipse'
if (options_id == 5) then !Manual input
    read(options_id,*) fault_shape
else
    read(options_id,*) str
    if (str == 'circle') then                    
        fault_shape = 1
        circle = 1
    else if (str == 'ellipse') then
        fault_shape = 2
        circle = 0
    else
        backspace(options_id)
        fault_shape = 2 !Default from before I added the circle option, since this question wasn't asked then.
        circle = 0
    end if    
end if
if (fault_shape==2) then
    do
        print*,'How should the ellipse geometry be parameterized?'
        print*,'(Note: Base and top elevation are where the fault becomes straight.)'
        print*,'    (1) base elevation, horizontal semiaxis, vertical semiaxis'
        print*,'    (2) base elevation, top elevation, eccentricity'
        if (options_id == 5) then !Manual input
            read(options_id,*) ellipse_parameterization
        else
            read(options_id,*) str
            if (str=='ellipse_parameterization') then
                read(options_id,*) ellipse_parameterization
            else
                backspace(options_id)
                ellipse_parameterization = 1 !Default, initially the only option.
            end if
        end if
        if (ellipse_parameterization==1) then
            exit
        else if (ellipse_parameterization==2) then
            print*,'For the purpose of this program, ellipticity should be positive if the horizontal semiaxis is longer'
            print*,'and negative if the vertical semiaxis is longer.'
            exit
        else
            print*,'Error: Unrecognized ellipse parameterization.'
        end if
    end do
else
    ellipse_parameterization=1 !For a circle, parameterization in terms of the radius of curvature is the only one allowed.
end if
if (TipToSolve==1) then
    !See if the fault should begin propagating from the detachment or above it.
    !This is only an option if fitting for the initial fault tip position, since in that case tipy will also be the detachment depth so it's easier that way.
    do
        print*,'Must detachment depth be at initial fault tip depth? (0/1)'
        if (options_id == 5) then !Manual input
            read(options_id,*) decoly_at_tipy
        else
            read(options_id,*) str
            if (str=='decoly_at_tipy') then
                decoly_at_tipy = 1
            else
                backspace(options_id)
                decoly_at_tipy = 0 !Default, initially the only option.
            end if
        end if
        if (decoly_at_tipy==0 .or. decoly_at_tipy==1) then
            exit
        else
            print*,'Error: Unrecognized option.'
        end if
    end do
else
    decoly_at_tipy = 0
end if
!Otherwise, this uses the same options as trishear_multibend.
print*,'Note: number of fault segments is the number of segments the listric fault will be divided into.'
print*,'This includes the horizontal detachment and the upper straight segment.'
call trishear_multibend_options
nparams = nparams-2*nangles+6-circle-decoly_at_tipy !Correct the number of parameters from what trishear_multibend_options will have assigned to it.
end subroutine trishear_elliptic_options

subroutine run_model_elliptic(output,params)
use constants
use parameters, only: nparams,nparam_start_growth,nparam_start_terr,nparam_start_restored_fit
use trishear_multibend,only: nangles,run_model_multibend,bend_at_tipy,bend_start
use options, only: errs_id,ResultType
use data_uncertainties, only: nparam_start_sigma_bed,nparam_start_sigma2_bed,nparam_start_sigma_dip,nparam_start_sigma_terr,&
    nparam_start_sigma_fault,nparam_start_sigma_bed_restored,nparam_lc
!use options, only: TipToSolve
implicit none
double precision,intent(out) :: output !the output result (RMS or probability)
double precision,dimension(nparams),intent(inout) :: params !values of the parameters
double precision :: tipx, tipy !Fault tip position (may be initial or final; doesn't really matter).
double precision :: ramp_angle_max,ramp_angle_max_deg !Maximum ramp angle, above which the fault becomes straight.
double precision :: ramp_angle_min,ramp_angle_min_deg !Minimum ramp angle, below which the fault becomes straight.
double precision :: y_base,x_base !x and y coordinates of the horizontal detachment or lowest point of the listric fault at which it becomes straight.
integer :: ramp_dir !1 = fault dips left, -1 = fault dips right.
double precision :: a,b !Horizontal and vertical axes of the ellipse.
double precision :: e,ba_ratio !eccentricity, ratio b/a
double precision :: xc,yc !Coordinates of the center of the ellipse.
double precision :: x_str,y_str !Coordinates of the point at which the fault changes from elliptical to straight.
double precision :: theta,theta_max,theta_min !Polar coordinate of a fault bend point.
double precision,dimension(nangles-1) :: bend_x,bend_y !x and y coordinates of bends in the fault.
double precision,dimension(nangles) :: ramp_angle !Fault ramp angles of all segments of the fault.
double precision :: angle_increment !Increment by which the ramp angle changes between segments.
integer :: n !A counter of ramp bends or segments.
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
if (decoly_at_tipy == 0) then
    y_base = params(6) !Detachment depth
else
    y_base = tipy
end if
if (ellipse_parameterization==1) then
    a = params(7-decoly_at_tipy) !Horizontal semiaxis
    if (fault_shape==2) then !ellipse
        b = params(8-decoly_at_tipy) !Vertical semiaxis
    else !circle
        b = a;
    end if
    ba_ratio = b/a
else
    y_str = params(7-decoly_at_tipy)
    e = params(8-decoly_at_tipy)
    if (e>=0) then
        ba_ratio = sqrt(1-e**2)
    else
        ba_ratio = 1/sqrt(1-e**2)
    end if
end if
!Determine whether the fault dips left or right
if (ramp_angle_max<=pi/2) then
    ramp_dir = 1 !Dips to the left
else
    ramp_dir = -1 !Dips to the right.
end if
!Convert maximum and minimum ramp angles to theta values.
theta_max = ramp_dir*atan(-ba_ratio/tan(ramp_angle_max))
if (ramp_angle_min_deg == 0 .or. ramp_angle_min_deg == 180) then !In this case, tan(ramp_angle_min) = 0, so we would end up with 0 in the denominator of the theta_min equation.
    theta_min = -pi/2
else
    theta_min = ramp_dir*atan(-ba_ratio/tan(ramp_angle_min))
end if
!Calculate a and b if still needed.
if (ellipse_parameterization==2) then
    b = (y_str-y_base)/(sin(theta_max)-sin(theta_min))
    a = b/ba_ratio
end if
!Calculate the center coordinates (xc, yc).
yc = y_base-b*sin(theta_min)
!Calculate y_str if still needed.
if (ellipse_parameterization==1) then
    y_str = yc+b*sin(theta_max) !Elevation at which the elliptical part gives way to a straight fault.
end if
if (decoly_at_tipy == 1) then !This should be covered by the "else" case, but doing it this way avoids potential rounding errors.
    xc = tipx
    x_str = xc+ramp_dir*a*cos(theta_max)
    x_base = tipx
else if (tipy>y_str) then !Tip is above the elliptical part of the fault.
    x_str = tipx+(y_str-tipy)/tan(ramp_angle_max)
    xc = x_str-ramp_dir*sqrt((a**2)*(1.-((y_str-yc)**2)/b**2.))
    x_base = xc+ramp_dir*a*cos(theta_min)
else if (tipy<y_base) then !Tip is below the elliptical part of the fault.
    x_base = tipx+(y_base-tipy)/tan(ramp_angle_min)
    xc = x_base-ramp_dir*sqrt((a**2)*(1.-((y_base-yc)**2)/b**2.))
    x_str =  xc+ramp_dir*a*cos(theta_max)
else !Tip is within the elliptical part of the fault.
    xc = tipx-ramp_dir*sqrt((a**2)*(1.-((tipy-yc)**2)/b**2))
    x_str = xc+ramp_dir*a*cos(theta_max)
    x_base = xc+ramp_dir*a*cos(theta_min)
end if
ramp_angle(1) = ramp_angle_max_deg
bend_y(1) = y_str
bend_x(1) = x_str
ramp_angle(nangles) = ramp_angle_min_deg
bend_y(nangles-1) = y_base
bend_x(nangles-1) = x_base
angle_increment = abs(ramp_angle_max_deg-ramp_angle_min_deg)/(nangles-2)
do n = 2,(nangles-2) !Find the bend x and y coordinates
    theta = ramp_dir*atan(-b/(a*tan((ramp_angle_max_deg-ramp_dir*(n-1)*angle_increment)*pi/180.)))
    bend_x(n) = xc+ramp_dir*a*cos(theta)
    bend_y(n) = yc+b*sin(theta)
end do
do n = 2,(nangles-1) !Find the fault segment ramp angles.
    ramp_angle(n) = atan2(bend_y(n-1)-bend_y(n),bend_x(n-1)-bend_x(n))*180/pi
end do
if (decoly_at_tipy == 1) then
    bend_at_tipy = 1
    bend_start = nangles-1
end if
allocate(params_for_multibend(nparams+2*nangles-6+circle))
params_for_multibend(1:3) = params(1:3)
params_for_multibend(4:(3+nangles)) = ramp_angle
if (decoly_at_tipy == 1) then
    params_for_multibend((4+nangles):(2+2*nangles-1)) = bend_y(1:nangles-2)
else
    params_for_multibend((4+nangles):(2+2*nangles)) = bend_y
end if
params_for_multibend((3+2*nangles-decoly_at_tipy):(nparams+2*nangles-6+circle)) = params((9-circle-decoly_at_tipy):nparams)
if (nparam_vars_fixed .eqv. .false.) then !Fitting for restored dips or bed positions as parameters
    nparam_vars_fixed = .true.
    nparam_start_growth = nparam_start_growth+2*nangles-6+circle
    nparam_start_terr = nparam_start_terr+2*nangles-6+circle
    nparam_start_restored_fit = nparam_start_restored_fit+2*nangles-6+circle
    nparam_start_sigma_bed = nparam_start_sigma_bed+2*nangles-6+circle
    nparam_start_sigma2_bed = nparam_start_sigma2_bed+2*nangles-6+circle
    nparam_start_sigma_dip = nparam_start_sigma_dip+2*nangles-6+circle
    nparam_start_sigma_terr = nparam_start_sigma_terr+2*nangles-6+circle
    nparam_start_sigma_fault = nparam_start_sigma_fault+2*nangles-6+circle
    nparam_start_sigma_bed_restored = nparam_start_sigma_bed_restored+2*nangles-6+circle
    nparam_lc = nparam_lc+2*nangles-6+circle
end if
!Check for some things not allowed, and if okay, proceed to run_model_multibend.
if (ramp_dir*ramp_angle_min>ramp_dir*ramp_angle_max) then
    goodrun = .false.
else if ((ramp_angle_min_deg == 0 .or. ramp_angle_min_deg == 180) .and. tipy<y_base) then
    goodrun = .false.
else if (y_str < y_base) then
    goodrun = .false.
else
    goodrun = .true.
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
end subroutine run_model_elliptic

end module trishear_elliptic_approx