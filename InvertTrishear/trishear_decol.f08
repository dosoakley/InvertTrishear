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

module trishear_decol

integer :: decoly_at_tipy !Tells whether decollement depth and tipy must be the same (1) or not (0)

contains

subroutine trishear_decol_options
use parameters, only: nparams
use options, only: options_id
implicit none
do
    print*,'Must detachment depth be at initial fault tip depth? (0/1)'
    read(options_id,*) decoly_at_tipy
    if (decoly_at_tipy == 0) then
        nparams = 8
        exit
    else if (decoly_at_tipy == 1) then
        nparams = 7
        exit
    else
        print*,'Invalid Command'
    end if
end do
end subroutine trishear_decol_options

subroutine run_model_decol(output,params) !Subroutine to run a model and return the RMS error
use constants
use parameters, only: nparams,increment,v0,slip_sense,nparam_start_restored_fit, & 
    nrestored_fit_params,nparam_start_growth,nparam_start_terr
use options, only: ResultType,DatTypes,TipToSolve,ndatatypes_poss,DatType_beds,DatType_dips,DatType_faultpts,DatType_terr,FitType, &
    terr_age_order,DatType_beds_restored,beds_age_order,beds_age_order,limit_bed_dips,min_bed_slopes,max_bed_slopes,FitType
use math, only: xy_to_ze2D,ze_to_xy2D
use data_module
use err_and_prob
use data_uncertainties, only: GetSigmaBed,GetSigmaDip,GetSigmaTerr,GetSigmaFault,GetSigmaBedRestored,GetLc,GetSigma2Bed
use growth_strata, only: GetGrowthSlip
use beds_order, only: check_beds_order
implicit none
double precision,intent(out) :: output !the output result (RMS or probability)
double precision,dimension(ndatatypes_poss) :: outputs !A single vector holding outputs for beds, dips, and terraces (in that order).
double precision,dimension(nparams),intent(inout) :: params !values of the parameters
double precision :: tipx,tipy,total_slip,total_slip_mag,ramp_angle,ramp_angle_deg,ramp_acute,phi,phi_deg,PoverS,s,m !parameters
double precision :: slip !Slip to restore a bed, may be growth or pregrowth.
double precision,dimension(nterr) :: terr_slip !Slip needed to restore each terrace.
double precision,dimension(ngrowth) :: growth_slip !Slip needed to restore each growth strata bed.
double precision,dimension(nrestored_fit_params) :: restored_fit_params
double precision,dimension(2,1) :: decolxy,decolze !x,y and zeta,eta coordinates of where the decollement meets the ramp
double precision :: tipxinit, tipyinit, tipxfinal, tipyfinal !Initial and final (in terms of forward motion) tip positions
integer :: ramp_dir !dip direction of ramp
double precision :: tan_ramp,sin_ramp_acute,cos_ramp_acute !Trig functions of the ramp angle.
double precision :: NaN = 0 !For creating NaN by dividing by zero
integer :: bedlength !number of points in the bed being worked on
double precision,dimension(:,:),allocatable :: thisbed,thisbedze !the bed being worked on
double precision :: slope,b !Slope and intercept of a bed
double precision :: slope_last, b_last !Slope and intercept of the previous bed restored.
double precision,dimension(2) :: bed_xlims,bed_xlims_last !x limits of current and previous bed restored.
double precision,dimension(:,:),allocatable :: dipposze !dip positions in zeta, eta coordinates
double precision,dimension(ndips) :: dipsze !dips in zeta, eta coordinates
double precision,dimension(ndips) :: RestoredDips !values of dips after inverse fault movement is complete
double precision,dimension(nbeds) :: beds_output_indiv !Output from each of the beds individually.
double precision,dimension(n_restored_beds) :: restored_beds_output_indiv !Output from each of the beds individually.
double precision,dimension(nterr) :: terr_output_indiv !Output from each of the terraces individually.
double precision,dimension(nfaultpts) :: fault_errs !Errors from each fault point.
double precision,dimension(ndips) :: sigma_dip_array !Array of sigmas for dips
double precision,dimension(2,2) :: R, Refl !The rotation matrix for a rotation by ramp_acute and reflection matrix for faults dipping to the right.
double precision,dimension(:,:,:),allocatable :: Cbed,Cterr !Covariance matrix for uncertainties in bed point positions and terrace point positions.
double precision,dimension(:,:,:),allocatable :: Cbed2 !Covariance matrix for second (correlated) uncerrtainties in bed point positions if ResultType is 5.
double precision,dimension(nbeds) :: sigma_bed_x,sigma_bed_y !Uncertainty in bed points x and y coordinates before any propagation.
double precision,dimension(nbeds) :: sigma2_bed_x,sigma2_bed_y !Uncertainty in bed points x and y coordinates for second (correlated) uncertainty if ResultType is 5.
double precision :: Lc !Correlation length for bed data if correlated.
double precision :: sigma_bed_restored_x,sigma_bed_restored_y !Uncertainty in x and y coordinates of points on the restored beds, when these are used as independent data.
double precision :: sigma_dip,sigma_restored_dip !Uncertainty in dips in the deformed and restored states.
double precision,dimension(nterr) :: sigma_terr,sigma_orig_elev !Uncertainty in terrace point positions and in the original elevation of the terrace inner edge.
double precision :: sigma_fault !Uncertainty in fault points.
integer :: i,n,g ! Counters. n and g are for all beds and growth beds respectively.
logical :: goodrun !Tells if the run is good.
NaN = NaN/NaN
goodrun = .true. !Initialize it to true. So far, this is a good run.
do !This do loop lets us exit anytime it turns out we have an impossible combination of parameters
!Make assignments
tipx = params(1)
tipy = params(2)
total_slip_mag = params(3)
ramp_angle_deg = params(4)
phi_deg = params(5)
PoverS = params(6)
s = params(7)
if (DatTypes(DatType_beds) == 1 .and. ngrowth > 0) then !For growth strata, need additional parameters
    call GetGrowthSlip(params(nparam_start_growth:nparam_start_growth+ngrowth-1),total_slip_mag,slip_sense,growth_slip,goodrun)
    if (goodrun .eqv. .false.) then
        exit
    end if
end if
if (DatTypes(DatType_terr) == 1) then !For terraces, need additional parameters
    terr_slip = params(nparam_start_terr:nparam_start_terr+nterr-1)
    if (maxval(terr_slip) > total_slip_mag) then !Note: terr_slip at this point should be positive.
        goodrun = .false.
        exit
    end if
    if (terr_age_order==1 .and. nterr>1) then !Youngest to oldest
        if (minval(terr_slip(2:nterr)-terr_slip(1:(nterr-1)))<0) then
            goodrun = .false.
        end if
    else if(terr_age_order==2 .and. nterr>1) then !Oldest to youngest
        if (minval(terr_slip(1:(nterr-1))-terr_slip(2:nterr))<0) then
            goodrun = .false.
        end if
    end if
    terr_slip = terr_slip*slip_sense
end if
if (FitType==5 .or. FitType==6) then !Fitting for restored dips or bed positions as parameters
    restored_fit_params = params(nparam_start_restored_fit:nparam_start_restored_fit+nrestored_fit_params-1)
end if
!Take care of some degree-radian conversions, etc.
total_slip = total_slip_mag*slip_sense
ramp_angle = ramp_angle_deg*pi/180; !convert to radians
if (ramp_angle<=pi/2) then !1=dips left, -1 = dips right
    ramp_dir = 1
else
    ramp_dir = -1
end if
ramp_acute = min(ramp_angle,pi-ramp_angle) !Acute version of ramp angle
tan_ramp = tan(ramp_angle)
phi = phi_deg*pi/180 !convert to radians
m = tan(phi)
if (phi == pi/2) then
        m = tan(89*pi/180) !Since tan(pi) = inf, which causes problems
end if
if (phi >= max(ramp_angle,pi-ramp_angle)/2) then !Trishear boundary impinges on fold axis, which we don't want.
    goodrun = .false.
    exit
end if
!Calculate initial and final tipx and tipy (final = present-day, initial = when trishear motion began)
if (TipToSolve == 1) then !If tipx, tipy are initial
    tipxinit = tipx
    tipyinit = tipy
    tipxfinal = tipx+slip_sense*total_slip*PoverS*cos(ramp_angle) !Remember that if doing an inversion, total_slip and slip_sense will be negative
    tipyfinal = tipy+slip_sense*total_slip*PoverS*sin(ramp_angle)
else !If tipx, tipy are final
    tipxinit = tipx-slip_sense*total_slip*PoverS*cos(ramp_angle)
    tipyinit = tipy-slip_sense*total_slip*PoverS*sin(ramp_angle)
    tipxfinal = tipx
    tipyfinal = tipy
end if
!Get the appropriate decollement coordinates
if (decoly_at_tipy == 0) then
    decolxy(2,1) = params(8) !decol y
    decolxy(1,1) = tipxfinal-(tipyfinal-decolxy(2,1))/tan_ramp; !x coordinate of where decollement meets ramp
else
    decolxy(2,1) = tipyinit !decol y
    decolxy(1,1) = tipxinit !decol x
end if
if (tipyinit <= decolxy(2,1)) then !tipyinit isn't allowed to be below decol y.
    goodrun = .false.
    exit
end if
decolze = xy_to_ze2D(decolxy,tipxfinal,tipyfinal,ramp_acute,ramp_dir)
!Initialize the outputs array.
outputs = 0 !For RMS or chi squared it should be 0 anyway. For probabilities, p should be 1, so ln(p) = 0.
!Now restore data.
R = reshape((/cos(ramp_acute),sin(ramp_acute),-sin(ramp_acute),cos(ramp_acute)/),(/2,2/)) !Rotation matrix (counterclockwise by ramp angle). This should reshape to (in matlab syntax) [cos, -sin; sin, cos].
Refl = reshape((/ramp_dir,0,0,1/),(/2,2/))
if (DatTypes(DatType_beds) == 1) then !Beds
    g = 1 !A counter of growth strata
    call GetSigmaBed(params,sigma_bed_x,sigma_bed_y) !Get the uncertainties in bed points.
    call GetLc(params,Lc) !Get the correlation length.
    if (ResultType==5) call GetSigma2Bed(params,sigma2_bed_x,sigma2_bed_y) !Get the uncertainties in bed points.
    if(DatTypes(DatType_beds_restored) == 1) then
        call GetSigmaBedRestored(params,sigma_bed_restored_x,sigma_bed_restored_y) !Get the uncertainty in restored bed points.
    end if
    do n = 1,nbeds !Go through each of the beds
        bedlength = beds(n)%npts
        allocate(thisbed(2,bedlength))
        allocate(thisbedze(2,bedlength))
                allocate(Cbed(2,2,bedlength))
        Cbed = 0 !Initialize covariance matrix for bed points (x,y).
        Cbed(1,1,:) = sigma_bed_x(n)**2
        Cbed(2,2,:) = sigma_bed_y(n)**2
        do i = 1,bedlength !Reflect, then rotate.
            Cbed(:,:,i) = matmul(matmul(Refl,Cbed(:,:,i)),Refl) !No need to transpose Refl, because transpose(Refl) = Refl.
            Cbed(:,:,i) = matmul(matmul(transpose(R),Cbed(:,:,i)),R) !Change Cbed into the rotated (zeta,eta) coordinate system. I think this is only necessary if sigma_x /= sigma_y to start or if there is covariance to start out.
        end do
        thisbed(:,:) = beds(n)%pts
        thisbedze = xy_to_ze2D(thisbed,tipxfinal,tipyfinal,ramp_acute,ramp_dir); !Convert to eta, zeta
        if (beds(n)%growth .eqv. .true.) then !For growth strata
            slip = growth_slip(g)
            g = g+1
        else
            slip = total_slip
        end if
        call trishear_func_decol(thisbedze,v0,m,slip,increment,-slip_sense*PoverS,s,decolze,ramp_acute,Cbed) !Run the trishear on this bed
        thisbed = ze_to_xy2D(thisbedze,tipxfinal,tipyfinal,ramp_acute,ramp_dir) !Convert back to regular coordinate system.
        do i = 1,bedlength !Rotate, then Reflect.
            Cbed(:,:,i) = matmul(matmul(R,Cbed(:,:,i)),transpose(R)) !Rotate the covariance matrix from (zeta,eta) back into (x,y) coordinates.
            Cbed(:,:,i) = matmul(matmul(Refl,Cbed(:,:,i)),Refl) !No need to transpose Refl, because transpose(Refl) = Refl.
        end do
        if (ResultType==5) then !Correlated and uncorrelated probabilities.
            allocate(Cbed2(2,2,bedlength))
            Cbed2 = 0 !Initialize covariance matrix for bed points (x,y).
            Cbed2(1,1,:) = sigma2_bed_x(n)**2
            Cbed2(2,2,:) = sigma2_bed_y(n)**2
        end if
        call calc_bed_err(n,bedlength,thisbed,beds_output_indiv,slope,b,Cbed,Lc,restored_fit_params,Cbed2)
        if (DatTypes(DatType_beds_restored) == 1) then !Points on the restored beds.
            do i = 1,n_restored_beds
                if (beds(n)%ident == restored_beds(i)%ident) then
                    call calc_restored_bed_err(i,slope,b,restored_beds_output_indiv(i),sigma_bed_restored_x,sigma_bed_restored_y,&
                        restored_fit_params)
                end if
            end do
        end if
        if (beds_age_order==1 .or. beds_age_order == 2) then
            bed_xlims = (/minval(thisbed(1,:)),maxval(thisbed(1,:))/)
            if (check_beds_order(n,b,slope,b_last,slope_last,bed_xlims,bed_xlims_last) .eqv. .false.) then
                goodrun = .false.
            end if
            b_last = b !For the next bed, this will be be from the last bed for comparison.
            slope_last = slope
            bed_xlims_last = bed_xlims
        end if
        if (FitType==3 .and. limit_bed_dips==1) then
            if (slope<min_bed_slopes(n) .or. slope>max_bed_slopes(n)) then
                goodrun = .false.
            end if
        end if
        deallocate(thisbed,thisbedze,Cbed)
        if (ResultType==5) deallocate(Cbed2)
    end do
    call calc_beds_err(beds_output_indiv,outputs(1))
    if (DatTypes(DatType_beds_restored) == 1) then
        call calc_beds_err(restored_beds_output_indiv,outputs(5)) !Note: This won't work properly if the output is RMS error, but currently you're not allowed to use RMS error with more than one data type.
    end if
end if
if (goodrun .eqv. .false.) then !If the above numerical solutions didn't work out, then discard the run.
    exit !Leave the loop now before wasting any more time or risking issues with later parts of the run.
end if
if (DatTypes(DatType_dips) == 1) then !Dips
    call GetSigmaDip(params,sigma_dip,sigma_restored_dip)
    dipposze = xy_to_ze2D(dip_pos,tipxfinal,tipyfinal,ramp_acute,ramp_dir) !dip positions in zeta, eta
    dipsze = (ramp_dir*dips+ramp_acute) !dips in ze coord. system
    sigma_dip_array = sigma_dip*pi/180.
    call trishear_func_dip_decol(dipposze,dipsze,v0,m,total_slip,increment,-slip_sense*PoverS,s,decolze,ramp_acute,sigma_dip_array)
    RestoredDips = ramp_dir*(dipsze - ramp_acute) !convert back to world coordinate system
    RestoredDips = RestoredDips*180/pi !Convert back to degrees
    !Fix dips > 90 deg or < -90 deg
    where(RestoredDips > 90.) RestoredDips = RestoredDips - 180. !Fix dips > 90 deg
    where(RestoredDips < -90.) RestoredDips = RestoredDips + 180. !Fix dips < - 90 deg
    call calc_dip_err(RestoredDips,outputs(2),sigma_dip_array*180./pi,sigma_restored_dip,restored_fit_params)
end if
if (DatTypes(DatType_terr) == 1) then !Terraces
    sin_ramp_acute = sin(ramp_acute)
    cos_ramp_acute = cos(ramp_acute)
    call GetSigmaTerr(params,sigma_terr,sigma_orig_elev)
    do n = 1,nterr !Go through each terrace
        bedlength = terraces(n)%npts !I'm keeping the "bed" variable names since there's no point creating new ones. They do the same ting either way.
        allocate(thisbed(2,bedlength))
        allocate(thisbedze(2,bedlength))
        allocate(Cterr(2,2,bedlength))
        Cterr = 0 !Initialize covariance matrix for terrace points (x,y).
        Cterr(1,1,:) = sigma_terr(n)**2
        Cterr(2,2,:) = sigma_terr(n)**2
        do i = 1,bedlength !Reflect, then rotate.
            Cterr(:,:,i) = matmul(matmul(Refl,Cterr(:,:,i)),Refl) !No need to transpose Refl, because transpose(Refl) = Refl.
            Cterr(:,:,i) = matmul(matmul(transpose(R),Cterr(:,:,i)),R) !Change Cterr into the rotated (zeta,eta) coordinate system. I think this is only necessary if sigma_x /= sigma_y to start or if there is covariance to start out.
        end do
        thisbed(:,:) = terraces(n)%pts
        thisbedze = xy_to_ze2D(thisbed,tipxfinal,tipyfinal,ramp_acute,ramp_dir); !Convert to eta, zeta
        call trishear_func_decol(thisbedze,v0,m,terr_slip(n),increment,-slip_sense*PoverS,s,decolze,ramp_acute,Cterr) !Run the trishear on this terrace
        thisbed = ze_to_xy2D(thisbedze,tipxfinal,tipyfinal,ramp_acute,ramp_dir) !Convert back to regular coordinate system.
        do i = 1,bedlength !Rotate, then reflect.
            Cterr(:,:,i) = matmul(matmul(R,Cterr(:,:,i)),transpose(R)) !Rotate the covariance matrix from (zeta,eta) back into (x,y) coordinates.
            Cterr(:,:,i) = matmul(matmul(Refl,Cterr(:,:,i)),Refl) !No need to transpose Refl, because transpose(Refl) = Refl.
        end do
        call calc_terr_err(n,bedlength,thisbed,terr_output_indiv,Cterr,sigma_orig_elev)
        deallocate(thisbed,thisbedze,Cterr)
    end do
    call calc_beds_err(terr_output_indiv,outputs(3)) 
end if
if (DatTypes(DatType_faultpts) == 1) then !Fault Points
    call GetSigmaFault(params,sigma_fault)
    do n=1,nfaultpts
        if (faultpts(2,n)>-(faultpts(1,n)-decolxy(1,1))*ramp_dir*tan((pi-ramp_acute)/2.)+decolxy(2,1)) then !Closest to ramp
            fault_errs(n) = abs(((faultpts(1,n)-decolxy(1,1))*tan_ramp-faultpts(2,n)+decolxy(2,1))/sqrt(1+tan_ramp**2))
        else !Closest to flat
            fault_errs(n) = abs(faultpts(2,n)-decolxy(2,1))
        end if
        if (faultpts(2,n)>tipyfinal) then
            if (sqrt((faultpts(2,n)-tipyfinal)**2+(faultpts(1,n)-tipxfinal)**2)>fault_errs(n)) then
                fault_errs(n) = sqrt((faultpts(2,n)-tipyfinal)**2+(faultpts(1,n)-tipxfinal)**2)
            end if
        end if
    end do
    call calc_fault_pts_err(fault_errs,outputs(4),sigma_fault)
end if
!Get the appropriate output into the single "output" variable
output = sum(outputs) !For chi squared sum. For probabilities, sum because output is ln(p). For RMS, there should only be one output.
exit
end do !end the do loop that lets us exit early
if (goodrun .eqv. .false.) then
    if (ResultType == 1 .or. ResultType == 3) then !RMS or chisq
        output = NaN
    else !Probability
        output = -huge(0.) !Because these models have 0 probability of being correct.
    end if
end if
end subroutine run_model_decol


subroutine trishear_func_decol(ptsze,v0,m,total_slip,increment,PoverS,s,decolze,ramp_angle,Cbed)
use constants
use analytic_trishear_module, only: analytic_trishear
use options, only: propagate
use err_and_prob, only: prop_bend_error
implicit none
double precision, dimension(1:,1:), intent(inout) :: ptsze
double precision, dimension(2,1), intent(in) :: decolze
double precision, intent(in) :: v0,m,total_slip,increment,PoverS,s,ramp_angle !Note: ramp_angle should be acute
double precision,dimension(1:,1:,1:),intent(inout) :: Cbed !Covariance matrices for x and y for all the points in the bed.
double precision,dimension(2,2) :: C,J !Covariance matrix for the specific point being worked on and Jacobian matrix for motion in the trishear velocity field.
integer :: k
double precision :: half_v0,xexp,yexp,slip,x,y,signy,decol_ang,xd,bisect_ang
double precision :: total_slip_abs
integer :: slip_sign,prop_sign
double precision :: tan_bisect,flatvx,flatvy
double precision :: tan_ramp,sin_ramp,cos_ramp
double precision :: hw_slip,fl_slip,fw_slip,rem_slip,slip2fl,slip2tri,tri_slip
double precision :: vx, vy
integer :: loc !tells which zone a point is in.

decol_ang = max(ramp_angle,pi-ramp_angle) !decollement angle from ramp
bisect_ang = 0.5*decol_ang !bisector angle
tan_bisect = tan(bisect_ang) !Just calculate once = faster
tan_ramp = tan(ramp_angle) !Just calculate once
sin_ramp = sin(ramp_angle)
cos_ramp = cos(ramp_angle)
half_v0 = v0/2; !v0/2
flatvx = v0*cos(ramp_angle)!vx in the flat hanging wall
flatvy = -v0*sin(ramp_angle) !vy in the flat hanging wall
xexp = 1/s; !exponent in vx term
yexp = (1+s)/s; !exponent in vy term
total_slip_abs = abs(total_slip)
if (v0 >= 0) then !determine the sign of the slip
    slip_sign = 1
else
    slip_sign = -1
end if
if (PoverS >= 0) then !determine the sign of P/S
    prop_sign = 1
else
    prop_sign = -1
end if
do k = 1,size(ptsze,2)
    x = ptsze(1,k)
    y = ptsze(2,k)
    slip = 0 !slip so far
    loc = 0 !location: ramp = 1, flat = 2, tri = 3, fw = 4    
    C = Cbed(:,:,k) !Covariance matrix for x and y for this point.
    do while (abs(slip)<total_slip_abs)
        rem_slip = total_slip - slip !remaining slip
        xd = decolze(1,1)-PoverS*slip !x for decollement start            
        if (loc == 0) then !need to find out location
            if (y>0 .and. y>=m*x .and. y<=tan_bisect*(x-xd)) then !ramp hanging wall
                loc = 1
            else if (y>=tan_ramp*(xd-x) .and. y>=m*x .and. y>0) then  !flat hanging wall
                loc = 2
            else if (x<xd .or. y<=-m*x) then !foot wall
                loc = 4
            else !trishear zone
                loc = 3
            end if
        end if
        if (loc == 1) then !ramp hanging wall
            if (slip_sign == prop_sign) then !Forward motion
                if (y <= m*(x+rem_slip*(1-PoverS))) then !goes into trishear zone
                    hw_slip = (y/m-x)/(1-PoverS)
                    loc = 3
                else !stays in ramp hw
                    hw_slip = rem_slip
                end if
            else !Backward motion / thrust inversion
                slip2fl = y/tan_bisect-x+xd !Slip required to reach the hanging wall flat
                slip2tri = (y/m-x)/(1-PoverS) !Slip required to reach the trishear zone
                if (abs(PoverS) > 1. .and. slip_sign*slip2tri < slip_sign*slip2fl .and. slip_sign*slip2tri &
                    < slip_sign*rem_slip) then !Goes into trishear zone.
                    hw_slip = slip2tri
                    loc = 3! Go to trishear zone
                else if (slip_sign*slip2fl < slip_sign*rem_slip) then !Goes into flat
                    hw_slip = slip2fl
                    call prop_bend_error(C,(/v0,0.0d0/),(/v0*cos(ramp_angle),-v0*sin(ramp_angle)/),(/0.0d0,0.0d0/),tan_bisect)
                    loc = 2 !Go to flat
                else !Stays in ramp hanging wall
                    hw_slip = rem_slip
                end if
            end if
            x = x+hw_slip-hw_slip*PoverS
            slip = slip + hw_slip
        else if (loc == 2) then !flat hanging wall
            if (y>tan_bisect*(x+rem_slip*cos_ramp-xd)+rem_slip*sin_ramp) then !stays on flat
                fl_slip = rem_slip
                x = x+fl_slip*cos_ramp-fl_slip*PoverS
                y = y-fl_slip*sin_ramp
            else !goes into hanging wall ramp
                fl_slip = slip_sign*(y-tan_bisect*(x-xd))/(sin_ramp+cos_ramp*tan_bisect);
                x = x+fl_slip*cos_ramp-fl_slip*PoverS
                y = tan_bisect*(x-(xd-fl_slip*PoverS))
                call prop_bend_error(C,(/v0*cos(ramp_angle),-v0*sin(ramp_angle)/),(/v0,0.0d0/),(/0.0d0,0.0d0/),tan_bisect)
                loc = 1
            end if
            slip = slip+fl_slip;
        else if (loc == 4) then !foot wall               
            if (slip_sign == prop_sign) then !forward motion, can't leave footwall
                fw_slip = rem_slip !slip while in footwall
            else !inversion, can leave fw when trishear zone comes back to it
                if (x-rem_slip*PoverS>0 .and. y > -m*(x-rem_slip*PoverS)) then !enters trishear zone
                    fw_slip = (x+y/m)/PoverS
                    loc = 3
                else
                    fw_slip = rem_slip
                end if
            end if
            x = x-fw_slip*PoverS
            slip = slip+fw_slip
        else !trishear zone, loc = 3
            if (s==1) then !Use the analytic solution for s==1
                call analytic_trishear(x,y,C,m,PoverS,rem_slip,slip_sign,tri_slip)
                slip = slip+tri_slip
                loc = 0;
            else
                if (y >= 0) then
                    signy = 1
                else
                    signy = -1
                end if
                slip = slip+increment
                vx = half_v0*((signy*(abs(y)/(m*x))**xexp)+1)-PoverS*increment
                vy = half_v0*(m/(1+s))*(((abs(y)/(m*x))**yexp)-1)
                J(1,1) = 1.-(increment/(2*s*x))*signy*(abs(y)/(m*x))**xexp !1+dvx/dx
                J(1,2) = (increment/(2*m*s*x))*(abs(y)/(m*x))**((1.-s)/s) !dvx/dy
                J(2,1) = -((increment*m)/(2*s*x))*(abs(y)/(m*x))**yexp !dvy/dx
                J(2,2) = 1.+(increment/(2*s*x))*signy*(abs(y)/(m*x))**xexp !1+dvy/dy
                C = matmul(matmul(J,C),transpose(J))
                x = x+vx
                y = y+vy
                loc = 0;
            end if
        end if
    end do
    ptsze(1,k) = x+PoverS*total_slip
    ptsze(2,k) = y
    if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
        Cbed(:,:,k) = C
    end if
end do
end subroutine trishear_func_decol

subroutine trishear_func_dip_decol(dipposze,dipsze,v0,m,total_slip,increment,PoverS,s,decolze,ramp_angle,sigma_dip_array)
use data_module, only : ndips
use constants
use analytic_trishear_module, only: analytic_trishear_dip
use options, only: propagate
implicit none
double precision, dimension(1:,1:), intent(inout) :: dipposze
double precision, dimension(1:),intent(inout) :: dipsze
double precision, dimension(1:,1:), intent(in) :: decolze
double precision, intent(in) :: v0,m,total_slip,increment,PoverS,s,ramp_angle !Note: Ramp angle should be acute.
double precision,dimension(ndips),intent(inout) :: sigma_dip_array !Array of sigma dip
integer :: j
double precision :: half_v0,xexp,yexp,slip,x,y,dip,signy,decol_ang,xd,bisect_ang
double precision :: sigmadip !The uncertainty in the dip being moved.
double precision :: v0term,exp1,exp2,sqrtm
double precision :: tan_bisect,flatvx,flatvy
double precision :: tan_ramp,sin_ramp,cos_ramp
double precision :: total_slip_abs
integer :: slip_sign,prop_sign
double precision :: theta,phi,gamma1,gamma2,delta
double precision :: hw_slip,fl_slip,fw_slip,rem_slip,slip2fl,slip2tri,tri_slip
double precision :: vx, vy
integer :: loc !tells which zone a point is in.

decol_ang = max(ramp_angle,pi-ramp_angle) !decollement angle from ramp
bisect_ang = 0.5*decol_ang !bisector angle
tan_bisect = tan(bisect_ang) !Just calculate once = faster
tan_ramp = tan(ramp_angle) !Just calculate once
sin_ramp = sin(ramp_angle)
cos_ramp = cos(ramp_angle)
half_v0 = v0/2; !v0/2
v0term = half_v0/s !v0/(2s) constant term that goes in front in dip equation
exp1 = (1+s)/(2*s) !exponent in 1st term for change in dip
exp2 = (1-s)/(2*s) !exponent in 2nd term for change in dip
flatvx = v0*cos(ramp_angle)!vx in the flat hanging wall
flatvy = -v0*sin(ramp_angle) !vy in the flat hanging wall
xexp = 1/s; !exponent in vx term
yexp = (1+s)/s; !exponent in vy term
total_slip_abs = abs(total_slip)
sqrtm = sqrt(m) !square root of m
if (v0 >= 0) then !determine the sign of the slip
    slip_sign = 1
else
    slip_sign = -1
end if
if (PoverS >= 0) then !determine the sign of P/S
    prop_sign = 1
else
    prop_sign = -1
end if
do j = 1,ndips
    x = dipposze(1,j)
    y = dipposze(2,j)
    dip = dipsze(j)
    sigmadip = sigma_dip_array(j)
    slip = 0 !slip so far
    loc = 0 !location: ramp = 1, flat = 2, tri = 3, fw = 4
    do while (abs(slip)<total_slip_abs)
        rem_slip = total_slip - slip !remaining slip
        xd = decolze(1,1)-PoverS*slip !x for decollement start
        if (loc == 0) then !need to find out location
            if (y>0 .and. y>=m*x .and. y<=tan_bisect*(x-xd)) then !ramp hanging wall
                loc = 1
            else if (y>=tan_ramp*(xd-x) .and. y>=m*x .and. y>0) then  !flat hanging wall
                loc = 2
            else if (x<xd .or. y<=-m*x) then !foot wall
                loc = 4
            else !trishear zone
                loc = 3
            end if
        end if
        if (loc == 1) then !ramp hanging wall
            if (slip_sign == prop_sign) then !Forward motion
                if (y <= m*(x+rem_slip*(1-PoverS))) then !goes into trishear zone
                    hw_slip = (y/m-x)/(1-PoverS)
                    loc = 3
                else !stays in ramp hw
                    hw_slip = rem_slip
                end if
            else !Backward motion / thrust inversion
                slip2fl = y/tan_bisect-x+xd !Slip required to reach the hanging wall flat
                slip2tri = (y/m-x)/(1-PoverS) !Slip required to reach the trishear zone
                if (abs(PoverS) > 1. .and. slip_sign*slip2tri < slip_sign*slip2fl .and. slip_sign*slip2tri &
                    < slip_sign*rem_slip) then !Goes into trishear zone.
                    hw_slip = slip2tri
                    loc = 3! Go to trishear zone
                else if (slip_sign*slip2fl < slip_sign*rem_slip) then !Goes into flat
                    hw_slip = slip2fl
                    loc = 2 !Go to flat
                    !Solve for change in dip assuming slip conservation across syncline
                    phi = ramp_angle !don't really need 2 names for the same thing   
                    theta = -dip !Since going backwards, theta is dip with ramp as horizontal, but negative since phi is other way now
                    gamma1 = pi/2-phi/2-theta
                    gamma2 = atan((tan(phi/2)-sin(theta)/(cos(phi/2)*cos(phi/2+theta)))**(-1)) ! Precalculating at least tan(phi/2) and cos(phi/2) would help.
                    if (gamma2 < 0) then !gamma2 should be in the range [0,pi]
                        gamma2 = pi+gamma2
                    end if
                    delta = pi-gamma1-gamma2 !change in dip
                    dip = dip+delta;
                    sigmadip = abs(1.-4.*sin(theta)*(-sin(theta)+sin(theta+phi))/ &
                        (-3.+2*cos(2*theta)+cos(phi)-2*cos(2*theta+phi)))*sigmadip
                else !Stays in ramp hanging wall
                    hw_slip = rem_slip
                end if
                x = x+hw_slip-hw_slip*PoverS
                slip = slip + hw_slip
            end if
        else if (loc == 2) then !flat hanging wall
            if (y>tan_bisect*(x+rem_slip*cos_ramp-xd)+rem_slip*sin_ramp) then !stays on flat
                fl_slip = rem_slip
                x = x+fl_slip*cos_ramp-fl_slip*PoverS
                y = y-fl_slip*sin_ramp
            else !goes into hanging wall ramp
                fl_slip = slip_sign*(y-tan_bisect*(x-xd))/(sin_ramp+cos_ramp*tan_bisect);
                x = x+fl_slip*cos_ramp-fl_slip*PoverS
                y = tan_bisect*(x-(xd-fl_slip*PoverS))
                loc = 1
                !Solve for change in dip assuming slip conservation across syncline
                phi = ramp_angle !don't really need 2 names for the same thing   
                theta = dip-ramp_angle! !Just rotate back to xy coord. system.
                gamma1 = pi/2-phi/2-theta
                gamma2 = atan((tan(phi/2)-sin(theta)/(cos(phi/2)*cos(phi/2+theta)))**(-1)) ! Precalculating at least tan(phi/2) and cos(phi/2) would help.
                if (gamma2 < 0) then !gamma2 should be in the range [0,pi]
                    gamma2 = pi+gamma2
                end if
                delta = pi-gamma1-gamma2 !change in dip
                dip = dip-delta;
                sigmadip = abs(1.+4.*sin(theta)*(-sin(theta)+sin(theta+phi))/ &
                    (-3.+2*cos(2*theta)+cos(phi)-2*cos(2*theta+phi)))*sigmadip
            end if
            slip = slip+fl_slip;
        else if (loc == 4) then !foot wall               
            if (slip_sign == prop_sign) then !forward motion, can't leave footwall
                fw_slip = rem_slip !slip while in footwall
            else !inversion, can leave fw when trishear zone comes back to it
                if (x-rem_slip*PoverS>0 .and. y > -m*(x-rem_slip*PoverS)) then !enters trishear zone
                    fw_slip = (x+y/m)/PoverS
                    loc = 3
                else
                    fw_slip = rem_slip
                end if
            end if
            x = x-fw_slip*PoverS
            slip = slip+fw_slip
        else !trishear zone, loc = 3
            if (s==1) then
                call analytic_trishear_dip(x,y,dip,sigmadip,m,PoverS,rem_slip,slip_sign,tri_slip)
                slip = slip+tri_slip
                loc = 0;
            else
                if (y >= 0) then
                    signy = 1
                else
                    signy = -1
                end if
                slip = slip+increment
                dip = dip+(v0term/x)*(signy*sqrtm*((abs(y)/(m*x))**(exp1))*cos(dip)+(1/sqrtm)*((abs(y)/(m*x))**(exp2))*sin(dip))**2
                sigmadip = sigmadip*abs(1.+signy*(increment/(s*x))*((abs(y)/(m*x))**xexp)*(0.5*((x/y)-(y/x))*sin(2*dip)+cos(2*dip)))
                vx = half_v0*((signy*(abs(y)/(m*x))**xexp)+1)-PoverS*increment
                vy = half_v0*(m/(1+s))*(((abs(y)/(m*x))**yexp)-1)
                x = x+vx
                y = y+vy
                loc = 0;
            end if
        end if
    end do
    dipposze(1,j) = x+PoverS*total_slip
    dipposze(2,j) = y
    dipsze(j) = dip
    if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
        sigma_dip_array(j) = sigmadip
    end if
end do
end subroutine trishear_func_dip_decol

end module trishear_decol