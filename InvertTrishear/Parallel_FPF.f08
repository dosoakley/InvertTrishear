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

module parallel_FPF
!Parallel fault propagation fold, kink-band model, using equations of Suppe and Medwedeff (1990).
integer :: bend_at_tipy !Tells whether decollement depth and tipy must be the same (1) or not (0)
integer :: tipy_can_be_below_bend !Tells whether tipy (initial) can be below the bend

contains

subroutine parallel_FPF_options
    use options, only: options_id
    use parameters, only: nparams
        implicit none
    do
        print*,'Must depth of fault bend be at initial fault tip depth? (0/1)'
        print*,'Choosing 0 may not work properly.'
        read(options_id,*) bend_at_tipy
        if (bend_at_tipy == 0) then !Note: I'm not sure if this actually work properly if depth of fault bend is not at initial fault tip depth. I would have to look into Suppe and Medwedeff's equations more.
            nparams = 6
            exit
        else if(bend_at_tipy == 1) then
            nparams = 5
            exit
        else
            print*,'Invalid Command'
        end if
    end do
    tipy_can_be_below_bend = 0
    !For now, we're requiring the tip to start above the bend. It might be possible to have it start below the bend, but that would create a more complicated geometry.
end subroutine parallel_FPF_options

subroutine run_model_pFPF(output,params)
use constants
use parameters, only: nparams,v0,slip_sense,nparam_start_restored_fit,nrestored_fit_params,nparam_start_growth,nparam_start_terr
use options, only: TipToSolve,ndatatypes_poss,DatTypes,ResultType,DatType_beds,DatType_dips,DatType_faultpts,DatType_terr,FitType,&
    terr_age_order,DatType_beds_restored,beds_age_order,beds_age_order,limit_bed_dips,min_bed_slopes,max_bed_slopes,FitType
use math, only: csc,sec,init_random_seed
use data_module
use err_and_prob
use data_uncertainties, only: GetSigmaBed,GetSigmaDip,GetSigmaTerr,GetSigmaFault,GetSigmaBedRestored,GetLc,GetSigma2Bed
use growth_strata, only: GetGrowthSlip
use beds_order, only: check_beds_order
implicit none
double precision,intent(out) :: output !the output result (RMS or probability)
double precision,dimension(ndatatypes_poss) :: outputs !A single vector holding outputs for beds, dips, and terraces (in that order).
double precision,dimension(nparams),intent(inout) :: params !values of the parameters
double precision :: tipx,tipy,total_slip,total_slip_mag !parameters
double precision :: slip !Slip to restore a bed, may be growth or pregrowth.
double precision :: PoverS !P/S ratio
double precision,dimension(nterr) :: terr_slip !Slip needed to restore each terrace.
double precision,dimension(ngrowth) :: growth_slip !Slip needed to restore each growth strata bed.
double precision,dimension(nrestored_fit_params) :: restored_fit_params
double precision,dimension(2) :: ramp_angle, ramp_angle_deg,ramp_acute !1st value is higher segment, 2nd is lower segment
double precision :: phi !Difference between the two ramp angles.
double precision,dimension(2,1) :: bendxy !x,y coordinates of where the decollement meets the ramp
double precision :: gamma1, gamma !Backlimb axial angle and forelimb axial angle
double precision :: gamma_star,beta1,beta2,kappa !Angles within the fold used in the calculation.
double precision :: R0 !Ratio of slip on the upper segment to slip on the lower segment.
double precision :: dist !The distance between (tipx,tipy) and the bend
double precision :: tipxinit, tipyinit, tipxfinal, tipyfinal !Initial and final (in terms of forward motion) tip positions
integer :: tip_seg,tip_seg_init,tip_seg_final !tells which segment of the fault the fault tip is in. 1=upper segment; 2=lower segment
integer :: ramp_dir !dip direction of ramp
logical :: goodrun !tells if the run is good
double precision :: NaN = 0 !For creating NaN by dividing by zero
integer :: bedlength !number of points in the bed being worked on
double precision,dimension(:,:),allocatable :: thisbed !the bed being worked on
double precision :: slope,b !Slope and intercept of a bed
double precision :: slope_last, b_last !Slope and intercept of the previous bed restored.
double precision,dimension(2) :: bed_xlims,bed_xlims_last !x limits of current and previous bed restored.
double precision,dimension(2,ndips) :: dip_pos_temp !dip positions array that I can change.
double precision,dimension(ndips) :: dips_temp !dips array that I can change
double precision,dimension(ndips) :: RestoredDips !values of dips after inverse fault movement is complete
double precision,dimension(nbeds) :: beds_output_indiv !Output from each of the beds individually.
double precision,dimension(n_restored_beds) :: restored_beds_output_indiv !Output from each of the beds individually.
double precision,dimension(nterr) :: terr_output_indiv !Output from each of the terraces individually.
double precision,dimension(nfaultpts) :: fault_errs !Errors from each fault point.
double precision,dimension(ndips) :: sigma_dip_array !Array of sigmas for dips
double precision,dimension(:,:,:),allocatable :: Cbed,Cterr !Covariance matrix for uncertainties in bed point positions and terrace point positions.
double precision,dimension(:,:,:),allocatable :: Cbed2 !Covariance matrix for second (correlated) uncerrtainties in bed point positions if ResultType is 5.
double precision,dimension(nbeds) :: sigma_bed_x,sigma_bed_y !Uncertainty in bed points x and y coordinates before any propagation.
double precision,dimension(nbeds) :: sigma2_bed_x,sigma2_bed_y !Uncertainty in bed points x and y coordinates for second (correlated) uncertainty if ResultType is 5.
double precision :: Lc !Correlation length for bed data if correlated.
double precision :: sigma_bed_restored_x,sigma_bed_restored_y !Uncertainty in x and y coordinates of points on the restored beds, when these are used as independent data.
double precision :: sigma_dip,sigma_restored_dip !Uncertainty in dips in the deformed and restored states.
double precision,dimension(nterr) :: sigma_terr,sigma_orig_elev !Uncertainty in terrace point positions and in the original elevation of the terrace inner edge.
double precision :: sigma_fault !Uncertainty in fault points.
double precision :: f !Result of function evaluation, used in Newton-Raphson parts.
double precision :: prec_limit !Precision limit that f must be lower than for the Newton-Raphson evaluations to stop.
integer :: niter,max_iter !Number of iterations performed and maximum number to use in Newton-Raphson method.
integer :: i,n,g ! counters, n is for beds; g is for growth strata.
double precision :: rn !Holds a random number
double precision,dimension(2) :: bounds,fbounds !Bounding values and corresponding values of the objective function f used in the bisection method.
NaN = NaN/NaN
goodrun = .true. !Initialize it to true. So far, this is a good run.
do !This do loop lets us exit anytime it turns out we have an impossible combination of parameters
!Make assignments
tipx = params(1)
tipy = params(2) !Note: The model of Suppe and Medwedeff (1990) seems to assume tipy at bendy. I will have to check whether or not that has to be the case.
total_slip_mag = params(3)
ramp_angle_deg(1) = params(4)
ramp_angle_deg(2) = params(5)
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
if (ramp_angle(1)<=pi/2 .and. ramp_angle(2)<=pi/2) then !1=dips left, -1 = dips right
    ramp_dir = 1
else if(ramp_angle(1)>pi/2 .and. ramp_angle(2)>pi/2) then
    ramp_dir = -1
else !This should only occur if ramp_dir is different for the two ramp angles.
    goodrun = .false.
    exit
end if
do i = 1,2
    ramp_acute(i) = min(ramp_angle(i),pi-ramp_angle(i)) !Acute version of ramp angle
end do
if (ramp_acute(1) < ramp_acute(2)) then !The fault must get steeper, not shallower, as it gets higher
    goodrun = .false.
    exit
end if
phi = ramp_acute(1)-ramp_acute(2) !Phi should be acute as well.
!Calculate gamma angles. I can't figure out how to solve the gamma1 or gamma_star equations analytically. (No answer from Mathematica.) So I'll have to use Newton-Raphson.
prec_limit = 1e-5 !Necessary precision at which we stop trying to improve the value of gamma1 or gamma star. Possibly allow the user to choose, or base it on number of degrees. e.g. difference between tan(phi) and tan(phi+1deg).
max_iter = 1e3 !Maximum number of iterations after which we stop searching with Newton_Raphson.
!Solve for gamma1 with Suppe and Medwedeff eq. 9
!gamma1 should be >pi/4, or the fold will have an overturned backlimb
!call init_random_seed()
niter = 0 !Number of iterations
call random_number(rn)
gamma1 = rn*pi/4+pi/4 !Initial random gamma1 in the range [pi/4,pi/2). This should avoid any values on the wrong side of the asymptote.
f = (-sin(gamma1 + ramp_acute(2))*(sin(2*gamma1 + ramp_acute(2)) + sin(ramp_acute(2))))/(cos(gamma1 + ramp_acute(2))* &
    (sin(2*gamma1 + ramp_acute(2)) + sin(ramp_acute(2))) - sin(gamma1)) - tan(phi)
do while (abs(f) > prec_limit .or. gamma1<0 .or. gamma1>pi/2) !Newton-Raphson
    gamma1 = gamma1 - ((cos(gamma1+phi)-cos(gamma1)*cos(2*gamma1+2*ramp_acute(2)+phi))*sec(phi)*(-sin(gamma1)+cos(gamma1) &
        *sin(2*(gamma1+ramp_acute(2)))))/(cos(2*gamma1)-cos(2*ramp_acute(2))+4*sin(gamma1+ramp_acute(2))**2)
    if (gamma1 > pi/2) then !Keep gamma1 in the first quadrant.
        gamma1 = pi-gamma1
    else if (gamma1 < 0) then
        gamma1 = -gamma1
    end if
    f = (-sin(gamma1 + ramp_acute(2))*(sin(2*gamma1 + ramp_acute(2)) + sin(ramp_acute(2))))/(cos(gamma1 + ramp_acute(2))* &
        (sin(2*gamma1 + ramp_acute(2)) + sin(ramp_acute(2))) - sin(gamma1)) - tan(phi)
    niter = niter+1
    if (niter == max_iter) then
        print*,'Error: No solution for gamma1 found within desired precision.'
        print*,'gamma1 = ',gamma1
        print*,'error = ',f
        goodrun = .false.
        exit
    end if
end do
!Calculate beta1
beta1 = ramp_acute(1)-(pi-2.*gamma1)
!Solve for gamma_star with Suppe and Medwedeff eq. 6
!Note: There are 2 possible solutions for gamma_star, but one is gamma_star = gamma1, which means gamma = pi/2 and 2*gamma = pi, so there is no fold. We want to be sure to avoid that solution.
!Newton-Raphson too often gives us a value close to the gamma_star = gamma1 solution (but just below it so that gamma_star < gamma1 still). I'm not sure how best to avoid this.
!Therefore, I'm just using the bisection method. It tends to be a little slower, but it should avoid the problems that Newton-Raphson was causing.
niter = 0 !Number of iterations
!Start with beta1 and gamma1. gammastar must be less than gamma1 b/c gamma = pi/2 + gamma_star-gamma1 and gamma must be < pi/2. gamma_star must be > beta1 so distance ef and P/S are positive.
bounds(1) = beta1
bounds(2) = gamma1
fbounds = (sin(bounds)*sin(gamma1-beta1))/(sin(gamma1-bounds)+(sin(gamma1)*sin(bounds-beta1))/ &
    sin(2*bounds-beta1))-sin(ramp_acute(1)) !gamma_star = 0 should give a negative f, while the other bound, once less than gamma1 (for which it's 0) should give a positive one.
do
    gamma_star = bounds(1)+0.5*(bounds(2)-bounds(1)) !Take the bisector of the two bounds.
    f = (sin(gamma_star)*sin(gamma1-beta1))/(sin(gamma1-gamma_star)+(sin(gamma1)*sin(gamma_star-beta1))/ &
        sin(2*gamma_star-beta1))-sin(ramp_acute(1))
    if (abs(f) < prec_limit) then !If we've found a gamma_star that is good enough to proceed with.
        exit
    else
        if (f<0) then !Replace the lower bound
            bounds(1) = gamma_star
            fbounds(1) = f
        else !Replace the upper bound
            bounds(2) = gamma_star
            fbounds(2) = f
        end if
    end if
    niter = niter+1
    if (niter == max_iter) then
        print*,'Error: No solution for gamma_star found within desired precision.'
        print*,'gamma_star = ',gamma_star
        print*,'error = ',f
        goodrun = .false.
        exit
    end if
end do
if (goodrun .eqv. .false.) then !If the above numerical solutions didn't work out, then discard the run.
    exit !Leave the loop now before wasting any more time or risking issues with later parts of the run.
end if
!There are cases where there is no 0 between beta1 and gamma1. (There is another 0 at a values > gamma1, but that's not allowed.) So really, those are impossible models. We should discard them.
!One way would be to test derivative of f at gamma_star = gamma. If it's positive (f decreases as gamma_star decreases), there is no 0 between beta1 and gamma1, and the model is impossible.
!if (gamma1-gamma_star<0.01) then
!    print*,gamma1,gamma_star
!end if
!Calculate gamma from gamma_star and gamma1
gamma = pi/2 + gamma_star-gamma1
!Calculate beta2
beta2 = pi-2.*gamma_star+beta1
!Calculate kappa: The slope of the line between the fault tip and point e.
kappa = pi - beta2 + ramp_acute(1)
!Calculate P/S ratio
PoverS = 1./(1.-sin(ramp_acute(1))/sin(beta2)); !Note that this is the ratio of fault propagation to slip on the upper fault segment. So the total propagation is total_slip*R0*PoverS
R0 = sin(gamma1+ramp_acute(2))/sin(gamma1+ramp_acute(1));
if (PoverS < 0) then
    !This seems to occur when there are two possible solutions for gamma_star and we pick the smaller one, resulting in too small a sin(beta2).
    !It looks like the f(gamma_star) function has a minimum at beta1, and I should be finding a value of gamma_star between beta1 and gamma1.
    !This also makes sense from consideration of the fact that for distance ef to be positive, gamma_star must be > beta1 (see eqn. for ef_constant below).
    print*,'Error: P/S<0'
    print*,'gamma1 = ',gamma1,' gamma_star = ',gamma_star
    goodrun = .false.
    exit
end if
!Get or calculate y coordinate of the fault bend.
if (bend_at_tipy == 0) then
    bendxy(2,1) = params(6)
else
    if (TipToSolve == 1) then !tipy is initial tip
        bendxy(2,1) = tipy
    else !tipy = tipy is final tip
        bendxy(2,1) = tipy + R0*total_slip*PoverS*sin(ramp_angle(1)) !Note: This assumes that the final tip must be in the upper fault segment (both by using R0 and by using ramp_angle(1)).
    end if
end if
!Determine which segment the tip (tipx,tipy) is in. Note: this can be initial or final tip depending which one tipx,tipy represent.
if (tipy >= bendxy(2,1)) then
    tip_seg = 1
else
    tip_seg = 2
end if
!Calculate the decollement x coordinate
if (bend_at_tipy == 0) then
    bendxy(1,1) = tipx-(tipy-bendxy(2,1))/tan(ramp_angle(tip_seg)); !x coordinate of where decollement meets ramp
else
    if (TipToSolve == 1) then !tipx is initial tip
        bendxy(1,1) = tipx
    else !tipx is final tip
        bendxy(1,1) = tipx + R0*total_slip*PoverS*cos(ramp_angle(1)) !Note: This assumes that the final tip must be in the upper fault segment (both by using R0 and by using ramp_angle(1)).
    end if
end if
!Calculate initial and final tipx and tipy (final = present-day, initial = when trishear motion began). Final must be higher than initial.
if (TipToSolve == 1) then !If tipx, tipy are initial
    tipxinit = tipx
    tipyinit = tipy
    tip_seg_init = tip_seg
    !Note: If we allowed tipxfinal to be lower than tipxinit, this wouldn't be so simple. Same for if TipToSolve = 2 below
    if (tip_seg==1) then !Starts and stays in upper segment
        tipxfinal = tipx+R0*total_slip_mag*PoverS*cos(ramp_angle(1)) !Remember that if doing an inversion, total_slip and slip_sense will be negative
        tipyfinal = tipy+R0*total_slip_mag*PoverS*sin(ramp_angle(1))
        tip_seg_final = 1
    else
        dist = sqrt((tipy-bendxy(2,1))**2+(tipx-bendxy(1,1))**2) !Distance the point propagates to reach the bend.
        if (total_slip_mag*PoverS <= dist) then !Starts and stays in lower segment
            tipxfinal = tipx+total_slip_mag*PoverS*cos(ramp_angle(2))
            tipyfinal = tipy+total_slip_mag*PoverS*sin(ramp_angle(2))
            tip_seg_final = 2
        else !Goes from lower into upper segment
            tipxfinal = tipx+dist*cos(ramp_angle(2)) + R0*(total_slip_mag*PoverS-dist)*cos(ramp_angle(1))
            tipyfinal = tipy+dist*sin(ramp_angle(2)) + R0*(total_slip_mag*PoverS-dist)*sin(ramp_angle(1))
            tip_seg_final = 1
        end if
    end if
else !If tipx, tipy are final
    tipxfinal = tipx
    tipyfinal = tipy
    tip_seg_final = tip_seg
    if (tip_seg == 2) then !Starts and stays in lower segment
        tipxinit = tipx-total_slip_mag*PoverS*cos(ramp_angle(2)) !I'm not completely sure if P/S is the same if the tip is in the lower segment.
        tipyinit = tipy-total_slip_mag*PoverS*sin(ramp_angle(2))
        tip_seg_init = 2
    else
        dist = sqrt((tipy-bendxy(2,1))**2+(tipx-bendxy(1,1))**2) !Distance tip will have to propagate to change segments.
        if (R0*total_slip_mag*PoverS <= dist) then !Starts and stays in upper segment
            tipxinit = tipx-R0*total_slip_mag*PoverS*cos(ramp_angle(1))
            tipyinit = tipy-R0*total_slip_mag*PoverS*sin(ramp_angle(1))
            tip_seg_init = 1
        else !tip_init in lower segment, tip_final in upper segment
            tipxinit = tipx-dist*cos(ramp_angle(1)) - (total_slip_mag*PoverS-dist/R0)*cos(ramp_angle(2))
            tipyinit = tipy-dist*sin(ramp_angle(1)) - (total_slip_mag*PoverS-dist/R0)*sin(ramp_angle(2))
            tip_seg_init = 2
        end if
    end if
end if
!Check that tipy_init > bend_y if required
if (tipy_can_be_below_bend == 0 .and. tipyinit < bendxy(2,1)) then
    goodrun = .false.
    exit
end if
!Initialize the outputs array.
outputs = 0 !For RMS or chi squared it should be 0 anyway. For probabilities, p should be 1, so ln(p) = 0.
!Now restore data.
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
        allocate(Cbed(2,2,bedlength))
        Cbed = 0 !Initialize covariance matrix for bed points (x,y).
        Cbed(1,1,:) = sigma_bed_x(n)**2
        Cbed(2,2,:) = sigma_bed_y(n)**2
        thisbed(:,:) = beds(n)%pts !Still need to reflect if it's backwards.
        if (ramp_dir == -1) then !Reflect around the y axis, so that the fault now dips left. tipxfinal also gets reflected via the ramp_dir below.
            thisbed(1,:) = -thisbed(1,:)
            bendxy(1,:) = -bendxy(1,:)
            Cbed(1,2,:)=-Cbed(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cbed(2,1,:)=-Cbed(2,1,:)
        end if
        if (beds(n)%growth .eqv. .true.) then !For growth strata
            slip = growth_slip(g)
            g = g+1
        else
            slip = total_slip
        end if
        call PFPF_points(thisbed,v0,slip,PoverS,bendxy,ramp_acute,gamma,gamma1,gamma_star,kappa,beta1,beta2,&
            (/ramp_dir*tipxfinal,tipyfinal/),Cbed) !Run the model on this bed
        if (ramp_dir == -1) then !Reflect back.
            thisbed(1,:) = -thisbed(1,:)
            bendxy(1,:) = -bendxy(1,:)
            Cbed(1,2,:)=-Cbed(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cbed(2,1,:)=-Cbed(2,1,:)
        end if
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
        deallocate(thisbed,Cbed)
        if (ResultType==5) deallocate(Cbed2)
    end do
    call calc_beds_err(beds_output_indiv,outputs(1))
    if (DatTypes(DatType_beds_restored) == 1) then
        call calc_beds_err(restored_beds_output_indiv,outputs(5)) !Note: This won't work properly if the output is RMS error, but currently you're not allowed to use RMS error with more than one data type.
    end if
end if
if (goodrun .eqv. .false.) then
    exit !Leave the loop now before wasting any more time or risking issues with later parts of the run.
end if
if (DatTypes(DatType_dips) == 1) then !Dips
    call GetSigmaDip(params,sigma_dip,sigma_restored_dip)
    dip_pos_temp = dip_pos !Temporary dip positions for this inversion, so we don't mess with the original ones
    dips_temp = ramp_dir*dips !Rotate if necessary
    sigma_dip_array = sigma_dip*pi/180.
    if (ramp_dir == -1) then !Reflect around the y axis, so that the fault now dips left. tipxfinal also gets reflected via the ramp_dir below.
        dip_pos_temp(1,:) = -dip_pos_temp(1,:)
        bendxy(1,:) = -bendxy(1,:)
    end if
    call PFPF_dips(dip_pos_temp,dips_temp,v0,total_slip,PoverS,bendxy,ramp_acute,gamma,gamma1,gamma_star,kappa,beta1,&
        beta2,(/ramp_dir*tipxfinal,tipyfinal/),sigma_dip_array)
    if (ramp_dir == -1) then !Reflect back.
        dip_pos_temp(1,:) = -dip_pos_temp(1,:)
        bendxy(1,:) = -bendxy(1,:)
    end if
    RestoredDips = ramp_dir*(dips_temp) !convert back to world coordinate system
    RestoredDips = RestoredDips*180/pi !Convert back to degrees
    !Fix dips > 90 deg or < -90 deg
    where(RestoredDips > 90.) RestoredDips = RestoredDips - 180. !Fix dips > 90 deg
    where(RestoredDips < -90.) RestoredDips = RestoredDips + 180. !Fix dips < - 90 deg
    call calc_dip_err(RestoredDips,outputs(2),sigma_dip_array*180./pi,sigma_restored_dip,restored_fit_params)
end if
if (DatTypes(DatType_terr) == 1) then !Terraces
    call GetSigmaTerr(params,sigma_terr,sigma_orig_elev)
    do n = 1,nterr !Go through each terrace
        bedlength = terraces(n)%npts !I'm keeping the "bed" variable names since there's no point creating new ones. They do the same ting either way.
        allocate(thisbed(2,bedlength))
        allocate(Cterr(2,2,bedlength))
        Cterr = 0 !Initialize covariance matrix for bed points (x,y).
        Cterr(1,1,:) = sigma_terr(n)**2
        Cterr(2,2,:) = sigma_terr(n)**2
        thisbed(:,:) = terraces(n)%pts
        if (ramp_dir == -1) then !Reflect around the y axis, so that the fault now dips left. tipxfinal also gets reflected via the ramp_dir below.
            thisbed(1,:) = -thisbed(1,:)
            bendxy(1,:) = -bendxy(1,:)
            Cterr(1,2,:)=-Cterr(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cterr(2,1,:)=-Cterr(2,1,:)
        end if
        call PFPF_points(thisbed,v0,terr_slip(n),PoverS,bendxy,ramp_acute,gamma,gamma1,gamma_star,kappa,beta1,beta2,&
            (/ramp_dir*tipxfinal,tipyfinal/),Cterr) !Run the model on this terrace
        if (ramp_dir == -1) then !Reflect back.
            thisbed(1,:) = -thisbed(1,:)
            bendxy(1,:) = -bendxy(1,:)
            Cterr(1,2,:)=-Cterr(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cterr(2,1,:)=-Cterr(2,1,:)
        end if
        call calc_terr_err(n,bedlength,thisbed,terr_output_indiv,Cterr,sigma_orig_elev)
        deallocate(thisbed,Cterr)
    end do
    call calc_beds_err(terr_output_indiv,outputs(3))
end if
if (DatTypes(DatType_faultpts) == 1) then !Fault Points
    call GetSigmaFault(params,sigma_fault)
    do n=1,nfaultpts
        if (faultpts(2,n)>-(faultpts(1,n)-bendxy(1,1))*ramp_dir*tan((pi-ramp_acute(1)-ramp_acute(2))/2.)+bendxy(2,1)) then !Closest to upper ramp. Note: This slope might be different if we let the upper ramp be shallower than the lower.
            fault_errs(n) = abs(((faultpts(1,n)-bendxy(1,1))*tan(ramp_angle(1))-faultpts(2,n)+bendxy(2,1))/ &
                sqrt(1+tan(ramp_angle(1))**2))
        else !Closest to lower ramp
            fault_errs(n) = abs(((faultpts(1,n)-bendxy(1,1))*tan(ramp_angle(2))-faultpts(2,n)+bendxy(2,1))/ &
                sqrt(1+tan(ramp_angle(2))**2))
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
end subroutine run_model_pFPF

subroutine PFPF_points(pts,v0,total_slip,PoverS,bendxy,ramp_angle,gamma,gamma1,gamma_star,kappa,beta1,beta2, &
    tip_start,Cbed)
!Subroutine for moving points for a parallel fault propagation fold.
!Important notes: We're assuming that the beds start flat. If not, replace gamma, gamma1, and possibly some other angles with angles taking regional dip into account.
!Also, we're assuming the tip starts at the bend. If not, that will change things too.
!I'm not sure the best way to treat the uncertainty in this. Velocity within a domain isn't a function of position, so uncertainty won't change, but when a point reaches
!a boundary, uncertainty in position means there's an uncertainty in where it reaches the boundary and a correlated uncertainty in the amount of slip in each domain.
!Note: A lot of what is done in each zone is similar. Maybe I could write a more general function, applicable to any zone.
use options, only: propagate
use err_and_prob, only: prop_bend_error
implicit none
double precision,dimension(1:,1:), intent(inout) :: pts
double precision,dimension(2,1), intent(in) :: bendxy
double precision,intent(in) :: v0,total_slip,PoverS
double precision,dimension(2),intent(in) :: ramp_angle !Note: Ramp angle should be acute.
double precision,intent(in) :: gamma,gamma1,gamma_star,kappa,beta1,beta2
double precision,dimension(2),intent(in) :: tip_start
double precision,dimension(1:,1:,1:),intent(inout) :: Cbed !Covariance matrices for x and y for all the points in the bed.
double precision,dimension(2) :: bend !The position of the bend as a 1D array.
double precision :: ef_constant ! A constant used in calculating the position and velocity of the point e.
double precision :: slip_sign !Forwards = 1, backwards = -1
double precision,dimension(2) :: pte_start,pte !Point e. pte_start is where it is at the beginning of moving each point. pte can move as slip is restored for each point.
double precision :: R0, R1, R2 !Slip ratios.  All expresed relative to slip on the lowest segment.
double precision,dimension(2) :: v1,v2,v3,v4,v5 !Point velocities in each domain
double precision,dimension(2) :: vb,ve,vt !Boundary velocities (and velocities of the points: bend, point e, and tip that the boundary go through)
double precision :: gam1_slope,gam_slope,kap_slope !Boundary slopes (unchanging).
double precision,dimension(2) :: tip !The location of the fault tip. Moves as the fold is restored.
double precision,dimension(2) :: pt !The point being moved
integer :: k !A counter of points.
double precision,dimension(2,2) :: C !Covariance matrix for the specific point being worked on. Note: This isn't currently being used.
integer :: loc !tells which zone a point is in.
double precision :: slip,rem_slip !Slip so far and slip remaining
double precision :: slip1,slip2,slip3,slip4,slip5 !Slip amounts needed to reach different boundaries.
logical,dimension(5) :: loc_passed !Tells whether or not (true or false) a point has passed through each of the five domains.
bend = bendxy(:,1) !Only need this to be a 1D array.
ef_constant = sin(ramp_angle(1))*(sin(gamma_star-beta1)/(sin(gamma_star)*sin(2.*gamma_star-beta1))) !To be used in calculating the distance ef or the velocity of point e.
if (v0 >= 0) then !determine the sign of the slip
    slip_sign = 1
else
    slip_sign = -1
end if
!Calculate the slip ratios. All expresed relative to slip on the lowest segment.
R0 = sin(gamma1+ramp_angle(2))/sin(gamma1+ramp_angle(1));
R1 = R0*sin(gamma1+ramp_angle(1))/sin(gamma1+gamma);
R2=R0*sin(beta2)/sin(beta2-ramp_angle(1)+gamma);
!Calculate the starting position of point e.
if (slip_sign == 1) then !Forwards slip
    pte_start = tip_start !Note: This would not work for normal faults, if they were allowed.
else !Backwards slip
    pte_start = tip_start-PoverS*R0*total_slip*ef_constant*(/cos(kappa),sin(kappa)/) !minus because total_slip is negative.
end if
!Calculate velocities for each domain
v1 = v0*(/cos(ramp_angle(2)),sin(ramp_angle(2))/)
v2 = R0*v0*(/cos(ramp_angle(1)),sin(ramp_angle(1))/)
v3 = R1*v0*(/cos(gamma),sin(gamma)/)
v4 = R2*v0*(/cos(gamma),sin(gamma)/)
v5 = (/0,0/)
!Calculate velocities for each boundary
vb = (/0,0/) !Velocity of the bend: doesn't move.
ve = R0*v0*PoverS*(/cos(ramp_angle(1))+ef_constant*cos(kappa),sin(ramp_angle(1))+ef_constant*sin(kappa)/) !Velocity of point e.
vt = R0*v0*PoverS*(/cos(ramp_angle(1)),sin(ramp_angle(1))/) !Velocity of the fault tip.
!Calculate slopes for each boundary
gam1_slope = -tan(gamma1) !Slope of axes with dip gamma, dipping to the right.
gam_slope = tan(gamma) !Slope of axes with dip gamma1, dipping to the left.
kap_slope = tan(kappa) !Slope of boundary with dip kappa.
do k = 1,size(pts,2)
    tip = tip_start !Where the tip starts moving from. (Generally the final tip position for restoring a fault-propagation fold.)
    pt = pts(:,k) !The point to be moved.
    pte = pte_start !Location e in Figure 11.6  of Allmendinger et al. (2012) and Fig. 16 of Suppe and Medwedeff (1990)
    C = Cbed(:,:,k) !Covariance matrix for x and y uncertainty.
    slip = 0 !slip so far
    rem_slip = total_slip !Remaining slip
    loc = 0 !location, zone numbers from Figure 11.6, Allmendinger et al. (2012), 5 = footwall, 0 = unknown
    loc_passed = .false.
    do while (rem_slip*slip_sign>0.)
        tip = tip + (slip/v0)*vt
        pte = pte + (slip/v0)*ve
        if (loc == 1) then !Domain 1
            loc_passed(1) = .true.
            if (slip_sign == 1) then !Forward motion
                slip2 = calc_slip(pt,bend,v0,v1,vb,gam1_slope) !Slip needed to reach domain 2
                if ((.not. loc_passed(2)) .and. rem_slip > slip2) then !Goes into domain 2
                    slip = slip2
                    call prop_bend_error(C,v1,v2,vb,gam1_slope)
                    loc = 2
                else
                    slip = rem_slip !Stays in domain 1
                end if
            else !Backwards motion
                slip = rem_slip !Stays in domain 1
            end if
            pt = pt+(slip/v0)*v1
        else if (loc == 2) then !Domain 2
            loc_passed(2) = .true.
            if (slip_sign == 1) then !Forward motion
                slip3 = calc_slip(pt,pte,v0,v2,ve,gam1_slope) !Slip needed to reach domain 3
                if ((.not. loc_passed(3)) .and. rem_slip > slip3 .and. pt(2)+(slip3/v0)*v2(2) > pte(2)+(slip3/v0)*ve(2) &
                    .and. slip3>0) then !Goes into domain 3
                    slip = slip3
                    call prop_bend_error(C,v2,v3,ve,gam1_slope)
                    loc = 3
                else
                    slip4 = calc_slip(pt,tip,v0,v2,vt,kap_slope) !Slip needed to reach domain 4
                    if ((.not. loc_passed(4)) .and. rem_slip>slip4 .and. pt(2)+(slip4/v0)*v2(2) < pte(2)+(slip4/v0)*ve(2) &
                        .and. slip4>0)then !Goes into domain 4
                        slip = slip4
                        call prop_bend_error(C,v2,v4,vt,kap_slope)
                        loc = 4
                    else
                        slip = rem_slip !Stays in domain 2
                    end if
                end if
            else !Backward motion
                slip1 = calc_slip(pt,bend,v0,v2,vb,gam1_slope)
                slip3 = calc_slip(pt,pte,v0,v2,ve,gam1_slope)
                slip4 = calc_slip(pt,tip,v0,v2,vt,kap_slope)
                if ((.not. loc_passed(4)) .and. rem_slip < slip4 .and. slip4<0 .and. pt(2)+(slip4/v0)*v2(2) &
                    < pte(2)+(slip4/v0)*ve(2) .and. slip1<slip4) then !Goes into domain 4.
                    slip = slip4
                    call prop_bend_error(C,v2,v4,vt,kap_slope)
                    loc = 4
                else if ((.not. loc_passed(3)) .and. rem_slip < slip3 .and. slip3<0 .and. pt(2)+(slip3/v0)*v2(2) &
                    > pte(2)+(slip3/v0)*ve(2) .and. slip1<slip3) then !Goes into domain 3.
                    slip = slip3
                    call prop_bend_error(C,v2,v3,ve,gam1_slope)
                    loc = 3
                else if ((.not. loc_passed(1)) .and. rem_slip < slip1) then !Goes into domain 1. < because slips are negative
                    slip = slip1
                    call prop_bend_error(C,v2,v1,vb,gam1_slope)
                    loc = 1
                else
                    slip = rem_slip !Stays in domain 2
                end if
            end if
            pt = pt+(slip/v0)*v2
        else if (loc == 3) then !Domain 3
            loc_passed(3) = .true.
            if (slip_sign == 1) then !Forward motion
                slip2 = calc_slip(pt,pte,v0,v3,ve,gam1_slope) !Slip needed to reach domain 2
                slip4 = calc_slip(pt,pte,v0,v3,ve,gam_slope) !Slip needed to reach domain 4.
                if ((.not. loc_passed(2)) .and. rem_slip > slip2 .and. slip2 > 0 .and. &
                        (slip4>slip2 .or. slip4<0 .or. loc_passed(4))) then !Goes into domain 2
                    slip = slip2
                    call prop_bend_error(C,v3,v2,ve,gam1_slope)
                    loc = 2
                else if ((.not. loc_passed(4)) .and. rem_slip > slip4 .and. slip4>0) then !Goes into domain 4
                    slip = slip4;
                    call prop_bend_error(C,v3,v4,ve,gam_slope)
                    loc = 4;
                else                
                    slip = rem_slip !Stays in domain 3
                end if
                !Slip in domain 3 is parallel to the domain 3-4 boundary, so the point cannot enter domain 4 if slip (and thus propagation) is forward.        
            else !Backward motion
                slip2 = calc_slip(pt,pte,v0,v3,ve,gam1_slope) !Slip needed to reach domain 2.
                slip4 = calc_slip(pt,pte,v0,v3,ve,gam_slope) !Slip needed to reach domain 4.
                if ((.not. loc_passed(2)) .and. rem_slip < slip2 .and. (slip4<slip2 .or. slip4>0 .or. loc_passed(4)) &
                        .and. slip2<0) then !Goes into domain 2. I'm pretty sure this can only happen if P/S <1, otherwise material flows the other way.
                    slip = slip2
                    call prop_bend_error(C,v3,v2,ve,gam1_slope)
                    loc = 2
                else if ((.not. loc_passed(4)) .and. rem_slip < slip4 .and. slip4<0) then !Goes into domain 4
                    slip = slip4
                    call prop_bend_error(C,v3,v4,ve,gam_slope)
                    loc = 4
                else
                    slip = rem_slip
                end if
            end if
            pt = pt+(slip/v0)*v3
        else if (loc == 4) then !Domain 4
            loc_passed(4) = .true.
            if (slip_sign == 1) then !Forward motion
                slip3 = calc_slip(pt,pte,v0,v4,ve,gam_slope) !Slip needed to reach domain 3.
                if ((.not. loc_passed(3)) .and. rem_slip>slip3 .and. pt(2)+(slip3/v0)*v4(2) > pte(2)+(slip3/v0)*ve(2) &
                    .and. slip3 > 0) then
                    slip = slip3
                    call prop_bend_error(C,v4,v3,ve,gam_slope)
                    loc = 3
                else
                    slip2 = calc_slip(pt,tip,v0,v4,vt,kap_slope) !Slip needed to reach domain 2.
                    if ((.not. loc_passed(2)) .and. rem_slip > slip2 .and. pt(2)+(slip2/v0)*v4(2) < pte(2)+(slip2/v0)*ve(2) &
                        .and. slip2 > 0) then
                        slip = slip2
                        call prop_bend_error(C,v4,v2,vt,kap_slope)
                        loc = 2
                    else
                        slip = rem_slip
                    end if                                
                end if
            else !Backward motion
                slip2 = calc_slip(pt,tip,v0,v4,vt,kap_slope) !Slip needed to reach domain 2.
                slip3 = calc_slip(pt,pte,v0,v4,ve,gam_slope); !Slip needed to reach domain 3.
                slip5 = calc_slip(pt,tip,v0,v4,vt,gam_slope) !Slip needed to reach domain 5 (footwall).
                if ((.not. loc_passed(2)) .and. rem_slip < slip2 .and. (slip5<slip2 .or. slip5>0  .or. loc_passed(5)) &
                    .and. pt(2)+(slip2/v0)*v4(2) < pte(2)+(slip2/v0)*ve(2) .and. slip2<0) then
                    slip = slip2
                    call prop_bend_error(C,v4,v2,vt,kap_slope)
                    loc = 2
                else if ((.not. loc_passed(3)) .and. rem_slip < slip3 .and. (slip5<slip3 .or. slip5>0  .or. loc_passed(5)) &
                    .and. pt(2)+(slip3/v0)*v4(2) > pte(2)+(slip3/v0)*ve(2) .and. slip3<0) then
                    slip = slip3;
                    call prop_bend_error(C,v4,v3,ve,gam_slope)
                    loc = 3;
                else if ((.not. loc_passed(5)) .and. rem_slip < slip5 .and. slip5<0) then
                    slip = slip5
                    call prop_bend_error(C,v4,v5,vt,gam_slope)
                    loc = 5
                else
                    slip = rem_slip
                end if                        
            end if
            pt = pt+(slip/v0)*v4
        else if (loc == 5) then !Footwall (domain 5)
            loc_passed(5) = .true.
            if (slip_sign == 1) then !Forward motion
                slip4 = calc_slip(pt,tip,v0,v5,vt,gam_slope) !Slip needed to reach domain 4.
                if ((.not. loc_passed(4)) .and. rem_slip > slip4 .and. pt(2)+(slip4/v0)*v5(2) > tip(2)+(slip4/v0)*vt(2)) then
                    slip = slip4
                    call prop_bend_error(C,v5,v4,vt,gam_slope)
                    loc = 4
                else
                    slip = rem_slip
                end if            
            else
                slip = rem_slip
            end if
        else !Need to find out where the point is.
            slip = 0 !So no more slip gets subtracted from rem_slip while we're just deciding.
            if ((.not. loc_passed(5)) .and. (pt(1)<bend(1) .and. pt(2)-bend(2)<(pt(1)-bend(1))*tan(ramp_angle(2))) .or. &
                (pt(1)>bend(1) .and. pt(1)<tip(1) .and. pt(2)-bend(2)<(pt(1)-bend(1))*tan(ramp_angle(1))) .or. &
                    (pt(1)>tip(1) .and. pt(2)-tip(2)<(pt(1)-tip(1))*gam_slope)) then !Footwall
                loc = 5
            else if ((.not. loc_passed(1)) .and. pt(2)-bend(2)<(pt(1)-bend(1))*gam1_slope) then !Domain 1
                loc = 1
            else if ((.not. loc_passed(2)) .and. (pt(2)<pte(2) .and. pt(1)-tip(1)<(pt(2)-tip(2))/kap_slope) .or. &
                (pt(2)>pte(2) .and. pt(1)-pte(1)<(pt(2)-pte(2))/gam1_slope)) then !Domain 2
                loc = 2
            else if ((.not. loc_passed(3)) .and. pt(2)>pte(2) .and. pt(1)-pte(1)<(pt(2)-pte(2))/gam_slope) then !Domain 3
                loc = 3
            else if ((.not. loc_passed(4)) .and. pt(1)-tip(1)<(pt(2)-tip(2))/gam_slope) then !Domain 4
                loc = 4
            else
                print*,'Error: point is not in any of the domains' !This shouldn't happen.
            end if
        end if
        rem_slip = rem_slip - slip !remaining slip
    end do
    pts(:,k) = pt !Store the moved point back in the pts array.
    if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
        Cbed(:,:,k) = C
    end if
end do
end subroutine PFPF_points

subroutine PFPF_dips(dippos,dips,v0,total_slip,PoverS,bendxy,ramp_angle,gamma,gamma1,gamma_star,kappa,beta1,beta2, &
    tip_start,sigma_dip_array)
!Subroutine for moving and rotating dips for a parallel fault propagation fold.
!Important notes: We're assuming that the beds start flat. If not, replace gamma, gamma1, and possibly some other angles with angles taking regional dip into account.
!Also, we're assuming the tip starts at the bend. If not, that will change things too.
!I'm not sure the best way to treat the uncertainty in this. Velocity within a domain isn't a function of position, so uncertainty won't change, but when a point reaches
!a boundary, uncertainty in position means there's an uncertainty in where it reaches the boundary and a correlated uncertainty in the amount of slip in each domain.
use options, only: propagate
implicit none
double precision,dimension(1:,1:), intent(inout) :: dippos
double precision,dimension(1:) :: dips
double precision,dimension(2,1), intent(in) :: bendxy
double precision,intent(in) :: v0,total_slip,PoverS
double precision,dimension(2),intent(in) :: ramp_angle !Note: Ramp angle should be acute.
double precision,intent(in) :: gamma,gamma1,gamma_star,kappa,beta1,beta2
double precision,dimension(2),intent(in) :: tip_start
double precision,dimension(1:),intent(inout) :: sigma_dip_array !Array of sigmas for dips (x in first row, y in second)
double precision,dimension(2) :: bend !The position of the bend as a 1D array.
double precision :: ef_constant ! A constant used in calculating the position and velocity of the point e.
double precision :: slip_sign !Forwards = 1, backwards = -1
double precision,dimension(2) :: pte_start,pte !Point e. pte_start is where it is at the beginning of moving each point. pte can move as slip is restored for each point.
double precision :: R0, R1, R2 !Slip ratios.  All expresed relative to slip on the lowest segment.
double precision,dimension(2) :: v1,v2,v3,v4,v5 !Point velocities in each domain
double precision,dimension(2) :: vb,ve,vt !Boundary velocities (and velocities of the points: bend, point e, and tip that the boundary go through)
double precision :: gam1_slope,gam_slope,kap_slope !Boundary slopes (unchanging).
double precision,dimension(2) :: tip !The location of the fault tip. Moves as the fold is restored.
double precision,dimension(2) :: pt !The point being moved
double precision :: dip !The dip of the point being moved.
integer :: k !A counter of points.
double precision :: sigma !Uncertainty (sigma_x,sigma_y) in position for the point being moved. Note: This isn't currently being used.
integer :: loc !tells which zone a point is in.
double precision :: slip,rem_slip !Slip so far and slip remaining
double precision :: slip1,slip2,slip3,slip4,slip5 !Slip amounts needed to reach different boundaries.
logical,dimension(5) :: loc_passed !Tells whether or not (true or false) a point has passed through each of the five domains.
bend = bendxy(:,1) !Only need this to be a 1D array.
ef_constant = sin(ramp_angle(1))*(sin(gamma_star-beta1)/(sin(gamma_star)*sin(2.*gamma_star-beta1))) !To be used in calculating the distance ef or the velocity of point e.
if (v0 >= 0) then !determine the sign of the slip
    slip_sign = 1
else
    slip_sign = -1
end if
!Calculate the slip ratios. All expresed relative to slip on the lowest segment.
R0 = sin(gamma1+ramp_angle(2))/sin(gamma1+ramp_angle(1));
R1 = R0*sin(gamma1+ramp_angle(1))/sin(gamma1+gamma);
R2=R0*sin(beta2)/sin(beta2-ramp_angle(1)+gamma);
if (slip_sign == 1) then !Forwards slip
    pte_start = tip_start !Note: This would not work for normal faults, if they were allowed.
else !Backwards slip
    pte_start = tip_start-PoverS*R0*total_slip*ef_constant*(/cos(kappa),sin(kappa)/) !minus because total_slip is negative.
end if
!Calculate velocities for each domain
v1 = v0*(/cos(ramp_angle(2)),sin(ramp_angle(2))/)
v2 = R0*v0*(/cos(ramp_angle(1)),sin(ramp_angle(1))/)
v3 = R1*v0*(/cos(gamma),sin(gamma)/)
v4 = R2*v0*(/cos(gamma),sin(gamma)/)
v5 = (/0,0/)
!Calculate velocities for each boundary
vb = (/0,0/) !Velocity of the bend: doesn't move.
ve = R0*v0*PoverS*(/cos(ramp_angle(1))+ef_constant*cos(kappa),sin(ramp_angle(1))+ef_constant*sin(kappa)/) !Velocity of point e.
vt = R0*v0*PoverS*(/cos(ramp_angle(1)),sin(ramp_angle(1))/) !Velocity of the fault tip.
!Calculate slopes for each boundary
gam1_slope = -tan(gamma1) !Slope of axes with dip gamma, dipping to the right.
gam_slope = tan(gamma) !Slope of axes with dip gamma1, dipping to the left.
kap_slope = tan(kappa) !Slope of boundary with dip kappa.
do k = 1,size(dippos,2)
    tip = tip_start !Where the tip starts moving from. (Generally the final tip position for restoring a fault-propagation fold.)
    pt = dippos(:,k) !The point to be moved.
    dip = dips(k) !The dip at the point to be moved.
    pte = pte_start !Location e in Figure 11.6  of Allmendinger et al. (2012) and Fig. 16 of Suppe and Medwedeff (1990)
    sigma = sigma_dip_array(k) !Uncertainty of the dip to be moved.
    slip = 0 !slip so far
    rem_slip = total_slip !Remaining slip
    loc = 0 !location, zone numbers from Figure 11.6, Allmendinger et al. (2012), 5 = footwall, 0 = unknown
    loc_passed = .false.
    do while (rem_slip*slip_sign>0.)
        tip = tip + (slip/v0)*vt
        pte = pte + (slip/v0)*ve
        if (loc == 1) then !Domain 1
            loc_passed(1) = .true.
            if (slip_sign == 1) then !Forward motion
                slip2 = calc_slip(pt,bend,v0,v1,vb,gam1_slope) !Slip needed to reach domain 2
                if ((.not. loc_passed(2)) .and. rem_slip > slip2) then !Goes into domain 2
                    slip = slip2
                    call calc_dip_change(dip,sigma,v1,v2,vb,gam1_slope)
                    loc = 2
                else
                    slip = rem_slip !Stays in domain 1
                end if
            else !Backwards motion
                slip = rem_slip !Stays in domain 1
            end if
            pt = pt+(slip/v0)*v1
        else if (loc == 2) then !Domain 2
            loc_passed(2) = .true.
            if (slip_sign == 1) then !Forward motion
                slip3 = calc_slip(pt,pte,v0,v2,ve,gam1_slope) !Slip needed to reach domain 3
                if ((.not. loc_passed(3)) .and. rem_slip > slip3 .and. pt(2)+(slip3/v0)*v2(2) > pte(2)+(slip3/v0)*ve(2) &
                    .and. slip3>0) then !Goes into domain 3             
                    slip = slip3
                    call calc_dip_change(dip,sigma,v2,v3,ve,gam1_slope)
                    loc = 3
                else
                    slip4 = calc_slip(pt,tip,v0,v2,vt,kap_slope) !Slip needed to reach domain 4
                    if ((.not. loc_passed(4)) .and. rem_slip>slip4 .and. pt(2)+(slip4/v0)*v2(2) < pte(2)+(slip4/v0)*ve(2) &
                        .and. slip4>0)then !Goes into domain 4
                        slip = slip4
                        call calc_dip_change(dip,sigma,v2,v4,vt,kap_slope)
                        loc = 4
                    else
                        slip = rem_slip !Stays in domain 2
                    end if
                end if
            else !Backward motion
                slip1 = calc_slip(pt,bend,v0,v2,vb,gam1_slope)
                slip3 = calc_slip(pt,pte,v0,v2,ve,gam1_slope)
                slip4 = calc_slip(pt,tip,v0,v2,vt,kap_slope)
                if ((.not. loc_passed(4)) .and. rem_slip < slip4 .and. slip4<0 .and. pt(2)+(slip4/v0)*v2(2) &
                    < pte(2)+(slip4/v0)*ve(2) .and. slip1<slip4) then !Goes into domain 4.
                    slip = slip4
                    call calc_dip_change(dip,sigma,v2,v4,vt,kap_slope)
                    loc = 4
                else if ((.not. loc_passed(3)) .and. rem_slip < slip3 .and. slip3<0 .and. pt(2)+(slip3/v0)*v2(2) &
                    > pte(2)+(slip3/v0)*ve(2) .and. slip1<slip3) then !Goes into domain 3.
                    slip = slip3
                    call calc_dip_change(dip,sigma,v2,v3,ve,gam1_slope)
                    loc = 3
                else if ((.not. loc_passed(1)) .and. rem_slip < slip1) then !Goes into domain 1. < because slips are negative
                    slip = slip1
                    call calc_dip_change(dip,sigma,v2,v1,vb,gam1_slope)
                    loc = 1
                else
                    slip = rem_slip !Stays in domain 2
                end if
            end if
            pt = pt+(slip/v0)*v2
        else if (loc == 3) then !Domain 3
            loc_passed(3) = .true.
            if (slip_sign == 1) then !Forward motion
                slip2 = calc_slip(pt,pte,v0,v3,ve,gam1_slope) !Slip needed to reach domain 2
                slip4 = calc_slip(pt,pte,v0,v3,ve,gam_slope) !Slip needed to reach domain 4.
                if ((.not. loc_passed(2)) .and. rem_slip > slip2 .and. slip2 > 0  .and. &
                        (slip4>slip2 .or. slip4<0  .or. loc_passed(4))) then !Goes into domain 2
                    slip = slip2
                    call calc_dip_change(dip,sigma,v3,v2,ve,gam1_slope)
                    loc = 2
                else if ((.not. loc_passed(4)) .and. rem_slip > slip4 .and. slip4>0) then !Goes into domain 4
                    slip = slip4;
                    call calc_dip_change(dip,sigma,v3,v4,ve,gam_slope)
                    loc = 4;
                else                
                    slip = rem_slip !Stays in domain 3
                end if
                !Slip in domain 3 is parallel to the domain 3-4 boundary, so the point cannot enter domain 4 if slip (and thus propagation) is forward.        
            else !Backward motion
                slip2 = calc_slip(pt,pte,v0,v3,ve,gam1_slope) !Slip needed to reach domain 2.
                slip4 = calc_slip(pt,pte,v0,v3,ve,gam_slope) !Slip needed to reach domain 4.
                 if ((.not. loc_passed(2)) .and. rem_slip < slip2 .and. (slip4<slip2 .or. slip4>0  .or. loc_passed(4)) &
                        .and. slip2<0) then !Goes into domain 2. I'm pretty sure this can only happen if P/S <1, otherwise material flows the other way, but I don't think it's quite that simple.
                    slip = slip2
                    call calc_dip_change(dip,sigma,v3,v2,ve,gam1_slope)
                    loc = 2
                else if ((.not. loc_passed(4)) .and. rem_slip < slip4 .and. slip4<0) then !Goes into domain 4
                    slip = slip4
                    call calc_dip_change(dip,sigma,v3,v4,ve,gam_slope)
                    loc = 4
                else
                    slip = rem_slip
                end if
            end if
            pt = pt+(slip/v0)*v3
        else if (loc == 4) then !Domain 4
            loc_passed(4) = .true.
            if (slip_sign == 1) then !Forward motion
                slip3 = calc_slip(pt,pte,v0,v4,ve,gam_slope) !Slip needed to reach domain 3.
                if ((.not. loc_passed(3)) .and. rem_slip>slip3 .and. pt(2)+(slip3/v0)*v4(2) > pte(2)+(slip3/v0)*ve(2) &
                    .and. slip3 > 0) then
                    slip = slip3
                    call calc_dip_change(dip,sigma,v4,v3,ve,gam_slope)
                    loc = 3
                else
                    slip2 = calc_slip(pt,tip,v0,v4,vt,kap_slope) !Slip needed to reach domain 2.
                    if ((.not. loc_passed(2)) .and. rem_slip > slip2 .and. pt(2)+(slip2/v0)*v4(2) < pte(2)+(slip2/v0)*ve(2) &
                        .and. slip2 > 0) then
                        slip = slip2
                        call calc_dip_change(dip,sigma,v4,v2,vt,kap_slope)
                        loc = 2
                    else
                        slip = rem_slip
                    end if                                
                end if
            else !Backward motion
                slip2 = calc_slip(pt,tip,v0,v4,vt,kap_slope) !Slip needed to reach domain 2.
                slip3 = calc_slip(pt,pte,v0,v4,ve,gam_slope); !Slip needed to reach domain 3.
                slip5 = calc_slip(pt,tip,v0,v4,vt,gam_slope) !Slip needed to reach domain 5 (footwall).
                if ((.not. loc_passed(2)) .and. rem_slip < slip2 .and. (slip5<slip2 .or. slip5>0 .or. loc_passed(5)) &
                    .and. pt(2)+(slip2/v0)*v4(2) < pte(2)+(slip2/v0)*ve(2) .and. slip2<0) then
                    slip = slip2
                    call calc_dip_change(dip,sigma,v4,v2,vt,kap_slope)
                    loc = 2
                else if ((.not. loc_passed(3)) .and. rem_slip < slip3 .and. (slip5<slip3 .or. slip5>0 .or. loc_passed(5)) &
                    .and. pt(2)+(slip3/v0)*v4(2) > pte(2)+(slip3/v0)*ve(2) .and. slip3<0) then
                    slip = slip3;
                    call calc_dip_change(dip,sigma,v4,v3,ve,gam_slope);
                    loc = 3;
                else if ((.not. loc_passed(5)) .and. rem_slip < slip5 .and. slip5<0) then
                    slip = slip5
                    call calc_dip_change(dip,sigma,v4,v5,vt,gam_slope)
                    loc = 5
                else
                    slip = rem_slip
                end if                        
            end if
            pt = pt+(slip/v0)*v4
        else if (loc == 5) then !Footwall (domain 5)
            loc_passed(5) = .true.
            if (slip_sign == 1) then !Forward motion
                slip4 = calc_slip(pt,tip,v0,v5,vt,gam_slope) !Slip needed to reach domain 4.
                if ((.not. loc_passed(4)) .and. rem_slip > slip4 .and. pt(2)+(slip4/v0)*v5(2) > tip(2)+(slip4/v0)*vt(2)) then
                    slip = slip4
                    call calc_dip_change(dip,sigma,v5,v4,vt,gam_slope)
                    loc = 4
                else
                    slip = rem_slip
                end if            
            else
                slip = rem_slip
            end if
        else !Need to find out where the point is.
            slip = 0 !So no more slip gets subtracted from rem_slip while we're just deciding.
            if ((.not. loc_passed(5)) .and. (pt(1)<bend(1) .and. pt(2)-bend(2)<(pt(1)-bend(1))*tan(ramp_angle(2))) .or. &
                (pt(1)>bend(1) .and. pt(1)<tip(1) .and. pt(2)-bend(2)<(pt(1)-bend(1))*tan(ramp_angle(1))) .or. &
                    (pt(1)>tip(1) .and. pt(2)-tip(2)<(pt(1)-tip(1))*gam_slope)) then !Footwall
                loc = 5
            else if ((.not. loc_passed(1)) .and. pt(2)-bend(2)<(pt(1)-bend(1))*gam1_slope) then !Domain 1
                loc = 1
            else if ((.not. loc_passed(2)) .and. (pt(2)<pte(2) .and. pt(1)-tip(1)<(pt(2)-tip(2))/kap_slope) .or. &
                (pt(2)>pte(2) .and. pt(1)-pte(1)<(pt(2)-pte(2))/gam1_slope)) then !Domain 2
                loc = 2
            else if ((.not. loc_passed(3)) .and. pt(2)>pte(2) .and. pt(1)-pte(1)<(pt(2)-pte(2))/gam_slope) then !Domain 3
                loc = 3
            else if ((.not. loc_passed(4)) .and. pt(1)-tip(1)<(pt(2)-tip(2))/gam_slope) then !Domain 4
                loc = 4
            else
                print*,'Error: point is not in any of the domains' !This shouldn't happen.
            end if
        end if
        rem_slip = rem_slip - slip !remaining slip
    end do
    dippos(:,k) = pt !Store the moved point back in the pts array.
    dips(k) = dip !Store the changed dip.
    if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
        sigma_dip_array(k) = sigma !Store the changed uncertainties.
    end if
end do
end subroutine PFPF_dips

function calc_slip(pt,ptb,v0,v1,vb,slope)
!This function calculates the slip necessary for a point moving in a constant velocity domain to reach a boundary.
!Everything with dimension(2) is x coordinate or x component followed by y coordinate or y component.
double precision,dimension(2),intent(in) :: pt !The point (x,y) that has to move to the boundary.
double precision,dimension(2),intent(in) :: ptb !A point that the boundary passes through and moves with.
double precision,intent(in) :: v0 !Velocity relative to which slip is calculated. S = v0*time
double precision,dimension(2),intent(in) :: v1 !Velocity at which pt is moving
double precision,dimension(2),intent(in) :: vb !Velocity at which ptb (and therefore the boundary) is moving
double precision,intent(in) :: slope !The slope of the boundary. This is not allowed to change. Equals tan(boundary dip).
double precision :: calc_slip !The slip necessary for pt to reach the boundary
calc_slip = ((pt(1)-ptb(1))*slope-pt(2)+ptb(2))/((1./v0)*(v1(2)-vb(2)-(v1(1)-vb(1))*slope))
end function calc_slip

subroutine calc_dip_change(dip,sigma_dip,v1,v2,vb,slope)
!This function calculates the changed dip and uncertainty in dip for a point crossing a boundary between two constant velocity domains.
!Everything with dimension(2) is x coordinate or x component followed by y coordinate or y component.
use math, only: sec
implicit none
double precision,intent(inout) :: dip !The dip at the point.
double precision,intent(inout) :: sigma_dip !Uncertainty in the dip.
double precision,dimension(2),intent(in) :: v1 !Velocity in domain point is coming from.
double precision,dimension(2),intent(in) :: v2 !Velocity in domain point is moving into.
double precision,dimension(2),intent(in) :: vb !Velocity at which the domain boundary is moving.
double precision,intent(in) :: slope !The slope of the boundary. This is not allowed to change. Equals tan(boundary dip).
double precision :: A,B,C,W,X,Y,Z !Terms to be used in the calculations.
double precision :: tan_dip,secsq_dip !Tangent and secant^2 of the (original) dip
dip = -dip !An unfortunate consequence of the sign convention I have been using, with dips positive down to the right.
tan_dip = tan(dip)
secsq_dip = (sec(dip))**2
A = v1(1)*tan_dip-v1(2)
B = vb(1)*slope-vb(2)
C = slope-tan_dip
W = A*slope-B*tan_dip+C*v2(2)
X = A-B+C*v2(1)
Y = v1(1)*secsq_dip*slope-B*secsq_dip-v2(2)*secsq_dip
Z = v1(1)*secsq_dip-v2(1)*secsq_dip
dip = atan(W/X) !New dip
sigma_dip = abs((Y*X-W*Z)/(W**2+X**2))*sigma_dip !New uncertainty in dip.
dip = -dip !Switch back to the sign convention used in the rest of the program.
end subroutine calc_dip_change

end module parallel_FPF
