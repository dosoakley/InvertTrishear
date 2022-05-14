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

module trishear_listric
!Trishear for a circular, listric fault, using fault parallel flow.

!Note: This doesn't seem to work properly.
!Also it doesn't do full error propagation; it only does it for the trishear zone.

integer :: decoly_at_tipy !Tells whether decollement depth and tipy must be the same (1) or not (0)

contains

subroutine trishear_listric_options
use parameters, only: nparams
use options, only: options_id,TipToSolve
implicit none
do
    if (TipToSolve==1) then !If solving for the initial tip
        print*,'Must detachment depth be at initial fault tip depth? (0/1)'
        read(options_id,*) decoly_at_tipy
    else
        decoly_at_tipy = 0
    end if    
    if (decoly_at_tipy == 0) then
        nparams = 9
        exit
    else if (decoly_at_tipy == 1) then
        nparams = 8
        exit
    else
        print*,'Invalid Command'
    end if
end do
end subroutine trishear_listric_options

subroutine run_model_listric(output,params)
use constants
use parameters, only: nparams,increment,v0,slip_sense,nparam_start_restored_fit,&
    nrestored_fit_params,nparam_start_growth,nparam_start_terr
use options, only: ResultType,DatTypes,TipToSolve,ndatatypes_poss,DatType_beds,DatType_dips,DatType_faultpts,DatType_terr,FitType,&
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
double precision :: tipx,tipy,total_slip,total_slip_mag,theta_max,theta_max_deg,phi,phi_deg,PoverS,s,m,yd,Rc !parameters
double precision :: slip !Slip to restore a bed, may be growth or pregrowth.
double precision,dimension(nterr) :: terr_slip !Slip needed to restore each terrace
double precision,dimension(ngrowth) :: growth_slip !Slip needed to restore each growth strata bed.
double precision,dimension(nrestored_fit_params) :: restored_fit_params
double precision,dimension(2) :: cc !Position (x,y) of the center of the circle
double precision :: thetat,theta0 !Angle to the tip up from vertical. [0,180], and lowest thetat
double precision :: theta_max_acute !Theta_max as an acute angle
integer :: ramp_dir !dip direction of ramp
integer :: i,n,g ! Counters. n and g are for all beds and growth beds respectively.
logical :: goodrun !tells if the run is good
double precision :: NaN = 0 !For creating NaN by dividing by zero
integer :: tip_loc,tip_loc_init,tip_loc_final !Tells which location the tip is in: 1 = curve, 2 = straight; flat is not allowed
double precision,dimension(2) :: list_str !point (x,y) where the fault change from listric to straight
double precision,dimension(:,:),allocatable :: thisbed !the bed being worked on
double precision :: slope,b !Slope and intercept of a bed
double precision :: slope_last, b_last !Slope and intercept of the previous bed restored.
double precision,dimension(2) :: bed_xlims,bed_xlims_last !x limits of current and previous bed restored.
double precision :: dist !The distance between (tipx,tipy) and the bend
double precision :: tipxinit, tipyinit, tipxfinal, tipyfinal !Initial and final (in terms of forward motion) tip positions
integer :: bedlength !number of points in the bed being worked on
double precision,dimension(:,:),allocatable :: dip_pos_temp !dip positions array that I can change.
double precision,dimension(ndips) :: dips_temp !dips array that I can change
double precision,dimension(ndips) :: RestoredDips !values of dips after inverse fault movement is complete
double precision,dimension(nbeds) :: beds_output_indiv !Output from each of the beds individually.
double precision,dimension(n_restored_beds) :: restored_beds_output_indiv !Output from each of the beds individually.
double precision,dimension(nterr) :: terr_output_indiv !Output from each of the terraces individually.
double precision,dimension(nfaultpts) :: fault_errs !The individual fault point errors
double precision,dimension(2,nfaultpts) :: faultpts_temp !Array of fault points that can be changed.
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

NaN = NaN/NaN
goodrun = .true. !Initialize it to true. So far, this is a good run.
do !This do loop lets us exit anytime it turns out we have an impossible combination of parameters
!Make assignments
tipx = params(1)
tipy = params(2)
total_slip_mag = params(3)
theta_max_deg = params(4)
phi_deg = params(5)
PoverS = params(6)
s = params(7)
Rc = params(8)
if (decoly_at_tipy == 0) then
    yd = params(9)
    if (tipy < yd) then !Whether it's tipy initial or final, it can't be below the detachment.
        goodrun = .false.
        exit
    end if
else
    yd = tipy !Only possible if tipy = tipy_init. Currently not allowing the other way.
end if
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
theta_max = theta_max_deg*pi/180; !convert to radians
theta_max_acute = min(theta_max,pi-theta_max) !Acute version of ramp angle
if (theta_max<=pi/2) then !1=dips left, -1 = dips right
    ramp_dir = 1
else
    ramp_dir = -1
end if
cc(2) = yd+Rc
if (tipy<cc(2)) then !Otherwise, the tip can't possibly be on the curve.
    cc(1) = tipx-ramp_dir*sqrt(Rc**2-(tipy-cc(2))**2) !ccx if the tip is in the curve
    thetat = atan(ramp_dir*(tipx-cc(1))/(cc(2)-tipy)) !Corrected to be the same whether the ramp dips left or right
else
    thetat = theta_max_acute !Should this be theta_max_acute?
end if
if (thetat >= theta_max_acute .or. thetat<0) then   !Then tip is in the straight part
    cc(1) = tipx-ramp_dir*Rc*sin(theta_max)+(cc(2)-ramp_dir*Rc*cos(theta_max)-tipy)/tan(theta_max)
    tip_loc = 2
else
    tip_loc = 1
end if
phi = phi_deg*pi/180 !convert to radians
if (phi >= (pi-theta_max_acute)/2) then !Check that the trishear boundary won't cross the kink axis
    goodrun = .false.
    exit
end if
m = tan(phi)
if (phi == pi/2) then
        m = tan(89*pi/180) !Since tan(pi) = inf, which causes problems
end if
!Calculate the position of the change from listric to straight
list_str(1) = cc(1)+ramp_dir*Rc*sin(theta_max)
list_str(2) = cc(2)-ramp_dir*Rc*cos(theta_max)
!list_str(2) = yd+Rc*(1.-cos(theta_max)) !Problem with this version was that it only works when theta_max is acute.
!Calculate initial and final tipx and tipy (final = present-day, initial = when trishear motion began)
if (TipToSolve == 1) then !If tipx, tipy are initial
    !Check that tipxinit, tipyinit don't conflict with detachment depth and circle position
    if (tipy < yd .or. ramp_dir*(tipx-cc(1))<0) then
        goodrun = .false.
        exit
    end if
    tipxinit = tipx
    tipyinit = tipy
    !Note: If we allowed tipxfinal to be lower than tipxinit, this wouldn't be so simple. Same for if TipToSolve = 2 below
    if (tip_loc==2) then !Starts and stays in straight segment
        tipxfinal = tipx+total_slip_mag*PoverS*cos(theta_max) !Remember that if doing an inversion, total_slip and slip_sense will be negative
        tipyfinal = tipy+total_slip_mag*PoverS*sin(theta_max)
        tip_loc_init = 2
        tip_loc_final = 2
    else !Starts in the listric segment
        dist = (theta_max_acute-thetat)*Rc !Should always be positive.
        if (total_slip_mag*PoverS <= dist) then !Starts and stays in listric segment
            tipxfinal = cc(1)+ramp_dir*Rc*sin(thetat+total_slip_mag*PoverS/Rc)
            tipyfinal = cc(2)-Rc*cos(thetat+total_slip_mag*PoverS/Rc)
            tip_loc_init = 1
            tip_loc_final = 1
        else !Goes from listric into straight segment
            tipxfinal = list_str(1) + (total_slip_mag*PoverS-dist)*cos(theta_max)
            tipyfinal = list_str(2) + (total_slip_mag*PoverS-dist)*sin(theta_max)
            tip_loc_init = 1
            tip_loc_final = 2
        end if
    end if
else !If tipx, tipy are final
    tipxfinal = tipx
    tipyfinal = tipy
    if (tip_loc == 1) then !Starts and stays in listric segment
        if (total_slip_mag*PoverS/Rc>thetat) then !Too much slip for the given fault geometry
            goodrun = .false.
            exit
        end if
        tipxinit = cc(1)+ramp_dir*Rc*sin(thetat-total_slip_mag*PoverS/Rc)
        tipyinit = cc(2)-Rc*cos(thetat-total_slip_mag*PoverS/Rc)
        tip_loc_init = 1
        tip_loc_final = 1
    else !tipfinal in the straight segment
        dist = sqrt((tipy-list_str(2))**2+(tipx-list_str(1))**2)
        if ((total_slip_mag*PoverS-dist)/Rc>theta_max_acute) then !Too much slip for the given fault geometry
            goodrun = .false.
            exit
        end if
        if (total_slip_mag*PoverS <= dist) then !Starts and stays in straight segment
            tipxinit = tipx-total_slip_mag*PoverS*cos(theta_max)
            tipyinit = tipy-total_slip_mag*PoverS*sin(theta_max)
            tip_loc_init = 2
            tip_loc_final = 2
        else !tip_init in listric segment, tip_final in straight segment
            tipxinit = cc(1)+ramp_dir*Rc*sin(theta_max_acute-(total_slip_mag*PoverS-dist)/Rc)
            tipyinit = cc(2) - Rc*cos(theta_max_acute-(total_slip_mag*PoverS-dist)/Rc)
            tip_loc_init = 1
            tip_loc_final = 2
        end if
    end if
end if
!Calculate the initial and final theta (tip position). These should always be acute angles. Note: This thetat is only right for inversion, not for forward motion.
if (tip_loc_final == 1) then !In listric segment
    thetat = atan(ramp_dir*(tipxfinal-cc(1))/(cc(2)-tipyfinal)) !Corrected to be the same whether the ramp dips left or right
else !In straight segment
    thetat = theta_max_acute
end if
if (tip_loc_init == 1) then !In listric segment
    theta0 = atan(ramp_dir*(tipxinit-cc(1))/(cc(2)-tipyinit)) !Corrected to be the same whether the ramp dips left or right
    if (tip_loc_final == 2) then
        if (tip_loc_init == 1 .and. abs(theta_max_acute-(total_slip_mag*PoverS-dist)/Rc-theta0)>1e-5) then
            print*,'Error: theta0 not right'
            print*,theta0,theta_max_acute-(total_slip_mag*PoverS-dist)/Rc
        end if
    end if
    if (theta0 < 0) then
        print*,'Error: theta0 < 0'
        print*,tipxinit-cc(1)
        print*,cc(2)-tipyinit
        print*,(total_slip_mag*PoverS-dist)/Rc !This is coming up much too large. Possibly these are models that should be rejected for having too high slip given the detachment depth.
    end if
else !In straight segment
    theta0 = theta_max_acute
end if
!Initialize the outputs array.
outputs = 0 !For RMS or chi squared it should be 0 anyway. For probabilities, p should be 1, so ln(p) = 0.
!Now do inversion of bed and dip data.
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
        thisbed(:,:) = beds(n)%pts
        if (ramp_dir == -1) then !Reflect everything around the y axis, so that the fault now dips left.
            thisbed(1,:) = -thisbed(1,:)
            tipxfinal = -tipxfinal !Note: I could just put ramp_dir*tipxfinal in the call trishear_func_listric line below.
            cc(1) = -cc(1)
            Cbed(1,2,:)=-Cbed(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cbed(2,1,:)=-Cbed(2,1,:)
        end if
        if (beds(n)%growth .eqv. .true.) then !For growth strata
            slip = growth_slip(g)
            g = g+1
        else
            slip = total_slip
        end if
        call trishear_func_listric(thisbed,v0,phi,m,slip,increment,-slip_sense*PoverS,s,theta_max_acute,theta0,thetat, &
            (/tipxfinal,tipyfinal/),cc,Rc,Cbed) !Run the trishear on this bed
        if (ramp_dir == -1) then
            thisbed(1,:) = -thisbed(1,:) !Reflect back
            tipxfinal = -tipxfinal
            cc(1) = -cc(1)
            Cbed(1,2,:)=-Cbed(1,2,:) !This could also be done by multiplying by an array of the form [-1,0;0,1].
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
if (goodrun .eqv. .false.) then !If the above numerical solutions didn't work out, then discard the run.
    exit !Leave the loop now before wasting any more time or risking issues with later parts of the run.
end if
if (DatTypes(DatType_dips) == 1) then !Dips
    call GetSigmaDip(params,sigma_dip,sigma_restored_dip)
    dip_pos_temp = dip_pos !Temporary dip positions for this inversion, so we don't mess with the original ones
    dips_temp = ramp_dir*dips !Rotate if necessary
        if (ramp_dir == -1) then !Reflect everything around the y axis, so that the fault now dips left.
            dip_pos_temp(1,:) = -dip_pos_temp(1,:)
            tipxfinal = -tipxfinal !Note: I could just put ramp_dir*tipxfinal in the call trishear_func_listric line below.
            cc(1) = -cc(1)
        end if
    sigma_dip_array = sigma_dip*pi/180.
    call trishear_func_dip_listric(dip_pos_temp,dips_temp,v0,phi,m,total_slip,increment,-slip_sense*PoverS,s,theta_max_acute,&
        theta0,thetat,(/tipxfinal,tipyfinal/),cc,Rc,sigma_dip_array)
    if (ramp_dir == -1) then !Reflect back. We don't need to bother with dip_pos_temp, but tipxfinal or cc(1) may be needed still.
        tipxfinal = -tipxfinal
        cc(1) = -cc(1)
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
        thisbed(:,:) = terraces(n)%pts
        Cterr = 0 !Initialize covariance matrix for bed points (x,y).
        Cterr(1,1,:) = sigma_terr(n)**2
        Cterr(2,2,:) = sigma_terr(n)**2
        if (ramp_dir == -1) then !Reflect everything around the y axis, so that the fault now dips left.
            thisbed(1,:) = -thisbed(1,:)
            tipxfinal = -tipxfinal !Note: I could just put ramp_dir*tipxfinal in the call trishear_func_listric line below.
            cc(1) = -cc(1)
            Cbed(1,2,:)=-Cbed(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cbed(2,1,:)=-Cbed(2,1,:)
        end if
        call trishear_func_listric(thisbed,v0,phi,m,terr_slip(n),increment,-slip_sense*PoverS,s,theta_max_acute,theta0,thetat, &
            (/tipxfinal,tipyfinal/),cc,Rc,Cterr) !Run the trishear on this terrace
        if (ramp_dir == -1) then
            thisbed(1,:) = -thisbed(1,:) !Reflect back
            tipxfinal = -tipxfinal
            cc(1) = -cc(1)
            Cbed(1,2,:)=-Cbed(1,2,:) !This could also be done by multiplying by a reflection array, as in the others.
            Cbed(2,1,:)=-Cbed(2,1,:)
        end if
        call calc_terr_err(n,bedlength,thisbed,terr_output_indiv,Cterr,sigma_orig_elev)
        deallocate(thisbed,Cterr)
    end do
    call calc_beds_err(terr_output_indiv,outputs(3))
end if
if (DatTypes(DatType_faultpts) == 1) then !Fault points
    call GetSigmaFault(params,sigma_fault)
    faultpts_temp = faultpts
    if (ramp_dir == -1) then
        faultpts_temp = -faultpts_temp !Reflect
        cc(1) = -cc(1)
    end if
    call fault_pts_err_listric(faultpts_temp,thetat,cc,Rc,fault_errs) !Note: This is only right for inversion where thetat = thetat final. It won't work properly for forward motion.
    if (ramp_dir == -1) then !Reflect back
        cc(1) = -cc(1) !May not be needed at present, since we don't use this again.
    end if
    do n = 1,nfaultpts
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

end subroutine run_model_listric

subroutine trishear_func_listric(pts,v0,phi,m,total_slip,increment,PoverS,s,theta_max,theta0,thetat_start,tip_pos,cc,Rc, &
    Cbed)
!Trishear for a listric fault-propagation fold rising from a horizontal detachment.
!This version assumes fault-parallel slip outside the trishear zone, thus conserving velocity.
!Points are initially in Cartesian world coordinates, not trishear coordinates.
!Fault is right-vergent in the cross section.
use constants
use math, only: xy_to_ze2D,ze_to_xy2D
use options, only: propagate
implicit none
double precision, dimension(1:,1:), intent(inout) :: pts
double precision, intent(in) :: v0,phi,m,total_slip,increment,PoverS,s
double precision,intent(in),dimension(2) :: tip_pos !Initial position (x,y) of the fault tip
double precision,intent(in),dimension(2) :: cc !Position (x,y) of the center of the circle
double precision,intent(in) :: theta_max !Final (maximum) thetat
double precision,intent(inout) :: thetat_start !Angle between radius to the fault tip and downwards vertical. Starting value.
double precision,intent(in) :: theta0 !Angle between radius to initial fault tip (i.e. lowest position it reaches) and downwards vertical. Note: I think what we read in is only correct for backwards slip.
double precision,intent(in) :: Rc !Radius of the circle
double precision,dimension(1:,1:,1:),intent(inout) :: Cbed !Covariance matrices for x and y for all the points in the bed.
double precision,dimension(2,2) :: C,J !Covariance matrix for the specific point being worked on and Jacobian matrix for motion in the trishear velocity field.
double precision,dimension(2,2) :: R !A rotation matrix.
double precision :: thetat !Angle between radius to the fault tip and downwards vertical.
integer :: tip_loc,tip_loc_start !Tells which location the tip is in: 1 = curve, 2 = straight; flat is not allowed
integer :: pt_loc !Tells where the point is: 1 = flat, 2 = curve, 3 = straight, 4 = trishear, 5 = footwall
integer :: slip_sign, prop_sign !slip_sign: 1 = forwards motion, -1 = backwards motion
double precision :: half_v0,xexp,yexp,signy
double precision :: total_slip_abs !absolute value of the total slip
double precision :: slip !the slip produced by a specific step
double precision :: Sremstr !remaining slip while tip is in straight zone
double precision :: Sremcurve !remaining slip while tip is in straight zone
double precision :: total_slip_curve,total_slip_str !Total slip with the tip in each part.
double precision :: xt,yt !Cartesian coordinates of the tip
double precision :: xt0,yt0 !Initial cartesian coordinates of the tip (as in lowest; initial in forwards slip)
double precision :: xchange !x tip coordinate at which it changes from listric to planar
double precision :: vx, vy !x and y components of the velocity
double precision :: x,y !Cartesian coordinates of the point in question
double precision,dimension(2,1) :: pt,ptze !(x,y) in one array, and its trishear equivalent (z,e)
double precision :: ccx,ccy !Coordinates of the center of the circle that the fault forms an arc of
integer :: k !A counter of beds.
double precision :: theta_prop,sin_theta_prop,cos_theta_prop !Theta for one increment of propagation, and its sine and cosine.

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
ccx = cc(1)
ccy = cc(2)
xt0 = ccx+Rc*sin(theta0) !Lowest fault tip coordinates. Note: This is only true if tip does enter curved segment, but if not these values don't get used anyway.
yt0 = ccy-Rc*cos(theta0)
xchange = ccx+Rc*sin(theta_max) !x tip coordinate at which it changes from listric to planar.
theta_prop = PoverS*increment/Rc !theta_prop will be in radians; propagation in one increment of slip
sin_theta_prop = sin(theta_prop) !Calculate these once for efficiency
cos_theta_prop = cos(theta_prop)
if (thetat_start >= theta_max) then !Find out which part of the fault the tip starts in
    tip_loc_start = 2 !In the straight part
else
    tip_loc_start = 1 !In the curved part
end if
if (tip_loc_start == 1) then  !Find out how much slip occurs with the tip in each part
    if (slip_sign == prop_sign) then
        total_slip_curve = min((theta_max-theta0)*Rc/PoverS,total_slip)
        total_slip_str = total_slip - total_slip_curve
    else
        total_slip_curve = total_slip
        total_slip_str = 0
    end if
else
    if (slip_sign == prop_sign) then
        total_slip_curve = 0
        total_slip_str = total_slip
    else
        total_slip_str = max(-sqrt((tip_pos(1)-(ccx+Rc*cos(theta_max-pi/2)))**2 + &
            (tip_pos(2)-(ccy+Rc*sin(theta_max-pi/2)))**2)/PoverS,total_slip) !Max b/c the slips should be negative in this case
        total_slip_curve = total_slip - total_slip_str
        if (atan((ccy-Rc*cos(theta_max)-tip_pos(2))/(ccx+Rc*sin(theta_max)-tip_pos(1)))-theta_max >= 1.0e-5) then
            print*,'Error'
            print*,atan((ccy-Rc*cos(theta_max)-tip_pos(2))/(ccx+Rc*sin(theta_max)-tip_pos(1)))
            print*,atan((ccy-Rc*cos(theta_max)-tip_pos(2))/(ccx+Rc*sin(theta_max)-tip_pos(1)))*180/pi
            print*,theta_max
            print*,theta_max*180./pi
        end if
        if (abs(theta_max+total_slip_curve*PoverS/Rc-theta0)>1e-5) then !Somehow this keeps happening, which doesn't make sense.
            print*,'Error: Tip does not end at theta0'
            print*,theta0,theta_max+total_slip_curve*PoverS/Rc
            print*,ccx+Rc*cos(theta_max-pi/2),ccy+Rc*sin(theta_max-pi/2)
        end if
    end if
end if
!print*,total_slip
do k = 1,size(pts,2)
    !print*,'k = ',k
    x = pts(1,k)
    y = pts(2,k)
    xt = tip_pos(1)
    yt = tip_pos(2)
    Sremcurve = total_slip_curve
    Sremstr = total_slip_str
    slip = 0
    pt_loc = 0 !To start
    C = Cbed(:,:,k) !Covariance matrix for x and y for this point.
    tip_loc = tip_loc_start
    thetat = thetat_start !Really, thetat_start will be theta_max if slip_sign == -1 and theta0 if slip_sign == 1
    do while (slip_sign*Sremcurve>0. .or. slip_sign*Sremstr>0.) !The decrease in precision is to deal with problems with rounding
        slip = 0 !Zero it out for next time. (Especially important if pt_loc == 0)
        !print*,pt_loc
        if (pt_loc == 4 .and. x<(y-yt)/tan(thetat-pi/2)+xt) then
            print*,'Error: Point should not be going into the trishear zone.'
        end if
        select case(pt_loc)
        case (5) !Point is in the footwall
            call footwall_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,xt0,yt0,theta0,&
                cc,Rc,phi,PoverS,theta_max,prop_sign)
        case(3) !Point is in the straight part
            call straight_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,xt0,yt0,theta0,&
                cc,Rc,phi,m,PoverS,theta_max,prop_sign)
        case(2)!Point is in listric part
            call listric_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,cc,Rc,phi,PoverS,theta_max,&
                prop_sign,C)
        case(1)!Point is in the flat part.
            call flat_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,cc,Rc,PoverS,theta_max,&
                prop_sign)
        case(4) !Point is in the trishear zone
            call trishear_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,theta0,cc,Rc,phi,m,PoverS,&
                theta_max,v0,s,theta_prop,sin_theta_prop,cos_theta_prop,C)
        case default
            !Need to figure out where we are
            if ((y <= ccy-Rc .and. x<=ccx) .or. &
                    (x>=ccx .and. y<=(x-ccx)*tan(thetat-pi/2)+ccy .and. (x-ccx)**2+(y-ccy)**2 >= Rc**2) &
                    .or. (x>=xchange .and. x<=xt .and. y<=(x-xt)*tan(theta_max)+yt) .or. &
                    (x>xt .and. y<=(x-xt)*tan(thetat-phi)+yt)) then
                pt_loc = 5 !Footwall
            else if ((x<=ccx .and. y<=ccy) .or. (y>=ccy .and. y<=-(x-ccx)*tan((pi-thetat)/2)+ccy)) then
                pt_loc = 1 !Flat
            else if (y<=ccy .and. atan((x-ccx)/(ccy-y)) <= thetat) then
                pt_loc = 2 !Curve
            else if (x<=(y-yt)/tan(thetat+phi)+xt) then
                pt_loc = 3 !Straight
            else
                pt_loc = 4 !Trishear
                if (pt_loc == 4 .and. ((x<(y-yt)/tan(thetat+phi)+xt) .or. (y<(x-xt)*tan(thetat-phi)+yt))) then
                    print*,'Error: Should not be in trishear zone.'
                    print*,'From unknown'
                end if
                if (x < (y-yt)/tan(thetat-pi/2)+xt) then !Compare to a line perpendicular to the fault. Corresponds to x<0 in trishear coordinates.
                    print*,'Error: Point is at x<0 in trishear coordinate system.'
                end if
            end if
        end select
        if (slip_sign == -1 .and. slip > 0) then !This shouldn't be the case.
            print*,'Error: slip > 0 on inversion.'
            print*,'slip = ',slip
        end if
        if (tip_loc == 1) then !Update the remaining amount of slip with the tip in the appropriate segment
            Sremcurve = Sremcurve - slip
            if (slip_sign*Sremcurve <= 0) then
                tip_loc = 2
            end if        
        else
            Sremstr = Sremstr - slip
            if (slip_sign*Sremstr <=0) then
                tip_loc = 1
            end if
        end if
    end do
    pts(1,k) = x !Store the new x and y position values
    pts(2,k) = y
    if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
        Cbed(:,:,k) = C !Store the new covariance matrix.
    end if
end do
end subroutine trishear_func_listric

subroutine trishear_func_dip_listric(pts,dips,v0,phi,m,total_slip,increment,PoverS,s,theta_max,theta0,thetat_start,tip_pos,cc,Rc,&
    sigma_dip_array)
!Trishear for a listric fault-propagation fold rising from a horizontal detachment.
!This version assumes fault-parallel slip outside the trishear zone, thus conserving velocity.
!Points are initially in Cartesian world coordinates, not trishear coordinates.
!Fault is right-vergent in the cross section.
use constants
use math, only: xy_to_ze2D,ze_to_xy2D
use data_module, only: ndips
use options, only: propagate
implicit none
double precision, dimension(1:,1:), intent(inout) :: pts !The positions of the dips.
double precision,dimension(1:),intent(inout) :: dips !The dips
double precision, intent(in) :: v0,phi,m,total_slip,increment,PoverS,s
double precision,intent(in),dimension(2) :: tip_pos !Initial position (x,y) of the fault tip
double precision,intent(in),dimension(2) :: cc !Position (x,y) of the center of the circle
double precision,intent(in) :: theta_max !Final (maximum) thetat
double precision,intent(inout) :: thetat_start !Angle between radius to the fault tip and downwards vertical. Starting value.
double precision,intent(in) :: theta0 !Angle between radius to initial fault tip (i.e. lowest position it reaches) and downwards vertical
double precision,intent(in) :: Rc !Radius of the circle
double precision,dimension(ndips),intent(inout) :: sigma_dip_array !Array of sigma dip
double precision :: sigmadip !The uncertainty in the dip being moved.
double precision :: thetat !Angle between radius to the fault tip and downwards vertical
integer :: tip_loc,tip_loc_start !Tells which location the tip is in: 1 = curve, 2 = straight; flat is not allowed
integer :: pt_loc !Tells where the point is: 1 = flat, 2 = curve, 3 = straight, 4 = trishear, 5 = footwall
integer :: slip_sign,prop_sign !slip_sign: 1 = forwards motion, -1 = backwards motion
double precision :: half_v0,xexp,yexp,signy
double precision :: total_slip_abs !absolute value of the total slip
double precision :: slip !the slip produced by a specific step
double precision :: theta_rot !theta that the tip has rotated through.
double precision :: Sremstr !remaining slip while tip is in straight zone
double precision :: Sremcurve !remaining slip while tip is in straight zone
double precision :: total_slip_curve,total_slip_str !Total slip with the tip in each part.
double precision :: xt,yt !Cartesian coordinates of the tip
double precision :: xt0,yt0 !Initial cartesian coordinates of the tip (as in lowest; initial in forwards slip)
double precision :: xchange !x tip coordinate at which it changes from listric to planar
double precision :: vx, vy !x and y components of the velocity
double precision :: x,y,dip !Cartesian coordinates of the point in question, and dip at (x,y)
double precision :: r,theta !Polar coordinates of the point in question
double precision,dimension(2,1) :: pt,ptze !(x,y) in one array, and its trishear equivalent (z,e)
double precision :: ccx,ccy !Coordinates of the center of the circle that the fault forms an arc of
integer :: k !A counter of beds.
double precision :: theta_prop,sin_theta_prop,cos_theta_prop !Theta for one increment of propagation, and its sine and cosine.
double precision :: theta_d,phi_d,gamma1,gamma2,delta !Terms used in calculating the change in dip.
double precision :: v0term,exp1,exp2,sqrtm !Terms for changing dip in the trishear zone.

half_v0 = v0/2; !v0/2
xexp = 1/s; !exponent in vx term
yexp = (1+s)/s; !exponent in vy term
exp1 = (1+s)/(2*s) !exponent in 1st term for change in dip
exp2 = (1-s)/(2*s) !exponent in 2nd term for change in dip
half_v0 = v0/2; !v0/2
v0term = half_v0/s !v0/(2s) constant term that goes in front in dip equation
sqrtm = sqrt(m) !square root of m
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
ccx = cc(1)
ccy = cc(2)
xt0 = ccx+Rc*sin(theta0) !Lowest fault tip coordinates. Note: This is only true if tip does enter curved segment, but if not these values don't get used anyway.
yt0 = ccy-Rc*cos(theta0)
xchange = ccx+Rc*sin(theta_max) !x tip coordinate at which it changes from listric to planar.
theta_prop = PoverS*increment/RC !theta_prop will be in radians; propagation in one increment of slip
sin_theta_prop = sin(theta_prop) !Calculate these once for efficiency
cos_theta_prop = cos(theta_prop)
if (thetat_start >= theta_max) then !Find out which part of the fault the tip starts in
    tip_loc_start = 2 !In the straight part
else
    tip_loc_start = 1 !In the curved part
end if
if (tip_loc_start == 1) then  !Find out how much slip occurs with the tip in each part
    if (slip_sign == prop_sign) then
        total_slip_curve = min((theta_max-theta0)*Rc/PoverS,total_slip)
        total_slip_str = total_slip - total_slip_curve
    else
        total_slip_curve = total_slip
        total_slip_str = 0
    end if
else
    if (slip_sign == prop_sign) then
        total_slip_curve = 0
        total_slip_str = total_slip
    else
        total_slip_str = max(-sqrt((tip_pos(1)-(ccx+Rc*cos(theta_max-pi/2)))**2 + &
            (tip_pos(2)-(ccy+Rc*sin(theta_max-pi/2)))**2)/PoverS,total_slip) !Max b/c the slips should be negative in this case
        total_slip_curve = total_slip - total_slip_str
        if (atan((ccy-Rc*cos(theta_max)-tip_pos(2))/(ccx+Rc*sin(theta_max)-tip_pos(1)))-theta_max >= 1.0e-5) then
            print*,'Error'
            print*,atan((ccy-Rc*cos(theta_max)-tip_pos(2))/(ccx+Rc*sin(theta_max)-tip_pos(1)))
            print*,atan((ccy-Rc*cos(theta_max)-tip_pos(2))/(ccx+Rc*sin(theta_max)-tip_pos(1)))*180/pi
            print*,theta_max
            print*,theta_max*180/pi
        end if
    end if
end if
!print*,total_slip
do k = 1,size(pts,2)
    !print*,'k = ',k
    x = pts(1,k)
    y = pts(2,k)
    xt = tip_pos(1)
    yt = tip_pos(2)
    dip = dips(k)
    sigmadip = sigma_dip_array(k)
    Sremcurve = total_slip_curve
    Sremstr = total_slip_str
    !print*,Sremstr,xt,yt
    !print*,Rc,theta_max,ccx,ccy,PoverS
    slip = 0
    pt_loc = 0 !To start
    tip_loc = tip_loc_start
    thetat = thetat_start
    do while (slip_sign*Sremcurve > 0 .or. slip_sign*Sremstr > 0)
        slip = 0 !Zero it out for next time. (Especially important if pt_loc == 0)
        !print*,pt_loc
        select case(pt_loc)
        case (5) !Point is in the footwall
            call footwall_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,xt0,yt0,theta0,&
                cc,Rc,phi,PoverS,theta_max,prop_sign)
        case(3) !Point is in the straight part
            call straight_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,xt0,yt0,theta0,&
                cc,Rc,phi,m,PoverS,theta_max,prop_sign)
            if (pt_loc == 1) then !Point moves backwards into the flat zone. We need to change the dip.
                phi_d = thetat !Change in ramp angle. Since tip is in straight part and flat is flat, it will be thetat, which is also theta_max if tip is in the straight part.
                theta_d = -dip-thetat
                gamma1 = pi/2-phi_d/2-theta_d
                gamma2 = atan((tan(phi_d/2)-sin(theta_d)/(cos(phi_d/2)*cos(phi_d/2+theta_d)))**(-1))
                if (gamma2 < 0) then !gamma2 should be in the range [0,pi]
                    gamma2 = pi+gamma2
                end if
                delta = pi-gamma1-gamma2 !change in dip
                dip = dip+delta;
                sigmadip = (1.-4.*sin(theta_d)*(-sin(theta_d)+sin(theta_d+phi_d))/ &
                    (-3.+2*cos(2*theta_d)+cos(phi_d)-2*cos(2*theta_d+phi_d)))*sigmadip
            end if
        case(2)!Point is in listric part
            theta = atan((x-cc(1))/(cc(2)-y))!Convert point to polar coordinates, with theta measured from downwards vertical
            r = sqrt((x-cc(1))**2+(y-cc(2))**2)
            call listric_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,cc,Rc,phi,PoverS,theta_max,&
                prop_sign)
            dip = -(theta+slip/r)+atan(1./((-slip/r)+1./tan(theta+dip)))
            sigmadip = sigmadip*(1./sin(theta+dip)**2)/(1.+((-slip/r)+1./tan(theta+dip))**2)
        case(1)!Point is in the flat part.
            call flat_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,cc,Rc,PoverS,theta_max,&
                prop_sign)
            if (pt_loc == 3) then !Point goes into the straight zone. (Slip must be forwards for this to happen.)
                phi_d = thetat!Change in ramp angle; since tip is in straight part and flat is flat, it's thetat.
                theta_d = dip !I think this is right
                gamma1 = pi/2-phi_d/2-theta_d
                gamma2 = atan((tan(phi_d/2)-sin(theta_d)/(cos(phi_d/2)*cos(phi_d/2+theta_d)))**(-1)) ! Precalculating at least tan(phi/2) and cos(phi/2) would help.
                if (gamma2 < 0) then !gamma2 should be in the range [0,pi]
                    gamma2 = pi+gamma2
                end if
                delta = pi-gamma1-gamma2 !change in dip
                dip = dip+delta;
                sigmadip = (1.+4.*sin(theta_d)*(-sin(theta_d)+sin(theta_d+phi_d))/ &
                    (-3.+2*cos(2*theta_d)+cos(phi_d)-2*cos(2*theta_d+phi_d)))*sigmadip
            end if
        case(4) !Point is in the trishear zone
            !Note: Right now I have ramp_dir as 1 in xy_to_ze2D and ze_to_xy2D, assuming any reflecting will be done earlier.
            pt(:,1) = (/x,y/)
            ptze = xy_to_ze2D(pt,xt,yt,thetat,1); !Rotate into trishear coordinates
            x = ptze(1,1)
            y = ptze(2,1)
            if (x < 0) then
                print*,'Error: x < 0'
            end if
            slip = 0
            theta_rot = 0
            do while (slip == 0 .or. (y<=m*x .and. y>=-m*x))
                !A point should only be sent here if it is going to have at least some slip in the trishear zone.
                !However, small rounding errors can make y>m*x or <-m*x when it should be = to one of these.
                !Thus the slip == 0 possibility in the do while loop.
                if (y >= 0) then
                    signy = 1
                else
                    signy = -1
                end if
                !First slip
                slip = slip+increment
                if (tip_loc == 1) then
                    theta_rot = theta_rot + theta_prop !This will keep track of the total rotation of the fault tip
                end if
                dip = dip+(v0term/x)*(signy*sqrtm*((abs(y)/(m*x))**(exp1))*cos(dip)+(1/sqrtm)*((abs(y)/(m*x))**(exp2))*sin(dip))**2
                sigmadip = sigmadip*(1.+signy*(increment/(s*x))*((abs(y)/(m*x))**xexp)*(0.5*((x/y)-(y/x))*sin(2*dip)+cos(2*dip)))
                vx = half_v0*((signy*(abs(y)/(m*x))**xexp)+1)
                vy = half_v0*(m/(1+s))*(((abs(y)/(m*x))**yexp)-1)
                x = x+vx
                y = y+vy
                !Then propagate
                if (tip_loc == 1) then !Combined rotation and translation. Here we rotate and then translate.                                
                    !We're rotating the axis by theta_prop, so the point rotates by negative theta_prop.
                    x = x*cos_theta_prop + y*sin_theta_prop - Rc*sin_theta_prop !Last term will end up being positive if theta_prop is negative, as it is for inversion.
                    y = -x*sin_theta_prop + y*cos_theta_prop - slip_sign*Rc*(1.-cos_theta_prop) !I think slip_sign is needed here, because cos will be positive whether or not theta_prop is.
                else
                    x = x-PoverS*increment
                end if
                if (tip_loc == 1 .and. abs(slip) > abs(Sremcurve)) then !See if we've reached the end of the slip in this section
                    exit
                else if (tip_loc == 2 .and. abs(slip) > abs(Sremstr)) then
                    exit
                end if
            end do
            !print*,slip
            if (tip_loc == 1) then!Update tip location.
                theta_rot = slip*PoverS/Rc !The amount that the tip rotates during this slip
                if (slip_sign == -1 .and. abs(theta_rot) > abs(thetat-theta0)) then !Note that for forwards slip we need something similar if it reaches theta_max
                    theta_rot = theta0-thetat
                end if
                thetat = thetat+theta_rot !Update thetat
                xt = cc(1)+Rc*sin(thetat)
                yt = cc(2)+slip_sign*Rc*cos(thetat)
            else
                xt = xt+PoverS*slip*cos(theta_max)
                yt = yt+PoverS*slip*sin(theta_max)
            end if
            ptze(:,1) = (/x,y/)
            pt = ze_to_xy2D(ptze,xt,yt,thetat,1)!Rotate out of trishear coordinates.
            x = pt(1,1)
            y = pt(2,1)
            pt_loc = 0;
        case default
            !Need to figure out where we are
            if ((y <= ccy-Rc .and. x<=ccx) .or. &
                    (x>=ccx .and. y<=(x-ccx)*tan(thetat-pi/2)+ccy .and. (x-ccx)**2+(y-ccy)**2 >= Rc**2) &
                    .or. (x>=xchange .and. x<=xt .and. y<=(x-xt)*tan(theta_max)+yt) .or. &
                    (x>xt .and. y<=(x-xt)*tan(thetat-phi)+yt)) then
                pt_loc = 5 !Footwall
            else if ((x<=ccx .and. y<=ccy) .or. (y>=ccy .and. y<=-(x-ccx)*tan((pi-thetat)/2)+ccy)) then
                pt_loc = 1 !Flat
            else if (y<=ccy .and. atan((x-ccx)/(ccy-y)) <= thetat) then
                pt_loc = 2 !Curve
            else if (x<=(y-yt)/tan(thetat+phi)+xt) then
                pt_loc = 3 !Straight
            else
                pt_loc = 4 !Trishear
            end if
 !           if (pt_loc /= 4 .or. (Sremcurve==total_slip_curve .and. Sremstr==total_slip_str)) then
 !               print*,pt_loc
 !           end if
        end select
        if (tip_loc == 1) then !Update the remaining amount of slip with the tip in the appropriate segment
            Sremcurve = Sremcurve - slip
            if (slip_sign*Sremcurve <= 0) then
                tip_loc = 2
            end if        
        else
            Sremstr = Sremstr - slip
            if (slip_sign*Sremstr <=0) then
                tip_loc = 1
            end if        
        end if
    end do
    pts(1,k) = x !Store the new x and y position values
    pts(2,k) = y
    dips(k) = dip !Store the new dip to pass back to the rest of the program.
    if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
        sigma_dip_array(k) = sigmadip
    end if
end do
end subroutine trishear_func_dip_listric

subroutine footwall_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,xt0,yt0,theta0,cc,Rc,phi,PoverS,&
    theta_max,prop_sign)
use constants
implicit none
double precision,intent(in) :: x,y !The point
double precision,intent(inout) :: slip,xt,yt,thetat !Definitions as in trishear_func_listric_above
integer,intent(inout) :: pt_loc
integer,intent(in) :: tip_loc,slip_sign,prop_sign
double precision,intent(in) :: Sremcurve,Sremstr
double precision,intent(in) :: xt0,yt0,theta0
double precision,dimension(2),intent(in) :: cc
double precision,intent(in) :: Rc,phi,PoverS,theta_max
double precision :: r,theta !Polar coordinates of the point in question
double precision :: Stri,theta_tri !slip and tip rotation necessary to reach the trishear zone
double precision :: theta_rot !theta that the tip has rotated through.
if (slip_sign == prop_sign) then !Forwards slip
    if (tip_loc == 1) then !Tip is in listric segment
        slip = Sremcurve
    else !Tip is in straight segment
        slip = Sremstr
    end if
else !Backwards slip
    if (tip_loc == 2) then !Tip is in straight segment
        Stri = (y-yt+(xt-x)*tan(thetat-phi))/(PoverS*(sin(thetat)-cos(thetat)*tan(thetat-phi)))
        if (Stri>0 .and. x>cc(1)+Rc*sin(theta_max)) then
            print*,'Stri > 0'
            print*,(y <= cc(2)-Rc .and. x<=cc(1))
            print*,(x>=cc(1) .and. y<=(x-cc(1))*tan(thetat-pi/2)+cc(2) .and. (x-cc(1))**2+(y-cc(2))**2 >= Rc**2)
            print*,(x>=cc(1)+Rc*sin(theta_max) .and. x<=xt .and. y<=(x-xt)*tan(theta_max)+yt)
            print*,(x>xt .and. y<=(x-xt)*tan(thetat-phi)+yt)
            print*,x<=(y-cc(2))/tan(thetat-pi/2)+cc(1)
            print*,-sqrt((xt-(cc(1)+Rc*cos(theta_max-pi/2)))**2 + (yt-(cc(2)+Rc*sin(theta_max-pi/2)))**2)/PoverS
        end if
        if (abs(Stri)<abs(Sremstr) .and. x>cc(1)+Rc*sin(theta_max)) then !Point goes into trishear zone
            slip = Stri
            pt_loc = 4
        else !Tip goes into curved segment (if there is slip where it is in the curved part)
            slip = Sremstr
        end if
    else !Tip is in curved segment
        if (x<xt0 .or. (y-yt0) < (x-xt0)*tan(theta0-phi)) then !Point stays in footwall. There can be points in the footwall above the line of the trishear zone boundary, but they will be at x<xt.
            slip = Sremcurve
        else !Point goes into trishear zone
            theta = atan2((x-cc(1)),(cc(2)-y))!Convert to polar coordinates, with theta measured from downwards vertical
            r = sqrt((x-cc(1))**2+(y-cc(2))**2)
            theta_tri = thetat - (-acos((Rc/r)*cos(phi))+theta+phi)
            if (theta_tri < 0) then
                print*,'theta_tri<0'
                print*,Sremcurve,Sremstr
                print*,theta_tri,thetat,theta0
                print*,acos((Rc/r)*cos(phi))+theta+phi
                print*,thetat - (-acos((Rc/r)*cos(phi))+theta+phi),thetat - (acos((Rc/r)*cos(phi))+theta+phi)
                print*,y-yt,(x-xt)*tan(thetat-phi),(y-yt0),(x-xt0)*tan(theta0-phi)
            end if
            if (theta_tri > thetat-theta0) then
                print*,'Error: theta_tri > thetat-theta0'
                print*,theta_tri,thetat,theta0
                print*,thetat - (-acos((Rc/r)*cos(phi))+theta+phi),thetat - (acos((Rc/r)*cos(phi))+theta+phi)
            end if
            slip = -theta_tri*Rc/PoverS !Negative since backwards slip
            pt_loc = 4
        end if
    end if
end if
if (pt_loc==4 .and. x<xt+PoverS*slip*cos(theta_max) .and. tip_loc==2) then
    print*,(y <= cc(2)-Rc .and. x<=cc(1))
    print*,(x>=cc(1) .and. y<=(x-cc(1))*tan(thetat-pi/2)+cc(2) .and. (x-cc(1))**2+(y-cc(2))**2 >= Rc**2)
    print*,(x>=cc(1)+Rc*sin(theta_max) .and. x<=xt .and. y<=(x-xt)*tan(theta_max)+yt)
    print*,(x>xt .and. y<=(x-xt)*tan(thetat-phi)+yt)
    print*,x<=(y-cc(2))/tan(thetat-pi/2)+cc(1)
    print*,x<=(y-(yt+PoverS*Sremstr*sin(thetat)))/tan(thetat-pi/2)+xt+PoverS*Sremstr*cos(thetat)
    print*,xt+PoverS*slip*cos(theta_max)
end if
if (tip_loc == 1) then!Update tip location.
    theta_rot = slip*PoverS/Rc !The amount that the tip rotates during this slip
    thetat = thetat+theta_rot !Update thetat
    xt = cc(1)+Rc*sin(thetat)
    yt = cc(2)-Rc*cos(thetat)
else
    xt = xt+PoverS*slip*cos(theta_max)
    yt = yt+PoverS*slip*sin(theta_max)
end if
if (pt_loc==4 .and. x<xt) then
    print*,'slip = ',slip
    if (tip_loc == 2) then
        print*,'Stri = ',Stri,'Sremstr = ',Sremstr
    end if
    print*,x-xt,(thetat-phi)*180/pi
    print*,'tip_loc = ',tip_loc
    !I think that at least when tip_loc = 2, there is a problem where the point shouldn't be in the footwall.
    print*,(x-xt)*tan(thetat)
    print*,(y-yt),(x-xt)*tan(thetat-phi)
    print*,(x-cc(1))**2+(y-cc(2))**2,Rc**2
end if
if (pt_loc == 4 .and. (abs(x-((y-yt)/tan(thetat+phi)+xt))>1.0e-5) .and. (abs(y-((x-xt)*tan(thetat-phi)+yt))>1.0e-5)) then !Allows small errors due to rounding.
    print*,'Error: Should not be in trishear zone.'
    print*,'From footwall'
    print*,y-yt,(x-xt)*tan(thetat+phi)
    print*,y-yt,(x-xt)*tan(thetat-phi)
end if
if (pt_loc == 4 .and. x<(y-yt)/tan(thetat-pi/2)+xt) then
    print*,'Error: Point should not be going into the trishear zone.'
    print*,'From footwall'
end if
!print*,pt_loc
end subroutine footwall_slip_listric

subroutine straight_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,xt0,yt0,theta0,cc,Rc,phi,m,&
    PoverS,theta_max,prop_sign)
use constants
use parameters,only: increment
implicit none
double precision,intent(inout) :: x,y !The point
double precision,intent(inout) :: slip,xt,yt,thetat !Definitions as in trishear_func_listric_above
integer,intent(inout) :: pt_loc
integer,intent(in) :: tip_loc,slip_sign,prop_sign
double precision,intent(in) :: Sremcurve,Sremstr
double precision,intent(in) :: xt0,yt0,theta0
double precision,dimension(2),intent(in) :: cc
double precision,intent(in) :: Rc,phi,m,PoverS,theta_max
double precision :: b !A variable used in calculating Stri for tip in curved segment, point in straight
double precision :: ccpx,ccpy !Coordinates of the center of a smaller circle that a point in the straight zone rotates around when tip is in listric zone.
double precision :: Sflat,Slist,Stri !Slip needed to move into different zones.
double precision :: theta_rot !theta that the tip has rotated through.
double precision :: f,flast !Function to be minimized by Newton-Raphson method and its previous iteration value.
integer :: niter !Number of iterations performed by the Newton-Raphson method.
if (tip_loc == 2) then !Tip in straight part too
    if (thetat+phi .ne. pi/2) then !Since in this case tan(thetat+phi) would be infinity
        Stri = ((yt-y)-(xt-x)*tan(thetat+phi))/((1-PoverS)*(sin(thetat)-cos(thetat)*tan(thetat+phi)))
    else
        Stri = (xt-x)/((1-PoverS)*cos(thetat))
    end if
    if (slip_sign == prop_sign) then !Forwards slip
        if (abs(Stri)<=abs(Sremstr) .and. Stri>0) then !Point goes into trishear zone
            pt_loc = 4
            slip = Stri
        else !Point statys in the straight part
            slip = Sremstr
        end if
    else !Backwards slip
        if (abs(Stri)<=abs(Sremstr) .and. Stri<0) then !Point goes into trishear zone
            pt_loc = 4
            slip = Stri
        else if (y<(x-cc(1))*tan(thetat)+cc(2)) then !Point moving towards listric zone
            Slist = -(x-cc(1))*cos(thetat)-(y-cc(2))*sin(thetat)
            if (abs(Slist)<=abs(Sremstr)) then !Point goes into listric zone
                pt_loc = 2
                slip = Slist
            else
                slip = Sremstr
            end if
        else !Point moving towards flat zone
            Sflat = (cc(2)-y)*tan(thetat/2)+(cc(1)-x)
            if (abs(Sflat) <= abs(Sremstr)) then
                pt_loc = 1
                slip = Sflat
            else
                slip = Sremstr
            end if
        end if
    end if
    x = x+slip*cos(theta_max) !Update point position. Uncertainty doesn't change, given these equations.
    y = y+slip*sin(theta_max)
    xt = xt+PoverS*slip*cos(theta_max) !Update tip position
    yt = yt+PoverS*slip*sin(theta_max)
else !Tip is in curved segment
    if (slip_sign == prop_sign) then !Forwards slip
        if (x+(Rc/PoverS)*(sin(theta_max)-sin(thetat))-xt0 > (y+(Rc/PoverS)*(-cos(theta_max)+cos(thetat)))/tan(theta_max+phi)) then !Point goes into trishear zone
            !I think there may be an analytic way to solve for Stri, but I haven't had much luck figuring it out.
            theta_rot = -Sremcurve*PoverS/Rc !Initial guess
            niter = 0 !Number of iterations
            do
                f = (y-cc(2)+Rc*(1.-1./PoverS)*cos(theta_rot+thetat)+(Rc/PoverS)*cos(thetat))&
                    -(x-cc(1)+Rc*((1./PoverS)-1.)*sin(theta_rot+thetat)-(Rc/PoverS)*sin(thetat))*tan(theta_rot+thetat+phi)
                if (niter == 1e4) then
                    print*,'Exceeded maximum number of iterations to find Stri'
                    print*,'Error = ',f
                    exit
                end if
                if (abs(f) < abs(increment)) then !Close enough.
                    exit
                else !Use Newton-Raphson method
                    theta_rot = theta_rot-(2*cos(theta_rot+thetat+phi)*((PoverS-1)*Rc*cos(phi)+Rc*cos(theta_rot+phi)+&
                        PoverS*((y-cc(2))*cos(theta_rot+thetat+phi)+(cc(1)-x)*sin(theta_rot+thetat+phi))))/&
                        (2*cc(1)*PoverS-2*PoverS*x+2*Rc*sin(thetat)+(PoverS-1)*Rc*sin(theta_rot+thetat)-&
                        Rc*sin(theta_rot+thetat+2*phi)+PoverS*Rc*sin(theta_rot+thetat+2*phi))
                end if
                niter = niter+1
            end do
            Stri = theta_rot*Rc/PoverS
            pt_loc = 4
            slip = Stri
        else
            Slist = (Rc/PoverS)*(atan2(x-cc(1)-(Rc/PoverS)*sin(thetat),cc(2)-y)-thetat)
            if (abs(Slist)<abs(Sremcurve)) then !Point goes into listric zone
                pt_loc = 2
                slip = Slist
            else !Point stays in straight zone
                slip = Sremcurve
            end if
        end if
    else !Backwards slip
        if (PoverS > 1. .and. x+(Rc/PoverS)*(sin(theta0)-sin(thetat))-xt0 > &
                (y+(Rc/PoverS)*(-cos(theta0)+cos(thetat))-yt0)/tan(theta0+phi)) then !Point can go into trishear zone
            !I think there may be an analytic way to solve for Stri, but I haven't been able to figure it out.
            theta_rot = Sremcurve*PoverS/Rc !Initial guess
            niter = 0 !Number of iterations
            do
                f = (y-cc(2)+Rc*(1.-1./PoverS)*cos(theta_rot+thetat)+(Rc/PoverS)*cos(thetat))&
                    -(x-cc(1)+Rc*((1./PoverS)-1.)*sin(theta_rot+thetat)-(Rc/PoverS)*sin(thetat))*tan(theta_rot+thetat+phi)
                if (niter >= 1e3) then
                    exit
                end if
                if (abs(f) < abs(increment)) then !Close enough.
                    if (abs(theta_rot)<abs(theta0-thetat)) then !Within acceptable range
                        exit
                    else !May be a solution, but it's not the solution we want.
                        theta_rot = theta_rot/2 !Try something else instead.
                    end if
                else !Use Newton-Raphson method
                    theta_rot = theta_rot-(2*cos(theta_rot+thetat+phi)*((PoverS-1)*Rc*cos(phi)+Rc*cos(theta_rot+phi)+&
                        PoverS*((y-cc(2))*cos(theta_rot+thetat+phi)+(cc(1)-x)*sin(theta_rot+thetat+phi))))/&
                        (2*cc(1)*PoverS-2*PoverS*x+2*Rc*sin(thetat)+(PoverS-1)*Rc*sin(theta_rot+thetat)-&
                        Rc*sin(theta_rot+thetat+2*phi)+PoverS*Rc*sin(theta_rot+thetat+2*phi))
                end if
                if (abs(theta_rot) >= 2*pi) then
                    theta_rot = mod(theta_rot,2*pi) !Extra factors of 2*pi will mess things up.
                end if
                if (theta_rot > 0) then
                    theta_rot = theta_rot - 2*pi !Should give same f but a theta_rot value that is negative as necessary for backwards slip
                end if
                niter = niter+1
            end do
            if (abs(f) > abs(increment) .or. theta_rot > 0 .or. theta_rot <theta0-thetat) then
                !print*,'Newton-Raphson answer not acceptable. Trying iterative slip.'
                Stri = 0.
                do while (Stri>Sremcurve) !Greater than because they both should be negative
                    Stri = Stri+increment !Increment it forwards.
                    theta_rot = PoverS*Stri/Rc
                    if ((x-cc(1)+Rc*((1./PoverS)-1.)*sin(theta_rot+thetat)-(Rc/PoverS)*sin(thetat)) > &
                        (y-cc(2)+Rc*(1.-1./PoverS)*cos(theta_rot+thetat)+(Rc/PoverS)*cos(thetat))/tan(theta_rot+thetat+phi)) then !Then you've entered the trishear zone
                        !print*,(y-cc(2)+Rc*(1.-1./PoverS)*cos(theta_rot+thetat)+(Rc/PoverS)*cos(thetat))&
                        !    -(x-cc(1)+Rc*((1./PoverS)-1.)*sin(theta_rot+thetat)-(Rc/PoverS)*sin(thetat))*tan(theta_rot+thetat+phi)
                        exit
                    end if
                    if (Stri <= Sremcurve) then
                        print*,'Error: Point does not enter trishear zone.'
                        exit
                    end if
                end do
            else
                Stri = theta_rot*Rc/PoverS !Newton-Raphson answer was good.
            end if
            Slist = (Rc/PoverS)*(atan2(x-cc(1)-(Rc/PoverS)*sin(thetat),cc(2)-y)-thetat) !Slip to go into listric zone. If this happens first, then we use it, not Stri. Note: For a given point, I think only one is possible. There might be a way to figure out which is which. However, points that pass through list-straight boundary can end up hitting back part of trishear zone.
            Sflat = (Rc/PoverS)*(2*atan2((Rc/PoverS)*(sin(thetat)-1)-x+cc(1),(Rc/PoverS)*cos(thetat)+y-cc(2))-thetat)
            if ((Stri >= Slist .or. Slist>0.) .and. (Stri>=Sflat .or. Sflat>0.)) then !Greater than because they should both be negative
                pt_loc = 4 !Goes into trishear zone
                slip = Stri
            else if (Slist > Sflat .or. Sflat>0.) then
                pt_loc = 2 !Goes into listric zone
                slip = Slist
                !print*,'a'
            else
                pt_loc = 1 !Goes into flat zone
                slip = Sflat
            end if
        else
            Slist = (Rc/PoverS)*(atan2(x-cc(1)-(Rc/PoverS)*sin(thetat),cc(2)-y)-thetat)
            Sflat = (Rc/PoverS)*(2*atan2((Rc/PoverS)*(sin(thetat)-1)-x+cc(1),(Rc/PoverS)*cos(thetat)+y-cc(2))-thetat)
            if (abs(Slist) < abs(Sremcurve) .and. (abs(Slist) < abs(Sflat) .or. Sflat>0.) .and. Slist<0.) then !Point goes into listric zone
                pt_loc = 2
                slip = Slist
                !print*,'b'
            else                
                if (abs(Sflat) < abs(Sremcurve) .and. Sflat<0.) then !Point goes into flat zone
                    pt_loc = 1
                    slip = Sflat
                else !Point stays in straight zone
                    slip = Sremcurve
                end if
            end if
        end if
    end if
    theta_rot = slip*PoverS/Rc !The amount that the tip rotates during this slip
    x = x+(Rc/PoverS)*(sin(theta_rot+thetat)-sin(thetat)) !Update point position. Again no change in uncertainty.
    y = y+(Rc/PoverS)*(-cos(theta_rot+thetat)+cos(thetat))
    thetat = thetat+theta_rot !Update thetat
    xt = cc(1)+Rc*sin(thetat)
    yt = cc(2)-Rc*cos(thetat)
end if
if (pt_loc == 4 .and. (abs(x-((y-yt)/tan(thetat+phi)+xt))>abs(increment*PoverS)) .and. &
        (abs(y-((x-xt)*tan(thetat-phi)+yt))>abs(increment*PoverS))) then !Allows small errors due to rounding and due to imprecision of methods above that use increment.
    print*,'Error: Should not be in trishear zone.'
    print*,'From straight'
    print*,'Error = ',(y-yt)/tan(thetat+phi)+xt-x
    print*,'or', (y-yt)-(x-xt)*tan(thetat+phi)
    print*,increment*PoverS
    print*,Stri
    print*,PoverS
end if
if (pt_loc == 4 .and. x<(y-yt)/tan(thetat-pi/2)+xt) then !Means point is hitting backward extension of trishear zone.
    print*,'Error: Point should not be going into the trishear zone.'
    print*,'From straight'
    print*,Stri,Slist
    print*,x,(y-yt)/tan(thetat-pi/2)+xt
end if
if (slip > 0) then
    print*,'Error: slip>0'
end if
if (pt_loc == 2 .and. x<cc(1)) then
    print*,'Error: Point is in flat zone, not listric zone.'
end if
!print*,pt_loc
end subroutine straight_slip_listric

subroutine listric_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,cc,Rc,phi,PoverS,theta_max,&
    prop_sign,C)
use constants
implicit none
double precision,intent(inout) :: x,y !The point
double precision,intent(inout) :: slip,xt,yt,thetat !Definitions as in trishear_func_listric_above
integer,intent(inout) :: pt_loc
integer,intent(in) :: tip_loc,slip_sign,prop_sign
double precision,intent(in) :: Sremcurve,Sremstr
double precision,dimension(2),intent(in) :: cc
double precision,intent(in) :: Rc,phi,PoverS,theta_max
double precision,dimension(2,2),intent(inout),optional :: C !Covariance matrix for uncertainties in x and y positions.
double precision,dimension(2,2) :: J !Jacobian matrix for propagating uncertainties. Will be used for more than one function.
double precision :: sigma_theta,sigma_r !Corresponding uncertainty in polar coordinates r and theta.
double precision :: r,theta !Polar coordinates of the point in question
double precision :: Sstr,Sflat !Slip needed to move into different zones.
double precision :: theta_rot !theta that the tip has rotated through.
theta = atan((x-cc(1))/(cc(2)-y))!Convert to polar coordinates, with theta measured from downwards vertical
r = sqrt((x-cc(1))**2+(cc(2)-y)**2)
if (theta < 0) then
    print*,'Error: Should not be in listric zone. Theta < 0.'
end if
if (present(C)) then
    J(1,1) = (x-cc(1))/r !dr/dx
    J(1,2) = -(cc(2)-y)/r !dr/dy
    J(2,1) = (cc(2)-y)/(r**2) !dtheta/dx
    J(2,2) = (x-cc(1))/(r**2) !dtheta/dy
    C = matmul(matmul(J,C),transpose(J)) !C is now for covariance of (r,theta)
end if
if (slip_sign == prop_sign) then !Forwards slip
    if (tip_loc == 1) then !Tip is in the listric part too
        Sstr = (thetat-theta)/(1/r-PoverS/Rc)
        if (abs(Sstr) <= abs(Sremcurve)) then !Point goes into straight zone
            pt_loc = 3
            slip = Sstr
        else !Point stays in straight zone
            slip = Sremcurve
        end if
    else !Tip is in straight part
        Sstr = r*(thetat - theta)
        if (abs(Sstr) <= abs(Sremstr)) then
            pt_loc = 3
            slip = Sstr
        else
            slip = Sremstr
        end if
    end if
else !Backwards slip
    if (tip_loc == 1) then !Tip is in the listric part too
        Sstr = (thetat-theta)/(1/r-PoverS/Rc)
        if (abs(Sstr) <= abs(Sremcurve) .and. Sstr<0.) then !Point goes into straight zone
            pt_loc = 3
            slip = Sstr
        else
            Sflat = -r*theta !Slip to get into flat part
            if (abs(Sflat) <= abs(Sremcurve)) then !Point goes into flat part
                pt_loc = 1
                slip = Sflat
            else !Point stays in listric zone
                slip = Sremcurve
            end if
        end if
    else !Tip is in the straight part
        Sflat = -r*theta !Slip to get into flat part
        if (abs(Sflat) <= abs(Sremstr)) then !Point goes into flat part
            pt_loc = 1
            slip = Sflat
        else !Point stays in listric zone. Tip goes into listric zone.
            slip = Sremstr
        end if
    end if
end if
theta = theta+slip/r !Update point position
x = r*sin(theta)+cc(1)
y = -r*cos(theta)+cc(2)
if (present(C)) then
    J = reshape((/0.0d0,-slip/(r**2),0.0d0,1.0d0/),(/2,2/)) !(dr/dr, dtheta/dr, dr/dtheta, dtheta/dtheta) go into J(1,1), J(1,2), J(2,1), and J(2,2) respectively.
    C = matmul(matmul(J,C),transpose(J)) !Propagate the errors in theta = theta+slip/r.
    J(1,1) = sin(theta) !dx/dr
    J(1,2) = r*cos(theta) !dx/dtheta
    J(2,1) = -cos(theta) !dy/dr
    J(2,2) = r*sin(theta) !dy/dtheta
    C = matmul(matmul(J,C),transpose(J)) !Propagate the errors in conversion from (r,theta) to (x,y).
    !Note: The new (x,y) could be written as a single pair of functions of old (x,y) and the Jacobian matrix found for that, so I'd only do J*C*J^T once.
end if
if (tip_loc == 1) then!Update tip location.
    theta_rot = slip*PoverS/Rc !The amount that the tip rotates during this slip
    thetat = thetat+theta_rot !Update thetat
    xt = xt+Rc*sin(theta_rot)
    yt = yt+slip_sign*Rc*(1-cos(theta_rot))
else
    xt = xt+PoverS*slip*cos(theta_max)
    yt = yt+PoverS*slip*sin(theta_max)
end if
if (pt_loc == 4 .and. abs(x-((y-yt)/tan(thetat+phi)+xt))>1.0e-5 .and. abs(y-((x-xt)*tan(thetat-phi)+yt))>1.0e-5) then
    print*,'Error: Should not be in trishear zone.'
    print*,'From listric'
end if
if (pt_loc == 4 .and. x<(y-yt)/tan(thetat-pi/2)+xt) then
    print*,'Error: Point should not be going into the trishear zone.'
    print*,'From listric'
end if
if (slip_sign == -1 .and. slip > 0) then
    print*,'Error: slip > 0'
end if
!print*,pt_loc
end subroutine listric_slip_listric

subroutine flat_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,cc,Rc,PoverS,theta_max,prop_sign)
use parameters,only:increment
use constants
implicit none
double precision,intent(inout) :: x,y !The point
double precision,intent(inout) :: slip,xt,yt,thetat !Definitions as in trishear_func_listric_above
integer,intent(inout) :: pt_loc
integer,intent(in) :: tip_loc,slip_sign,prop_sign
double precision,intent(in) :: Sremcurve,Sremstr
double precision,dimension(2),intent(in) :: cc
double precision,intent(in) :: Rc,PoverS,theta_max
double precision :: Slist,Sstr !Slip needed to move into different zones.
double precision :: theta_rot !theta that the tip has rotated through.
if (slip_sign == prop_sign) then !Forwards slip
    if (y < cc(2)) then !Heading towards listric zone
        Slist = cc(1) - x !Slip needed to enter listric zone
        if (tip_loc == 1) then !Tip is in the listric part                                    
            if (Slist <= Sremcurve) then !Point goes into listric zone
                slip = Slist
                pt_loc = 2
            else !Point stays in flat, tip completes propagation in listric zone
                slip = Sremcurve
            end if
        else !Tip is in the straight part
            if (Slist <= Sstr) then !Point goes into listric zone
                slip = Slist
                pt_loc = 2
            else !Point stays in flat, tip completes propagation in straight zone
                slip = Sremcurve
            end if
        end if
    else !Heading towards kink axis and straight zone
        if (tip_loc == 1) then !Tip is in the listric zone
            if (y > -(x-cc(1)+Sremcurve)*tan((pi-theta_max)/2)) then !Goes into straight zone
                Sstr = increment !To initialize it
                do while (Sstr <= Sremcurve) !This probably isn't the fastest way to solve this
                    if (y >= -(x-cc(1)+Sstr)*tan((pi-(thetat+Sstr*PoverS/Rc))/2)) then
                        exit !This is the (approximate) slip we want          
                    end if
                    Sstr = Sstr+increment
                end do
                slip = Sstr
                pt_loc = 3
            else !Stays in flat zone
                slip = Sremcurve
            end if
        else !Tip is in the straight zone
            Sstr = -(y-cc(2))/tan((pi-theta_max)/2)+cc(1)-x
            if (Sstr <= Sremstr) then !Point goes into the straight zone
                slip = Sstr
                pt_loc = 3
            else !Point stays in the flat zone
                slip = Sremstr
            end if
        end if
    end if
else !Backwards slip
    if (tip_loc == 1) then !Tip is in the listric part
        slip = Sremcurve
    else !Tip is in the straight part
        slip = Sremstr
    end if
end if
x = x+slip !Move point. y stays the same.
if (tip_loc == 1) then!Update tip location.
    theta_rot = slip*PoverS/Rc !The amount that the tip rotates during this slip
    thetat = thetat+theta_rot !Update thetat
    xt = xt+Rc*sin(theta_rot)
    yt = yt+slip_sign*Rc*(1-cos(theta_rot))
else
    xt = xt+PoverS*slip*cos(theta_max)
    yt = yt+PoverS*slip*sin(theta_max)
end if
!print*,pt_loc
if (pt_loc == 4 .and. x<(y-yt)/tan(thetat-pi/2)+xt) then
    print*,'Error: Point should not be going into the trishear zone.'
    print*,'From flat'
end if
end subroutine flat_slip_listric

subroutine trishear_slip_listric(x,y,slip,xt,yt,thetat,pt_loc,tip_loc,slip_sign,Sremcurve,Sremstr,theta0,cc,Rc,phi,m,PoverS,&
    theta_max,v0,s,theta_prop,sin_theta_prop,cos_theta_prop,C)
use math,only:xy_to_ze2D,ze_to_xy2D,rot2D
use parameters,only:increment
use constants
implicit none
double precision,intent(inout) :: x,y !The point
double precision,intent(inout) :: slip,xt,yt,thetat !Definitions as in trishear_func_listric_above
integer,intent(inout) :: pt_loc
integer,intent(in) :: tip_loc,slip_sign
double precision,intent(in) :: Sremcurve,Sremstr
double precision,intent(in) :: theta0
double precision,dimension(2),intent(in) :: cc
double precision,intent(in) :: Rc,phi,m,PoverS,theta_max,v0,s
double precision,dimension(2,2),intent(inout),optional :: C !Covariance matrix for uncertainties in x and y positions.
double precision,dimension(2,2) :: J !Jacobian matrix for propagating uncertainties.
double precision,dimension(2,2) :: R !The rotation matrix for a rotation into the trishear coordinate system.
double precision,dimension(2,1) :: pt,ptze !(x,y) in one array, and its trishear equivalent (z,e)
double precision :: theta_prop,sin_theta_prop,cos_theta_prop !Theta for one increment of propagation, and its sine and cosine.
double precision :: half_v0,xexp,yexp,signy
double precision :: vx, vy !x and y components of the velocity
double precision :: theta_rot !theta that the tip has rotated through.
half_v0 = v0/2; !v0/2
xexp = 1/s; !exponent in vx term
yexp = (1+s)/s; !exponent in vy term
!Note: Right now I have ramp_dir as 1 in xy_to_ze2D and ze_to_xy2D, assuming any reflecting will be done earlier.
pt(:,1) = (/x,y/)
ptze = xy_to_ze2D(pt,xt,yt,thetat,1); !Rotate into trishear coordinates.
x = ptze(1,1)
y = ptze(2,1)
R = reshape((/cos(thetat),sin(thetat),-sin(thetat),cos(thetat)/),(/2,2/))
C = matmul(matmul(transpose(R),C),R) !Change C into the rotated (zeta,eta) coordinate system.
slip = 0
theta_rot = 0
do while (slip == 0 .or. (y<=m*x .and. y>=-m*x)) 
    !A point should only be sent here if it is going to have at least some slip in the trishear zone.
    !However, small rounding errors can make y>m*x or <-m*x when it should be = to one of these.
    !Thus the slip == 0 possibility in the do while loop.
    if (x<0) then
        print*,'Error: x<0'
        print*,'x = ',x
    end if
    if (y >= 0) then
        signy = 1
    else
        signy = -1
    end if
    !First slip
    slip = slip+increment
    if (tip_loc == 1) then
        theta_rot = theta_rot + theta_prop !This will keep track of the total rotation of the fault tip
    end if
    if (y>=m*x+1. .or. y<=-m*x-1.) then
        print*,'Error: Point should not be in the trishear zone.'
        print*,y,m*x,-m*x
    end if
    vx = half_v0*((signy*(abs(y)/(m*x))**xexp)+1)
    vy = half_v0*(m/(1+s))*(((abs(y)/(m*x))**yexp)-1)
    J(1,1) = 1.-(increment/(2*s*x))*signy*(abs(y)/(m*x))**xexp !1+dvx/dx
    J(1,2) = (increment/(2*m*s*x))*(abs(y)/(m*x))**((1.-s)/s) !dvx/dy
    J(2,1) = -((increment*m)/(2*s*x))*(abs(y)/(m*x))**yexp !dvy/dx
    J(2,2) = 1.+(increment/(2*s*x))*signy*(abs(y)/(m*x))**xexp !1+dvy/dy
    C = matmul(matmul(J,C),transpose(J))
    x = x+vx
    y = y+vy
    !Then propagate
    if (tip_loc == 1) then !Combined rotation and translation. Here we rotate and then translate.
        !We're rotating the axis by theta_prop, so the point rotates by negative theta_prop.
        x = x*cos_theta_prop + y*sin_theta_prop - Rc*sin_theta_prop !Last term will end up being positive if theta_prop is negative, as it is for inversion.
        y = -x*sin_theta_prop + y*cos_theta_prop - slip_sign*Rc*(1.-cos_theta_prop) !I think slip_sign is needed here, because cos will be positive whether or not theta_prop is.
        R = reshape((/cos(theta_prop),-sin(theta_prop),sin(theta_prop),cos(theta_prop)/),(/2,2/))
        C = matmul(matmul(R,C),transpose(R)) !Transform C into the zeta, eta coordinate system after propagation and rotation.
    else
        x = x-PoverS*increment
    end if
    if (tip_loc == 1 .and. abs(slip) > abs(Sremcurve)) then !See if we've reached the end of the slip in this section
        exit
    else if (tip_loc == 2 .and. abs(slip) > abs(Sremstr)) then
        exit
    end if
end do
if (tip_loc == 1) then!Update tip location.
    theta_rot = slip*PoverS/Rc !The amount that the tip rotates during this slip
    if (slip_sign == -1 .and. abs(theta_rot) > abs(thetat-theta0)) then !Note that for forwards slip we need something similar if it reaches theta_max
        theta_rot = theta0-thetat
    end if
    thetat = thetat+theta_rot !Update thetat
    xt = cc(1)+Rc*sin(thetat)
    yt = cc(2)+slip_sign*Rc*cos(thetat)
else
    xt = xt+PoverS*slip*cos(theta_max)
    yt = yt+PoverS*slip*sin(theta_max)
end if
ptze(:,1) = (/x,y/)
pt = ze_to_xy2D(ptze,xt,yt,thetat,1)!Rotate out of trishear coordinates.
x = pt(1,1)
y = pt(2,1)
R = reshape((/cos(thetat),sin(thetat),-sin(thetat),cos(thetat)/),(/2,2/))
C = matmul(matmul(R,C),transpose(R)) !Transform C into the zeta, eta coordinate system after propagation and rotation.
pt_loc = 0;
end subroutine trishear_slip_listric

subroutine fault_pts_err_listric(pts,theta,cc,Rc,errs)
use data_module, only: nfaultpts
implicit none
double precision,dimension(1:,1:),intent(in) :: pts !The points on the fault
double precision,intent(in) :: theta !The angle that the fault tip makes with the vertical (in its final position)
double precision,dimension(2),intent(in) :: cc !The center of the circle.
double precision,intent(in) :: Rc !The radius of the circle
double precision,dimension(nfaultpts),intent(out) :: errs !The errors for each point
double precision :: x,y !Coordinates of each point
double precision :: tan_theta,tan_half_theta !Calculate these once for speed.
integer :: n !a counter
tan_theta = tan(theta) !Calculate once for simplicity
tan_half_theta = tan(theta/2.)
do n = 1,nfaultpts
    x = pts(1,n)
    y = pts(2,n)
    !Need to figure out which part of the fault it is closest to.
    !Note: This assumes fault verges right. If not, it must be reflected ahead of time.
    if ((x < cc(1) .and. y < cc(2)) .or. (y > cc(2) .and. y <= -tan_half_theta*(x-cc(1))+cc(2))) then !Closest to the flat part
        errs(n) = abs(y-(cc(2)-Rc))
    else if (y < cc(2) .and. y < -tan_theta*(x-cc(1))+cc(2)) then !Closest to the listric part
        errs(n) = abs(Rc-sqrt((x-cc(1))**2+(y-cc(2))**2))
    else !Closest to the straight part
        errs(n) = abs(((x-cc(1))*tan_theta-y+cc(2))/sqrt(1.+tan_theta**2))
    end if
end do
end subroutine fault_pts_err_listric

end module trishear_listric