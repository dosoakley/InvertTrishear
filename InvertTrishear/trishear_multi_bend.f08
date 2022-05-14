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

module trishear_multibend
!Trishear for a fault with multiple bends in it.
!It is assumed that the higher segment must dip more steeply than the lower one (This might not always have to be true.)
!and that both segments must slope  in the same direction.
integer :: bend_at_tipy !Tells whether one of the bend depths and tipy must be the same (1) or not (0)
integer :: bend_start !The bend at which the tip must start if bend_at_tipy == 1
!integer :: tipy_can_be_below_bend !Tells whether tipy (initial) can be below the bend
integer :: tipy_above_bend !Tells whether the initial tipy is forced to be above one of the bends, which is specified.
integer :: bend_tip_above !The number of the bend that the tip must be above if tipy_above_bend == 1.
integer :: nbends !Tells the total number of fault bends / points of change in ramp angle, P/S, or phi.
integer :: nsegs !Tells the total number of fault segments. This will be nbends+1
integer :: nangles,nPoverS,nphi !Numbers of different ramp angle, PoverS, and phi values. The sum of these three is nsegs.
integer :: backlimb_fold_type !(1) Fault parallel flow, (2) Fault bend folding, (3) Inclined Shear
integer :: phi_le_ramp !Tells whether or not phi is required to be less than the ramp angle in any segment that the tip passes through.
integer :: let_trishear_intersect_axes !Tells whether the trishear zone is allowed to intersect the fold axes above some specified elevation. (1 = yes, 0 = no)
double precision :: min_intersect_elevation !If trishear zone is allowed to intersect fold axes above some elevation, this is the elevation
integer :: final_tip_in_upper_segment !Requires the final position of the fault tip to be in the uppermost fault segment, so that there are no imaginary, unused fault segments above it.
integer :: limit_fault_x_pts !Tells whether or not to put limits on the x positions of fault points.
double precision,dimension(2) :: fault_x_limits !Minimum and maximum allowed value for fault points.
type multibend_fault
    double precision :: tipxinit, tipyinit, tipxfinal, tipyfinal !Initial and final (in terms of forward motion) tip positions
    integer :: tip_seg_init,tip_seg_final !tells which segment of the fault the fault tip is in. 1 is top segment.
    double precision,dimension(:),allocatable :: phi,m,PoverS !These values are now allowed to change.
    double precision :: s !Concentration factor.
    integer :: ramp_dir !dip direction of ramp: 1 for a fault that dips left, -1 for a fault that dips right.
    double precision,dimension(:,:),allocatable :: bendxy !x,y coordinates of the fault bends. Numbered from highest = 1 to lowest = nbends.
    double precision,dimension(:),allocatable :: ramp_angle !1st value is highest segment.
    double precision,dimension(:),allocatable :: axis_dips !Fold axis dips. (90 deg / pi/2 rad is vertical, 0 is horizontal.)
    double precision,dimension(:),allocatable :: R !Ratio of slip on a given segment relative to slip on the lowest segment.
    double precision,dimension(:),allocatable :: slipseg !The amoung of slip with the tip in each segment. (Absolute value.) In the order that the tip moves through them.
    double precision,dimension(:),allocatable :: growth_slip !Slip needed to restore each growth strata bed.
    double precision,dimension(:),allocatable :: terr_slip !Slip needed to restore each terrace.
    double precision,dimension(:),allocatable :: restored_fit_params
end type multibend_fault
integer :: constrain_backlimb_syncline !Tells whether to put constraints on x position of backlimb syncline (0 or 1).
double precision :: y_backlimb_syncline !Elevation at which to put the constraints on the x position of the backlimb syncline.
double precision,dimension(2) :: x_backlimb_syncline_limits !Minimum and maximum allowed x positions of backlimb syncline at y = y_backlimb_syncline
integer :: fault_concave_up !Tells whether the fault is required to be concave up (1) or can have convex bends (0).

contains

subroutine trishear_multibend_options
use parameters, only: nparams
use options, only: options_id,TipToSolve,FaultType
implicit none
character (len = 30) :: str

nparams = 4 !To start: tipx, tipy, slip, and s.
do
    print*,'Enter number of fault segments with different ramp angles'
    read(options_id,*) nangles
    print*,'Enter number of different phi values'
    read(options_id,*) nphi
    print*,'Enter number of different P/S values'
    read(options_id,*) nPoverS
    if (nangles>0 .and. nphi>0 .and. nPoverS>0) then
        exit
    else
        print*,'Error: There must be at least one each of ramp angles, phi values, and P/S values'
    end if
end do
nsegs = nangles+nPoverS+nphi-2
nbends = nsegs-1
do
    print*,'Choose type of backlimb deformation:'
    print*,'(1) Fault-Parallel Flow (Fold axes bisect fault bends)'
    print*,'(2) Fault Bend Folding (Preserves line length)'
    print*,'(3) Inclined Shear'
    read(options_id,*) backlimb_fold_type
    if (backlimb_fold_type==1) then
        exit
    else if(backlimb_fold_type==2) then
        print*,'Warning: Restored state beds are assumed to be flat for calculating fold axis orientation'
        exit
    else if(backlimb_fold_type==3) then
        !print*,'Enter shear angle in degrees (range: [0,90], 90 deg = vertical):'
        print*,'Shear angle should be in degrees (range: [0,180], <90 = antithetic, 90 = vertical, >90 = synthetic):'
        !print*,'(Only antithetic shear is allowed)'
        !read(options_id,*) shear_angle
        nparams = nparams+1
        !I think it makes more sense for this to be a parameter. It's hard to have an a priori intuition for it.
        exit
    else
        print*,'Invalid Command'
    end if
end do
do
    print*,'Must fault tip stay above a specified fault bend at all times? (0=no, 1=yes)'
    read(options_id,*) tipy_above_bend
    if (tipy_above_bend == 0) then
        exit
    else if (tipy_above_bend == 1) then
        do
            print*,'Enter bend number that tip must stay above: (1 is the highest bend.)'
            read(options_id,*) bend_tip_above
            if (bend_tip_above >=1 .and. bend_tip_above <= nangles-1) then
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        exit
    else
        print*,'Invalid Command'
    end if
end do
if (TipToSolve==1 .and. FaultType==6) then
    !Only do this for pure multibend fault, not for elliptic or spline that calls this since it would mess things up with those.
    !For elliptic faults, there is an option to have tipy at the detachment, which is similar to this.
    do
        print*,'Must initial fault tip position be at one of the bends? (0/1)'
        read(options_id,*) bend_at_tipy
        if (bend_at_tipy == 0) then
            nparams = nparams+nbends+nangles+nphi+nPoverS
            exit
        else if (bend_at_tipy == 1) then
            nparams = nparams+nbends+nangles+nphi+nPoverS-1
            do
                print*,'Enter bend number that tip must start at: (1 is the highest bend.)'
                read(options_id,*) bend_start
                if (bend_start >= 1 .and. bend_start <= nangles-1) then
                    exit
                else
                    print*,'Invalid Command'
                end if
            end do
            exit
        else
            print*,'Invalid Command'
        end if
    end do
else
    nparams = nparams+nbends+nangles+nphi+nPoverS
end if
do
    print*,'Must phi be less than the ramp angle in any fault segment the tip passes through? (0=no, 1=yes)'
    print*,'(This prevents lowering of the footwall due to a trishear zone extending downwards)'
    read(options_id,*) phi_le_ramp
    if (phi_le_ramp == 0 .or. phi_le_ramp == 1) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
do
    print*,'Can the trishear zone boundary intersect the backlimb fold axes above some elevation? (0 = no, 1 = yes)'
    if (options_id == 5) then !Manual input
        read(options_id,*) let_trishear_intersect_axes
    else !Reading from a file
        read(options_id,*) str
        if (str == 'let_trishear_intersect_axes') then
            let_trishear_intersect_axes = 1
        else
            backspace(options_id)
            let_trishear_intersect_axes = 0
        end if
    end if
    if (let_trishear_intersect_axes == 0) then
        exit
    else if(let_trishear_intersect_axes == 1) then
        print*,'Enter elevation above which trishear zone and fold axes can intersect: (must be above all data)'
        read(options_id,*) min_intersect_elevation
        exit
    else
        print*,'Invalid Command'
    end if
end do
do
    print*,'Must the final fault tip position be in the uppermost fault segment? (0 = no, 1= yes)'
    print*,'(If not, some fault segments may go unused.)'
    if (options_id == 5) then !Manual input
        read(options_id,*) final_tip_in_upper_segment
    else !Reading from a file.
        read(options_id,*) str
        if (str == 'final_tip_in_upper_segment') then
            final_tip_in_upper_segment = 1
        else
            final_tip_in_upper_segment = 0
            backspace(options_id)
        end if
    end if
    if (final_tip_in_upper_segment==0 .or. final_tip_in_upper_segment==1) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
do
    print*,'Do you want to limit the minimum and maximum x values of fault bends? (0/1)'
    if (options_id == 5) then !Manual input
        read(options_id,*) limit_fault_x_pts
    else !Reading from a file.
        read(options_id,*) str
        if (str == 'limit_fault_x') then
            limit_fault_x_pts = 1
        else
            limit_fault_x_pts = 0
            backspace(options_id)
        end if
    end if
    if (limit_fault_x_pts==1) then
        print*,'Enter minimum and maximum allowed x values for fault bends.'
        read(options_id,*) fault_x_limits
        exit
    else if (limit_fault_x_pts==0) then
        exit
    else
        print*,'Error: Invalid Response'
    end if
end do
do
    print*,'Must the fault be concave upward at all points? (0 = no, 1= yes (default))'
    if (options_id == 5) then !Manual input
        read(options_id,*) fault_concave_up
    else !Reading from a file.
        read(options_id,*) str
        if (str == 'fault_concave_up') then
            fault_concave_up = 1
        else if (str == 'fault_concave_up_not_required') then
            fault_concave_up = 0
        else
            fault_concave_up = 1
            backspace(options_id)
        end if
    end if
    if (fault_concave_up==0 .or. fault_concave_up==1) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
do 
    print*,'Do you want to put constraints on the x position of the backlimb syncline? (0/1)'
    if (options_id==5) then !Manual inptut
        read(options_id,*) constrain_backlimb_syncline
    else !Reading from a file
        read(options_id,*) str
        if (str == 'constrain_backlimb_syncline') then
            constrain_backlimb_syncline = 1
        else
            constrain_backlimb_syncline = 0
            backspace(options_id)
        end if
    end if
    if (constrain_backlimb_syncline == 1) then
        print*,'Enter elevation (y coordinate) at which to measure the backlimb syncline position.'
        read(options_id,*) y_backlimb_syncline
        print*,'Enter minimum and maximum allowed x position of backlimb syncline at y.'
        read(options_id,*) x_backlimb_syncline_limits
        exit
    else if (constrain_backlimb_syncline == 0) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
end subroutine trishear_multibend_options

subroutine run_model_multibend(output,params) !Subroutine to run a model and return the RMS error
use constants
use parameters, only: increment,v0,slip_sense
use options, only: ResultType,DatTypes,DatType_beds,DatType_dips,DatType_terr,DatType_faultpts,ndatatypes_poss,&
    DatType_beds_restored,beds_age_order,limit_bed_dips,min_bed_slopes,max_bed_slopes,FitType
use data_module
use err_and_prob
use data_uncertainties, only: GetSigmaBed,GetSigmaDip,GetSigmaTerr,GetSigmaFault,GetSigmaBedRestored,GetLc,GetSigma2Bed
use beds_order, only: check_beds_order
implicit none
double precision,intent(out) :: output !the output result (RMS or probability)
double precision,dimension(ndatatypes_poss) :: outputs !A single vector holding outputs for beds, dips, and terraces (in that order).
double precision,dimension(1:),intent(inout) :: params !values of the parameters
type(multibend_fault) :: fault !The fault.
double precision :: growth_slip_rem,terr_slip_rem !Remaining growth bed or terrace slip as we go down the segments.
double precision,dimension(nsegs) :: slipseg_bed,slipseg_terr !The amoung of slip with the tip in each segment. (Absolute value.) In the order that the tip moves through them.
integer :: tip_seg_init_bed,tip_seg_init_terr !tells which segment of the fault the fault tip is in. 1=uppermost segment
logical :: goodrun !tells if the run is good
double precision :: NaN = 0 !For creating NaN by dividing by zero
integer :: bedlength !number of points in the bed being worked on
double precision,dimension(:,:),allocatable :: thisbed !the bed being worked on
double precision :: slope,b !Slope and intercept of a bed
double precision :: slope_last, b_last !Slope and intercept of the previous bed restored.
double precision,dimension(2) :: bed_xlims,bed_xlims_last !x limits of current and previous bed restored.
double precision,dimension(:,:),allocatable :: dippos_temp !dip positions temporary
double precision,dimension(ndips) :: dips_temp !dips temporary
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
integer :: i,n,g ! counters, n is for beds, g is for growth beds
NaN = NaN/NaN

do !This do loop lets us exit if it turns out we have an impossible combination of parameters
    call interpret_multibend_fault_params(params,fault,goodrun)
    if (goodrun .eqv. .false.) then
        exit
    end if
    !Initialize the outputs array.
    outputs = 0 !For RMS or chi squared it should be 0 anyway. For probabilities, p should be 1, so ln(p) = 0.
    !Now restore data.
    !print*,params
    if (DatTypes(DatType_beds) == 1) then !Beds
        g = 1 !A counter of growth strata
        call GetSigmaBed(params,sigma_bed_x,sigma_bed_y) !Get the uncertainties in bed points.
        call GetLc(params,Lc) !Get the correlation length.
        if (ResultType==5) call GetSigma2Bed(params,sigma2_bed_x,sigma2_bed_y) !Get the uncertainties in bed points.
        if(DatTypes(DatType_beds_restored) == 1) then
            call GetSigmaBedRestored(params,sigma_bed_restored_x,sigma_bed_restored_y) !Get the uncertainty in restored bed points.
        end if
        do n = 1,nbeds !Go through each of the beds
            !print*,'bed =',n
            bedlength = beds(n)%npts
            allocate(thisbed(2,bedlength))
            allocate(Cbed(2,2,bedlength))
            Cbed = 0 !Initialize covariance matrix for bed points (x,y).
            Cbed(1,1,:) = sigma_bed_x(n)**2
            Cbed(2,2,:) = sigma_bed_y(n)**2
            thisbed(:,:) = beds(n)%pts
            if (beds(n)%growth .eqv. .true.) then !For growth strata
                !This should work for now. There may be a simpler way to do this (and the similar code for terraces), if we put all the calculating slipseg stuff back in the trishear functions
                !and set it up so tip final is the only tip position we need to pass to those functions. Then we can just pass them the slip, which in this case would be the growth_slip(g).
                growth_slip_rem = fault%growth_slip(g)
                slipseg_bed = 0
                do i = fault%tip_seg_final,fault%tip_seg_init
                    if (abs(growth_slip_rem)>abs(fault%slipseg(i))) then
                        slipseg_bed(i) = fault%slipseg(i)
                        growth_slip_rem = growth_slip_rem-fault%slipseg(i)
                    else
                        slipseg_bed(i) = growth_slip_rem
                        tip_seg_init_bed = i
                        exit
                    end if
                end do
                g = g+1
            else
                slipseg_bed = fault%slipseg
                tip_seg_init_bed = fault%tip_seg_init
            end if
            call trishear_func_multi_bend(thisbed,(/fault%tipxfinal,fault%tipyfinal/),fault%tip_seg_final,tip_seg_init_bed,v0,&
                fault%phi,fault%m,increment,-slip_sense*fault%PoverS,fault%s,fault%bendxy,fault%ramp_angle,fault%ramp_dir,&
                fault%axis_dips,fault%R,slipseg_bed,Cbed) !Run the trishear on this bed
            if (ResultType==5) then !Correlated and uncorrelated probabilities.
                allocate(Cbed2(2,2,bedlength))
                Cbed2 = 0 !Initialize covariance matrix for bed points (x,y).
                Cbed2(1,1,:) = sigma2_bed_x(n)**2
                Cbed2(2,2,:) = sigma2_bed_y(n)**2
            end if
            call calc_bed_err(n,bedlength,thisbed,beds_output_indiv,slope,b,Cbed,Lc,fault%restored_fit_params,Cbed2)
            if (DatTypes(DatType_beds_restored) == 1) then !Points on the restored beds.
                do i = 1,n_restored_beds
                    if (beds(n)%ident == restored_beds(i)%ident) then
                        call calc_restored_bed_err(i,slope,b,restored_beds_output_indiv(i),sigma_bed_restored_x,&
                        sigma_bed_restored_y,fault%restored_fit_params)
                    end if
                end do
            end if
            if (beds_age_order==1 .or. beds_age_order == 2) then
                bed_xlims = (/minval(thisbed(1,:)),maxval(thisbed(1,:))/)
                if (check_beds_order(n,b,slope,b_last,slope_last,bed_xlims,bed_xlims_last) .eqv. .false.) then
                    goodrun = .false.
                end if
                b_last = b !For the next bed, this will be from the last bed for comparison.
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
        sigma_dip_array = sigma_dip*pi/180.
        dippos_temp = dip_pos !Using the originals would mess up the data for other runs.
        dips_temp = dips
        call trishear_func_dip_multi_bend(dippos_temp,dips_temp,(/fault%tipxfinal,fault%tipyfinal/),fault%tip_seg_final,&
            fault%tip_seg_init,v0,fault%phi,fault%m,increment,-slip_sense*fault%PoverS,fault%s,fault%bendxy,fault%ramp_angle,&
            fault%ramp_dir,fault%axis_dips,fault%R,fault%slipseg,sigma_dip_array)
        RestoredDips = dips_temp*180/pi !Convert back to degrees
        !Fix dips > 90 deg or < -90 deg
        where(RestoredDips > 90.) RestoredDips = RestoredDips - 180. !Fix dips > 90 deg
        where(RestoredDips < -90.) RestoredDips = RestoredDips + 180. !Fix dips < - 90 deg
        call calc_dip_err(RestoredDips,outputs(2),sigma_dip_array*180./pi,sigma_restored_dip,fault%restored_fit_params)
    end if
    if (DatTypes(DatType_terr) == 1) then !Terraces
        call GetSigmaTerr(params,sigma_terr,sigma_orig_elev)
        do n = 1,nterr !Go through each terrace
            !Now restore the terrace
            bedlength = terraces(n)%npts !I'm keeping the "bed" variable names since there's no point creating new ones. They do the same thing either way.
            allocate(thisbed(2,bedlength))
            allocate(Cterr(2,2,bedlength))
            Cterr = 0 !Initialize covariance matrix for bed points (x,y).
            Cterr(1,1,:) = sigma_terr(n)**2
            Cterr(2,2,:) = sigma_terr(n)**2
            thisbed(:,:) = terraces(n)%pts
            terr_slip_rem = fault%terr_slip(n)
            slipseg_terr = 0
            do i = fault%tip_seg_final,fault%tip_seg_init
                if (abs(terr_slip_rem)>abs(fault%slipseg(i))) then
                    slipseg_terr(i) = fault%slipseg(i)
                    terr_slip_rem = terr_slip_rem-fault%slipseg(i)
                else
                    slipseg_terr(i) = terr_slip_rem
                    tip_seg_init_terr = i
                    exit
                end if
            end do
            call trishear_func_multi_bend(thisbed,(/fault%tipxfinal,fault%tipyfinal/),fault%tip_seg_final,tip_seg_init_terr,v0,&
            fault%phi,fault%m,increment,-slip_sense*fault%PoverS,fault%s,fault%bendxy,fault%ramp_angle,fault%ramp_dir,&
            fault%axis_dips,fault%R,slipseg_terr,Cterr) !Run the trishear on this bed
            call calc_terr_err(n,bedlength,thisbed,terr_output_indiv,Cterr,sigma_orig_elev)
            deallocate(thisbed,Cterr)
        end do
        call calc_beds_err(terr_output_indiv,outputs(3))
    end if
    if (DatTypes(DatType_faultpts) == 1) then !Fault Points
        call GetSigmaFault(params,sigma_fault)
        do n=1,nfaultpts
            do i = 1,nbends
            if (faultpts(2,n)>(faultpts(1,n)-fault%bendxy(1,i))*tan((pi*fault%ramp_dir+fault%ramp_angle(i)+&
                fault%ramp_angle(i+1))/2)+fault%bendxy(2,i)) then !Closest to ramp i.
                fault_errs(n) = abs(((faultpts(1,n)-fault%bendxy(1,i))*tan(fault%ramp_angle(i))-faultpts(2,n)+&
                    fault%bendxy(2,i))/sqrt(1+tan(fault%ramp_angle(i))**2))
                    exit
            else if(i == nbends) then !Closest to lower ramp
                fault_errs(n) = abs(((faultpts(1,n)-fault%bendxy(1,nbends))*tan(fault%ramp_angle(nsegs))-faultpts(2,n)&
                    +fault%bendxy(2,nbends))/sqrt(1+tan(fault%ramp_angle(nsegs))**2))
            end if
            end do
            if (faultpts(2,n)>fault%tipyfinal) then
                if (sqrt((faultpts(2,n)-fault%tipyfinal)**2+(faultpts(1,n)-fault%tipxfinal)**2)>fault_errs(n)) then
                    fault_errs(n) = sqrt((faultpts(2,n)-fault%tipyfinal)**2+(faultpts(1,n)-fault%tipxfinal)**2)
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
end subroutine run_model_multibend

subroutine interpret_multibend_fault_params(params,fault,goodrun)
use data_module, only: ngrowth,nterr
use options, only: DatTypes,DatType_beds,DatType_terr,FitType,terr_age_order
use fault_tip, only: find_tip_position,check_fault_tip
use parameters, only: nparam_start_growth,nparam_start_restored_fit,nparam_start_terr,increment,slip_sense,nrestored_fit_params
use constants
use math, only: maxind,sec,swap
use growth_strata, only: GetGrowthSlip
implicit none
double precision,dimension(1:),intent(in) :: params !Model parameters to be read and interpreted.
type(multibend_fault),intent(out) :: fault !The fault.
logical,intent(out) :: goodrun !tells if the parameters are acceptable or not.
double precision :: tipx,tipy,total_slip,total_slip_mag,s !parameters
double precision,dimension(nsegs) :: phi,phi_deg,m,PoverS !These values are now allowed to change.
double precision,dimension(nsegs) :: ramp_angle, ramp_angle_deg !1st value is highest segment.
double precision :: shear_angle !For inclined shear. 0 deg is horizontal. 90 deg is vertical. Allowed range is 0 deg to 180 deg.
double precision,dimension(nbends) :: axis_dips !Fold axis dips. (90 deg / pi/2 rad is vertical, 0 is horizontal.)
double precision :: axis_dip !The dip of the current axis being dealt with in a loop that iterates through the axes / fault bends.
double precision,dimension(2,nbends) :: bendxy !x,y coordinates of the fault bends. Numbered from highest = 1 to lowest = nbends.
double precision :: y_tip_above !If tipy_above_bend == 1, then this number tells the elevation of the bend that the tip has to stay above.
double precision,dimension(nsegs) :: R !Ratio of slip on a given segment relative to slip on the lowest segment.
double precision,dimension(nsegs) :: slipseg !The amoung of slip with the tip in each segment. (Absolute value.) In the order that the tip moves through them.
double precision :: tipxinit, tipyinit, tipxfinal, tipyfinal !Initial and final (in terms of forward motion) tip positions
integer :: tip_seg,tip_seg_init,tip_seg_final !tells which segment of the fault the fault tip is in. 1=uppermost segment
integer :: ramp_dir !dip direction of ramp
double precision :: slope_tri,slope_axis !Slopes of the trishear zone boundary and the fold axis, used for checking whether and at what elevation they intersect.
double precision :: x_backlimb_syncline !x position of backlimb syncline.
double precision,dimension(nsegs) :: ramp_acute !Acute versions of ramp angles
double precision :: phi_fault_bend !Difference between the two ramp angles, for fault bend folding.
double precision :: f !Result of function evaluation, used in Newton-Raphson parts. For fault-bend folding.
double precision :: prec_limit !Precision limit that f must be lower than for the Newton-Raphson evaluations to stop. For fault-bend folding.
integer :: niter,max_iter !Number of iterations performed and maximum number to use in Newton-Raphson method. For fault-bend folding.
double precision :: rn !Holds a random number. For fault-bend folding.
double precision :: gamma1 !Backlimb axial angle for fault-bend folding.
integer :: i ! a counter
integer :: n
double precision :: intersect_elevation !Elevation at which the trishear zone and a fold axis intersect.
double precision,dimension(nrestored_fit_params) :: restored_fit_params
double precision,dimension(ngrowth) :: growth_slip !Slip needed to restore each growth strata bed.
double precision,dimension(nterr) :: terr_slip !Slip needed to restore each terrace.
goodrun = .true. !Initialize it to true. So far, this is a good run.
do !This do loop lets us quickly exit anytime it turns out we have an impossible combination of parameters
    !Make assignments
    tipx = params(1)
    tipy = params(2)
    total_slip_mag = params(3)
    do i = 1,nangles-1
        ramp_angle_deg(i) = params(3+i) !First listed is highest ramp segment, last is lowest.
        phi_deg(i) = 999 !A place holder value to mean we need to fill this in.
        PoverS(i) = 999 !A place holder value to mean we need to fill this in.
        if (bend_at_tipy == 0) then
            bendxy(2,i) = params(3+nangles+i) !Again, numbered from top to bottom.
        else !Tip has to start at a bend.
            if (i < bend_start) then
                bendxy(2,i) = params(3+nangles+i)
            else if (i == bend_start) then
                bendxy(2,i) = tipy !Since TipToSolve must == 1 in this case.
            else
                bendxy(2,i) = params(2+nangles+i)
            end if
            !An alternative to all this would be to think of tipy as just one more bend point (the top one),
            !And if bend_at_tipy == 1, use it to calculate total slip or something like that.
            !One could even search just for fault points (x,y) and not ramp angles at all, though there are other problems there.
            !Another option is to replace forcing the tip to start at a bend with requiring it to start near a bend, with some uncertainty, making this be another model constraint.
        end if
        if (tipy_above_bend==1 .and. i == bend_tip_above) then !Get the depth that the tip has to stay above
            y_tip_above = bendxy(2,i)
        end if
        if (i>1 .and. bendxy(2,i)>bendxy(2,i-1)) then !Check that the bends get progressively lower in elevation.
            goodrun = .false.
        end if
        if (fault_concave_up==1 .and. i>1 .and. min(ramp_angle_deg(i),180.-ramp_angle_deg(i)) > &
                min(ramp_angle_deg(i-1),180.-ramp_angle_deg(i-1))) then !The fault must get steeper, not shallower, as it gets higher (lower i). Note this may change later.
            goodrun = .false.
        end if
    end do
    ramp_angle_deg(nsegs) = params(3+nangles)
    if (fault_concave_up==1 .and. nangles>1 .and. min(ramp_angle_deg(nsegs),180-ramp_angle_deg(nsegs)) > &
            min(ramp_angle_deg(nangles-1),180-ramp_angle_deg(nangles-1))) then !The fault must get steeper, not shallower, as it gets higher
        goodrun = .false.
    end if
    do i = 1,nphi-1 !Get phi values and y values at which phi changes
        phi_deg(i+nangles-1) = params(2+2*nangles+i-bend_at_tipy)
        ramp_angle_deg(i+nangles-1) = 999 !A place holder value to mean we need to fill this in.
        PoverS(i+nangles-1) = 999 !A place holder value to mean we need to fill this in.
        bendxy(2,i+nangles-1) = params(2+2*nangles+i+nphi-bend_at_tipy)
        if (i>1 .and. bendxy(2,i+nangles-1)>bendxy(2,i+nangles-2)) then !Check that the bends get progressively lower in elevation.
            goodrun = .false.
        end if
    end do
    phi_deg(nsegs) = params(2+2*nangles+nphi-bend_at_tipy) !The lowest one always goes to the last segment.
    do i = 1,nPoverS-1 !Get P/S values and y values at which P/S changes.
        PoverS(i+nangles+nphi-2) = params(1+2*nangles+2*nphi+i-bend_at_tipy)
        ramp_angle_deg(i+nangles+nphi-2) = 999 !A place holder value to mean we need to fill this in.
        phi_deg(i+nangles+nphi-2) = 999 !A place holder value to mean we need to fill this in.
        bendxy(2,i+nangles+nphi-2) = params(1+2*nangles+2*nphi+i+nPoverS-bend_at_tipy)
        if (i>1 .and. bendxy(2,i+nangles+nphi-2)>bendxy(2,i+nangles+nphi-3)) then !Check that the bends get progressively lower in elevation.
            goodrun = .false.
        end if
    end do
    PoverS(nsegs) = params(1+2*nangles+2*nphi+nPoverS-bend_at_tipy) !The lowest one always goes to the last segment.
    !Sort the ramp angle, phi, and P/S changes. This is a simple selection sort. I could probably use a better (faster) sorting algorithm.
    do i = 1,nbends
        n = maxind(bendxy(2,i:nbends))+i-1
        call swap(bendxy(2,i:i),bendxy(2,n:n))
        call swap(ramp_angle_deg(i:i),ramp_angle_deg(n:n))
        call swap(phi_deg(i:i),phi_deg(n:n))
        call swap(PoverS(i:i),PoverS(n:n))
    end do
    !Fill in values that don't change at each bend.
    do i = 1,nsegs
        if (ramp_angle_deg(i) == 999) then
            do n = i+1,nsegs
                if (ramp_angle_deg(n) /= 999) then
                    ramp_angle_deg(i) = ramp_angle_deg(n)
                    exit
                end if
            end do
        end if
        if (phi_deg(i) == 999) then
            do n = i+1,nsegs
                if (phi_deg(n) /= 999) then
                    phi_deg(i) = phi_deg(n)
                    exit
                end if
            end do
        end if
        if (PoverS(i) == 999) then
            do n = i+1,nsegs
                if (PoverS(n) /= 999) then
                    PoverS(i) = PoverS(n)
                    exit
                end if
            end do
        end if
    end do
    s = params(4+nangles+nphi+nPoverS+nbends-bend_at_tipy)
    if (backlimb_fold_type==3) then !Inclined shear
        shear_angle = params(5+nangles+nphi+nPoverS+nbends-bend_at_tipy) !For now, at least, I'm making this be the last of the fault parameters but before the data specific params.
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
    ramp_angle = ramp_angle_deg*pi/180; !convert to radians
    !The fault should not change dip direction. Note: for a flat segment, use 180 deg rather than 0 deg for dipping right.
    if (minval(ramp_angle)<=pi/2 .and. maxval(ramp_angle)<=pi/2) then !1=dips left, -1 = dips right
        ramp_dir = 1
    else if(minval(ramp_angle)>pi/2 .and. maxval(ramp_angle)>pi/2) then
        ramp_dir = -1
    else !This should only occur if ramp_dir is different for the two ramp angles.
        goodrun = .false.
        exit
    end if
    if (goodrun .eqv. .false.) then
        exit
    end if
    phi = phi_deg*pi/180 !convert to radians
    !Calculate fold axes.
    if (backlimb_fold_type==1) then !Fault parallel flow
        do i = 1,nbends
            axis_dips(i) = (pi*ramp_dir+ramp_angle(i)+ramp_angle(i+1))/2 !Measured from horizontal. Currently doing full [0,180], rather than using ramp_acute. This works if fault is concave or convex.
        end do
    else if(backlimb_fold_type==2) then !Fault bend folding
        !Note: This assumes that the beds are initially flat.
        !Also, this could be moved out into a function, which could be called by both this module and the Parallel_FPF one.
        do i = 1,nsegs
            ramp_acute(i) = min(ramp_angle(i),pi-ramp_angle(i)) !Acute version of ramp angle
        end do
        do i = 1,nbends
            phi_fault_bend = ramp_acute(i)-ramp_acute(i+1) !Phi should be acute as well.
            !Calculate gamma angles. I can't figure out how to solve the gamma1 or gamma_star equations analytically. (No answer from Mathematica.) So I'll have to use Newton-Raphson.
            prec_limit = 1e-5 !Necessary precision at which we stop trying to improve the value of gamma1 or gamma star. Possibly allow the user to choose, or base it on number of degrees. e.g. difference between tan(phi) and tan(phi+1deg).
            max_iter = 1e3 !Maximum number of iterations after which we stop searching with Newton_Raphson.
            !Solve for gamma1 with Suppe and Medwedeff eq. 9
            !gamma1 should be >pi/4, or the fold will have an overturned backlimb
            !call init_random_seed()
            niter = 0 !Number of iterations
            call random_number(rn)
            gamma1 = rn*pi/4+pi/4 !Initial random gamma1 in the range [pi/4,pi/2). This should avoid any values on the wrong side of the asymptote.
            f = (-sin(gamma1 + ramp_acute(i+1))*(sin(2*gamma1 + ramp_acute(i+1)) + sin(ramp_acute(i+1))))/ &
                (cos(gamma1 + ramp_acute(i+1))*(sin(2*gamma1 + ramp_acute(i+1)) + sin(ramp_acute(2))) - sin(gamma1)) &
                - tan(phi_fault_bend)
            do while (abs(f) > prec_limit .or. gamma1<0 .or. gamma1>pi/2) !Newton-Raphson
                gamma1 = gamma1 - ((cos(gamma1+phi_fault_bend)-cos(gamma1)*cos(2*gamma1+2*ramp_acute(i+1)+phi_fault_bend))* &
                    sec(phi_fault_bend)*(-sin(gamma1)+cos(gamma1)*sin(2*(gamma1+ramp_acute(i+1)))))/(cos(2*gamma1) &
                    -cos(2*ramp_acute(i+1))+4*sin(gamma1+ramp_acute(i+1))**2)
                if (gamma1 > pi/2) then !Keep gamma1 in the first quadrant.
                    gamma1 = pi-gamma1
                else if (gamma1 < 0) then
                    gamma1 = -gamma1
                end if
                f = (-sin(gamma1 + ramp_acute(i+1))*(sin(2*gamma1 + ramp_acute(i+1)) + sin(ramp_acute(i+1))))/ &
                    (cos(gamma1 + ramp_acute(i+1))*(sin(2*gamma1 + ramp_acute(i+1)) + sin(ramp_acute(i+1))) - sin(gamma1)) &
                    - tan(phi_fault_bend)
                niter = niter+1
                if (niter == max_iter) then
                    print*,'Error: No solution for gamma1 found within desired precision.'
                    goodrun = .false.
                    exit
                end if
            end do
            axis_dips(i) = ((ramp_dir+1)/2)*pi-ramp_dir*gamma1 !This could be changed to just ramp_dir*gamma1, if I wanted to use the range [-pi/2,pi/2] instead of [0,pi]. Same goes for inclined shear below, and shear_angle.
        end do
    else if(backlimb_fold_type==3) then !Inclined shear
        do i = 1,nbends
            axis_dips(i) = ((ramp_dir+1)/2)*pi-ramp_dir*shear_angle*pi/180. !Will be in the range [90,180] if ramp_dir==1, and [0,90] if ramp_dir==-1
        end do
        if (ramp_dir*axis_dips(1)<maxval(ramp_dir*ramp_angle)) then !This checks that synthetic shear doesn't result in shear planes dipping less steeply than the fault itself.
            goodrun = .false.
            exit
        end if
    else
        print*,'Invalid Command'
    end if
    m = tan(phi)
    do i = 1,nsegs
        if (phi(i) == pi/2) then
                m = tan(89*pi/180) !Since tan(pi) = inf, which causes problems
        end if
    end do
    !Determine which segment the tip (tipx,tipy) is in. Note: this can be initial or final tip depending which one tipx,tipy represent.
    tip_seg = 1
    do i = 1,nbends
        if (tipy <= bendxy(2,i)) then
            tip_seg = tip_seg+1
        else
            exit
        end if
    end do
    !Calculate the bend x coordinates
    if (bend_at_tipy == 0) then
        if (tip_seg > 1) then
            bendxy(1,tip_seg-1) = tipx-(tipy-bendxy(2,tip_seg-1))/tan(ramp_angle(tip_seg))
            do i = tip_seg-2,1,-1
                bendxy(1,i) = bendxy(1,i+1)-(bendxy(2,i+1)-bendxy(2,i))/tan(ramp_angle(i+1))
            end do
        end if
        if (tip_seg < nsegs) then
            bendxy(1,tip_seg) = tipx-(tipy-bendxy(2,tip_seg))/tan(ramp_angle(tip_seg))
            do i = tip_seg+1,nbends
                bendxy(1,i) = bendxy(1,i-1)-(bendxy(2,i-1)-bendxy(2,i))/tan(ramp_angle(i))
            end do
        end if
    else
        bendxy(1,bend_start) = tipx !Again, since TipToSolve==1 is required in order to have bend_at_tipy==1
        do i = bend_start-1,1,-1
            bendxy(1,i) = bendxy(1,i+1)-(bendxy(2,i+1)-bendxy(2,i))/tan(ramp_angle(i+1))
        end do
        do i = bend_start+1,nbends
            bendxy(1,i) = bendxy(1,i-1)-(bendxy(2,i-1)-bendxy(2,i))/tan(ramp_angle(i))
        end do
    end if
    if (limit_fault_x_pts==1) then
        if (minval(bendxy(1,:))<fault_x_limits(1) .or. maxval(bendxy(1,:))>fault_x_limits(2)) then
            goodrun = .false.
            exit
        end if
    end if
    !Calculate the ratio of slip on each segment relative to slip on the lowest segment.
    R(nsegs) = 1.
    do i = nsegs-1,1,-1
        R(i) = R(i+1)*(sin(ramp_angle(i+1))-cos(ramp_angle(i+1))*tan(axis_dips(i)))/ &
            (sin(ramp_angle(i))-cos(ramp_angle(i))*tan(axis_dips(i))) !The slip ratio for each segment compared to the lowest segment.
    end do
    !Calculate initial and final tipx and tipy (final = present-day, initial = when fault motion began). Final must be higher than initial.'
    call find_tip_position(tipx,tipy,tipxinit,tipyinit,tipxfinal,tipyfinal,tip_seg,tip_seg_init,tip_seg_final,&
        nsegs,slipseg,bendxy,total_slip,PoverS,R,ramp_angle)
    !Check that the fault tip obeys any additional restrictions on initial or final position.
    if (.not. check_fault_tip(tipxinit,tipyinit,tipxfinal,tipyfinal)) then
        goodrun = .false.
        exit
    end if
    !Check that the trishear boundary won't cross the kink axes.
    axis_dip = 999 !A placeholder value to mean that there is so far no fold axis below the point (if we're in a low segment).
    do n = nsegs,tip_seg_final,-1
        if (axis_dip /= 999 .and. n <= tip_seg_init) then !If it == 999, then we're in the lowest segment (as defined by changes in dip), so there are no axes below it to worry about, and if it's below the initial tip segment, the tip will never be there.
            if (let_trishear_intersect_axes==0) then !Any intersection between the trishear zone and the fold axis below it is a problem.
                if (ramp_dir*(ramp_angle(n)+ramp_dir*phi(n))>=ramp_dir*axis_dip) then
                    goodrun = .false.
                    exit
                end if
            else !The trishear zone can intersect the synclinal axes above some elevation
                if (n == tip_seg_init) then !In the initial tip segment.
                    slope_tri = tan(ramp_angle(tip_seg_init)+ramp_dir*phi(tip_seg_init))
                    if (axis_dip==pi/2) then !A vertical axis (such as in the vertical shear case) is a special case, since tan(pi/2) = infinity.
                        intersect_elevation = (bendxy(1,tip_seg_init)-tipxinit)*slope_tri+tipyinit
                    else
                        slope_axis = tan(axis_dip)
                        intersect_elevation = (((tipxinit-bendxy(1,tip_seg_init))*slope_tri*slope_axis-tipyinit*slope_axis&
                            +bendxy(2,tip_seg_init)*slope_tri)/(slope_tri-slope_axis))
                    end if
                    if (intersect_elevation<min_intersect_elevation .and. intersect_elevation>bendxy(1,tip_seg_init)) then !Above the min elevation is allowed, and below the bend means it's dipping the other way and is harmless.
                        goodrun = .false.
                        exit
                    end if
                else if (n == tip_seg_final) then
                    slope_tri = tan(ramp_angle(tip_seg_final)+ramp_dir*phi(tip_seg_final))
                    if (axis_dip==pi/2) then !A vertical axis (such as in the vertical shear case) is a special case, since tan(pi/2) = infinity.
                        intersect_elevation = (bendxy(1,tip_seg_final)-tipxfinal)*slope_tri+tipyfinal
                    else
                        slope_axis = tan(axis_dip) !Note: For now, I'm only comparing to the first axis below the fault tip. For inclined simple shear, or for two fault segments, this is sufficient. If there are multiple fault segments, and the axes intersect, this may not be, but that isn't yet dealt with below either.
                        intersect_elevation = (((tipxfinal-bendxy(1,tip_seg_final))*slope_tri*slope_axis-tipyfinal*slope_axis &
                                +bendxy(2,tip_seg_final)*slope_tri)/(slope_tri-slope_axis))
                    end if
                    if (intersect_elevation<min_intersect_elevation .and. intersect_elevation>bendxy(1,tip_seg_final)) then
                        goodrun = .false.
                        exit
                    end if
                end if
                if (n /= tip_seg_init .and. bendxy(2,n)<min_intersect_elevation .and. ramp_dir*(ramp_angle(n)+ramp_dir*phi(n)) & !This needs to be done for tip_seg_final too, since the tip will pass that bend, but not initial. (Note: for normal sense slip the opposite would be true.)
                        >=ramp_dir*axis_dip .and. ramp_angle(n)/=ramp_angle(n-1)) then !If the tip is at the bend, then the trishear zone must not go past the axis.
                    goodrun = .false.
                    exit
                end if
            end if
        end if
        if (n /= 1 .and. ramp_angle(n) /= ramp_angle(n-1)) then !Then this is a real bend in the fault coming up
            axis_dip = axis_dips(n-1) !New last axis dip behind the tip position. (i.e. axis above the current fault segment, but below the next one)
        end if
    end do
    !Check that the axis is within allowed x limits if that constraint is included.
    !Note if axes can intersect and refract this won't necessarily work, but that case isn't properly dealt with anyway.
    if (constrain_backlimb_syncline == 1) then
        do n = nbends,1,-1
            if (ramp_angle(n+1)/=ramp_angle(n)) then !We're looking for the farthest backlimb axis that is a real fault bend.
                x_backlimb_syncline = (y_backlimb_syncline-bendxy(2,n))/tan(axis_dips(n))+bendxy(1,n)
                if (x_backlimb_syncline <= x_backlimb_syncline_limits(1) .or. &
                    (x_backlimb_syncline >= x_backlimb_syncline_limits(2))) then
                    goodrun = .false.
                end if
                exit
            end if
        end do
    end if
    !Check that tipy_init > bend_y if required
    if (tipy_above_bend == 1 .and. tipyinit <= y_tip_above) then
        goodrun = .false.
        exit
    end if
    !Check that the final tip position is in the uppermost segment if required.
    if (final_tip_in_upper_segment == 1 .and. tip_seg_final /= 1) then
        goodrun = .false.
        exit
    end if
    !Check that phi > ramp angle if required. This is only for segments the tip spends time in.
    if (phi_le_ramp == 1) then
        do i = tip_seg_final,tip_seg_init
            if (phi(i) > min(ramp_angle(i),pi-ramp_angle(i))) then
                goodrun = .false.
                exit
            end if    
        end do
    end if
    fault%tipxinit = tipxinit
    fault%tipyinit = tipyinit
    fault%tipxfinal = tipxfinal
    fault%tipyfinal = tipyfinal
    fault%tip_seg_init = tip_seg_init
    fault%tip_seg_final = tip_seg_final
    fault%phi = phi
    fault%m = m
    fault%PoverS = PoverS
    fault%s = s
    fault%bendxy = bendxy
    fault%ramp_angle = ramp_angle
    fault%axis_dips = axis_dips
    fault%R = R
    fault%slipseg = slipseg
    fault%ramp_dir = ramp_dir
    fault%growth_slip = growth_slip
    fault%terr_slip = terr_slip
    fault%restored_fit_params = restored_fit_params
    exit
end do
end subroutine interpret_multibend_fault_params

subroutine trishear_func_multi_bend(pts,tip_start,tip_seg_start,tip_seg_end,v0,phi,m,increment,PoverS,s,&
    bendxy,ramp_angle,ramp_dir,axis_dips,R,slipseg,Cbed)
!Note: This one takes in points in the trishear coord system around the starting tip (usually tip_final) and returns them in coords around the ending tip.
use constants
use math, only: rot2D,xy_to_ze2D,ze_to_xy2D
use analytic_trishear_module, only: analytic_trishear
use options, only: propagate
use err_and_prob, only: prop_bend_error
implicit none
double precision, dimension(1:,1:), intent(inout) :: pts
double precision,dimension(2) :: tip_start !Starting tip position.
integer :: tip_seg_start,tip_seg_end !Starting and ending tip segments. An alternative would be to calculate tip_seg_end within this program.
double precision, dimension(2,nbends), intent(in) :: bendxy
double precision, intent(in) :: v0,increment,s
integer,intent(in) :: ramp_dir
double precision,dimension(nsegs),intent(in) :: R !Ratio of slip on a given segment relative to slip on the lowest segment.
double precision,dimension(nsegs),intent(in) :: slipseg !The amoung of slip with the tip in each segment. (Absolute value.) In the order that the tip moves through them.
double precision, dimension(nsegs),intent(in) :: ramp_angle,phi,m,PoverS
double precision,dimension(nbends),intent(in) :: axis_dips
double precision,dimension(1:,1:,1:),intent(inout) :: Cbed !Covariance matrices for x and y for all the points in the bed.
double precision,dimension(2) :: tip !The location of the fault tip. Moves as the fold is restored.
double precision,dimension(2) :: pt !The point being moved
double precision,dimension(2,1) :: ptze !The point in trishear coordinates
integer :: k,i
!double precision,dimension(nbends) :: axis_slopes !Slopes of fold axes.
double precision,dimension(nsegs) :: sin_ramp,cos_ramp !Just calculate these once for speed.
double precision :: half_v0,xexp,yexp,slip,x,y,signy
double precision,dimension(2,2) :: C,J !Covariance matrix for the specific point being worked on and Jacobian matrix for motion in the trishear velocity field.
double precision,dimension(2,2) :: Refl,Rot !Reflection and rotation matrixes for translations to or from the trishear coordinate system.
double precision :: rem_slip
integer :: slip_sign !1 for forwards, -1 for backwards.
double precision :: slip_tri,slip_next !Slip to get into a specific zone.
double precision :: vx, vy !Trishear velocities
double precision,dimension(2) :: v,vb,vnext !Velocity of a point and a boundary, and velocity in the next fault segment that the point is moving towards.
integer :: loc !tells which zone a point is in.
integer :: nseg_tip !number of the segments that the tip will cross through during the course of the model.
integer :: n !a counter
!Variables added or changed to accommodate the change to treat fault bends differently from changes in P/S or phi.
double precision,dimension(nangles) :: ramp_angle_loc,sin_ramp_loc,cos_ramp_loc,tan_ramp_loc,R_loc
double precision,dimension(nangles-1) :: axis_slopes !Slopes of fold axes.
double precision,dimension(2,nangles-1) :: bendxy_loc
integer,dimension(nsegs) :: tip_loc

!Separate changes in fault angle from other changes separating segments. This is added in quickly after I realized there is a problem with the way things are currently done.
!It would probably be better eventually to do this kind of thing up top or to more completely overhaul how ramp segments of different kinds are dealt with.
!The ending loc is because these go with the loc variable that denotes which fault ramp angle a point is above.
i = 1
ramp_angle_loc(1) = ramp_angle(1)
R_loc(1) = R(1)
tip_loc(1) = 1
do n = 2,nsegs
    if (ramp_angle(n) /= ramp_angle(n-1)) then
        ramp_angle_loc(i+1) = ramp_angle(n)
        R_loc(i+1) = R(n)
        bendxy_loc(:,i) = bendxy(:,n-1)
        axis_slopes(i) = tan(axis_dips(n-1))
        i = i+1
    end if
    tip_loc(n) = i !This lets us translate a tip segment to the corresponding value of the "loc" variable.
end do
sin_ramp_loc = sin(ramp_angle_loc)
cos_ramp_loc = cos(ramp_angle_loc)
tan_ramp_loc = tan(ramp_angle_loc)


!axis_slopes = tan(axis_dips)
sin_ramp = sin(ramp_angle)
cos_ramp = cos(ramp_angle)
!tan_ramp = tan(ramp_angle)
half_v0 = v0/2; !v0/2
xexp = 1/s; !exponent in vx term
yexp = (1+s)/s; !exponent in vy term
if (v0 >= 0) then !determine the sign of the slip
    slip_sign = 1
else
    slip_sign = -1
end if
nseg_tip = abs(tip_seg_end-tip_seg_start)+1 !Number of segments tip will pass through.
do n = tip_seg_start,tip_seg_end !Loop through each segment that the tip goes through. Note this is fine when moving down with higher segments having lower numbers, but it won't work going the other way.
    do k = 1,size(pts,2)
        !print*,'seg = ',n,'k = ',k
        pt = pts(:,k)
        if (n == tip_seg_start) then
            tip = tip_start
        else
            tip = bendxy(:,n+(slip_sign-1)/2)
        end if
        slip = 0 !slip so far
        loc = 0 !location: 0 = unknown, 1 to nangles = fault segments, nangles+1 = trishear, nangles+2 = footwall
        C = Cbed(:,:,k) !Covariance matrix for x and y for this point.
        rem_slip = slipseg(n)
        do while (rem_slip*slip_sign>0)
            !print*,loc
            if (loc == 0) then !need to find out location
                !I also need to deal with what happens if two fold axes intersect.
                if ((ramp_dir*(pt(1)-tip(1))>ramp_dir*(pt(2)-tip(2))/tan(ramp_angle(n)+ramp_dir*phi(n))) .and.&
                        ((pt(2)-tip(2))>(pt(1)-tip(1))*tan(ramp_angle(n)-ramp_dir*phi(n)))) then !Trishear zone
                    loc = nangles+1
                else
                    do i = tip_loc(n),(nangles-1) !Segment i is always the segment above bend i.
                        if (pt(1)*ramp_dir>bendxy_loc(1,i)*ramp_dir .and. pt(2)-bendxy_loc(2,i)<(pt(1)-bendxy_loc(1,i))&
                                *tan_ramp_loc(i)) then
                            loc = nangles+2 !Footwall, below segment i
                            exit
                        else if (ramp_dir*(pt(1)-bendxy_loc(1,i))>ramp_dir*(pt(2)-bendxy_loc(2,i))/axis_slopes(i)) then!In the hangingwall
                            loc = i !Above segment i
                            exit
                        end if
                    end do
                    if (loc ==0) then !If we still haven't assigned it to anywhere, the only remaining place is above or below the lowest ramp.
                        if ((nangles>1 .and. pt(2)-bendxy_loc(2,nangles-1)>(pt(1)-bendxy_loc(1,nangles-1))*tan_ramp_loc(nangles)) &
                            .or. (nangles==1 .and. pt(2)-tip(2)>(pt(1)-tip(1))*tan_ramp_loc(nangles))) then !If there's only one fault segment, then there are no bends, so no bendxy_loc
                            loc = nangles !Hanging wall, lowest segment.
                        else
                            loc = nangles+2 !Footwall, lowest segment.
                        end if
                    end if
                end if
            else if (loc == nangles+1) then !Trishear
                ptze(:,1) = pt
                ptze = xy_to_ze2D(ptze,tip(1),tip(2),min(ramp_angle(n),pi-ramp_angle(n)),ramp_dir); !Rotate into trishear coordinates.
                x = ptze(1,1)
                y = ptze(2,1)
                Refl = reshape((/ramp_dir,0,0,1/),(/2,2/))
                Rot = reshape((/ramp_dir*cos_ramp(n),sin_ramp(n),-sin_ramp(n),ramp_dir*cos_ramp(n)/),&
                    (/2,2/))
                !Note: If we bring in the points non-reflected, then we need to reflect here if necessary.
                C = matmul(matmul(Refl,C),Refl) !Reflect C. (This will only have an effect if ramp_dir == -1; otherwise Refl is an identity matrix.)
                C = matmul(matmul(transpose(Rot),C),Rot) !Change C into the rotated (zeta,eta) coordinate system. The rotation matrix is its own Jacobian.
                slip = 0
                if (s==1) then !Use the analytic solution for s==1
                    call analytic_trishear(x,y,C,m(n),PoverS(n),rem_slip*R(n),slip_sign,slip)
                    slip = slip/R(n) !Convert slip in the trishear segment to the slip on the lowest segment that we're using to keep track of the slip accrued.
                else
                    do while (slip == 0 .or. (y<=m(n)*x .and. y>=-m(n)*x))
                        if (x<0) then
                            print*,'Error: x<0, x = ',x,' slip = ',slip
                            print*,'y = ',y,'m*x = ',m*x
                            print*,'tip segment = ',n,' tip = ',tip
                            print*,'Data type = bed'
                        end if
                        if (y >= 0) then
                            signy = 1
                        else
                            signy = -1
                        end if
                        !First slip
                        slip = slip+increment
                        vx = R(n)*half_v0*((signy*(abs(y)/(m(n)*x))**xexp)+1)
                        vy = R(n)*half_v0*(m(n)/(1+s))*(((abs(y)/(m(n)*x))**yexp)-1)
                        J(1,1) = 1.-(R(n)*increment/(2*s*x))*signy*(abs(y)/(m(n)*x))**xexp !1+dvx/dx
                        J(1,2) = (R(n)*increment/(2*m(n)*s*x))*(abs(y)/(m(n)*x))**((1.-s)/s) !dvx/dy
                        J(2,1) = -((R(n)*increment*m(n))/(2*s*x))*(abs(y)/(m(n)*x))**yexp !dvy/dx
                        J(2,2) = 1.+(R(n)*increment/(2*s*x))*signy*(abs(y)/(m(n)*x))**xexp !1+dvy/dy
                        !print*,1.+(R(n)*increment/(2*s*x))*signy*(abs(y)/(m(n)*x))**xexp
                        C = matmul(matmul(J,C),transpose(J))
                        x = x+vx-R(n)*PoverS(n)*increment
                        y = y+vy
                        if (abs(slip)>abs(rem_slip)) then
                            exit
                        end if
                    end do
                end if
                if (y==m(n)*x) then
                    loc = tip_loc(n) !In the hanging wall, tip segment.
                else if (y==-m(n)*x) then
                    loc = nangles+2 !In the footwall
                else
                    loc = 0 !Unknown
                end if
                ptze(:,1) = (/x+R(n)*PoverS(n)*slip,y/) !The +R(n)*PoverS(n)*slip is because we are transforming back using the same tip position we used to transform into trishear coordinates.
                ptze = ze_to_xy2D(ptze,tip(1),tip(2),min(ramp_angle(n),pi-ramp_angle(n)),ramp_dir)!Rotate out of trishear coordinates.
                pt = ptze(:,1)
                C = matmul(matmul(Rot,C),transpose(Rot)) !Transform C into the zeta, eta coordinate system after propagation and rotation.
                C = matmul(matmul(Refl,C),Refl) !Reflect C. (This will only have an effect if ramp_dir == -1; otherwise Refl is an identity matrix.)
            else if(loc==nangles+2) then !In the footwall
                vb = v0*R(n)*PoverS(n)*(/cos_ramp(n),sin_ramp(n)/) !Velocity of trishear boundary.
                slip_tri = calc_slip(pt,tip,v0,(/0.0d0,0.0d0/),vb,tan(ramp_angle(n)-ramp_dir*phi(n)))
                if((n==nsegs .or. pt(1)*ramp_dir>bendxy(1,n)*ramp_dir) .and. slip_sign*slip_tri<slip_sign*rem_slip &
                    .and. slip_sign*slip_tri>0) then !Goes into trishear segment.
                    slip = slip_tri
                    loc = nangles+1
                else !Stays in FW.
                    slip = rem_slip
                end if
            else if (loc==tip_loc(n)) then !In the tip segment.
                v = v0*R(n)*(/cos_ramp(n),sin_ramp(n)/)
                vb = v0*R(n)*PoverS(n)*(/cos_ramp(n),sin_ramp(n)/)
                slip_tri = calc_slip(pt,tip,v0,v,vb,tan(ramp_angle(n)+ramp_dir*phi(n)))
                if (((slip_sign==-1 .and. PoverS(n)>1) .or. (slip_sign==1 .and. PoverS(n)<1)) .and. &
                    slip_sign*slip_tri<slip_sign*rem_slip .and. slip_sign*slip_tri>0) then  !Goes into the trishear segment.
                    slip = slip_tri
                    !print*,slip,PoverS(n)
                    loc = nangles+1
                else if(slip_sign==-1 .and. loc<nangles) then
                    slip_next = calc_slip(pt,bendxy_loc(:,loc),v0,v,(/0.0d0,0.0d0/),axis_slopes(loc)) !Having this only after testing slip_tri is fine because if the trishear zone can reach the point before slipping by rem_slip, then the point must still be before the fold axis when that happens. So if slip_tri*slip_sense<rem_slip*slip_sense, slip_tri*slip_sign must be < slip_next*slip_sign.
                    if(slip_next*slip_sign<rem_slip*slip_sign  .and. slip_sign*slip_next>0) then
                        slip = slip_next
                        vnext = v0*R_loc(loc+1)*(/cos_ramp_loc(loc+1),sin_ramp_loc(loc+1)/)
                        call prop_bend_error(C,v,vnext,(/0.0d0,0.0d0/),axis_slopes(loc))
                        loc = loc+1
                    else
                        slip = rem_slip
                    end if
                else
                    slip = rem_slip
                end if
                pt = pt+v*(slip/v0)
            else !In one of the fault segments
                v = v0*R_loc(loc)*(/cos_ramp_loc(loc),sin_ramp_loc(loc)/)
                if (loc == nangles .and. slip_sign == -1) then !In the last segment and moving backwards. Can't go anywhere else.
                    slip = rem_slip
                else
                    if(slip_sign==1) then !Forwards slip, can slip into segment above.
                        slip_next = calc_slip(pt,bendxy_loc(:,loc-1),v0,v,(/0.0d0,0.0d0/),axis_slopes(loc-1))
                    else    !Backwards slip, can slip into segment below.
                        slip_next = calc_slip(pt,bendxy_loc(:,loc),v0,v,(/0.0d0,0.0d0/),axis_slopes(loc))
                    end if
                    if(slip_next*slip_sign<rem_slip*slip_sign .and. slip_sign*slip_next>0) then
                        slip = slip_next
                        vnext = v0*R_loc(loc-slip_sign)*(/cos_ramp_loc(loc-slip_sign),sin_ramp_loc(loc-slip_sign)/)
                        call prop_bend_error(C,v,vnext,(/0.0d0,0.0d0/),axis_slopes(loc-(1+slip_sign)/2))
                        loc = loc-slip_sign
                    else
                        slip = rem_slip
                    end if
                end if
                pt = pt+v*(slip/v0)
            end if
            rem_slip = rem_slip - slip !remaining slip
            tip(1) = tip(1)+slip*R(n)*PoverS(n)*cos_ramp(n)
            tip(2) = tip(2)+slip*R(n)*PoverS(n)*sin_ramp(n)  
            slip = 0 !Otherwise, when we go to loc = 0, we reuse the same slip twice.
        end do
        pts(:,k) = pt
        if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
            Cbed(:,:,k) = C
        end if
    end do  
end do
end subroutine trishear_func_multi_bend

subroutine trishear_func_dip_multi_bend(pts,dips,tip_start,tip_seg_start,tip_seg_end,v0,phi,m,increment,PoverS,s,&
    bendxy,ramp_angle,ramp_dir,axis_dips,R,slipseg,sigma_dip_array)
!Note: This one takes in points in the trishear coord system around the starting tip (usually tip_final) and returns them in coords around the ending tip.
use data_module,only:ndips
use constants
use math, only: rot2D,xy_to_ze2D,ze_to_xy2D
use analytic_trishear_module, only: analytic_trishear_dip
use options, only: propagate
implicit none
double precision, dimension(1:,1:), intent(inout) :: pts
double precision,dimension(1:),intent(inout) :: dips
double precision,dimension(2) :: tip_start !Starting tip position.
integer :: tip_seg_start,tip_seg_end !Starting and ending tip segments. An alternative would be to calculate tip_seg_end within this program.
double precision, dimension(2,nbends), intent(in) :: bendxy
double precision, intent(in) :: v0,increment,s
integer,intent(in) :: ramp_dir
double precision, dimension(nsegs),intent(in) :: ramp_angle,phi,m,PoverS
double precision,dimension(ndips),intent(inout) :: sigma_dip_array !Array of sigma dip
double precision,dimension(nbends),intent(in) :: axis_dips
double precision,dimension(nsegs),intent(in) :: R !Ratio of slip on a given segment relative to slip on the lowest segment.
double precision,dimension(nsegs),intent(in) :: slipseg !The amoung of slip with the tip in each segment. (Absolute value.) In the order that the tip moves through them.
double precision,dimension(2) :: tip !The location of the fault tip. Moves as the fold is restored.
double precision,dimension(2) :: pt !The point being moved
double precision,dimension(2,1) :: ptze !The point in trishear coordinates
integer :: k,i
!double precision,dimension(nbends) :: axis_slopes !Slopes of fold axes.
double precision,dimension(nsegs) :: sin_ramp,cos_ramp !Just calculate these once for speed.
double precision :: half_v0,xexp,yexp,slip,x,y,signy
!double precision,dimension(2,2) :: C,J !Covariance matrix for the specific point being worked on and Jacobian matrix for motion in the trishear velocity field.
!double precision,dimension(2,2) :: Refl,Rot !Reflection and rotation matrixes for translations to or from the trishear coordinate system.
double precision :: rem_slip
integer :: slip_sign !1 for forwards, 2 for backwards.
double precision :: slip_tri,slip_next !Slip to get into a specific zone.
double precision :: vx, vy !Trishear velocities
double precision,dimension(2) :: v,vb,vnext !Velocity of a point and a boundary, and velocity in the next fault segment that the point is moving towards.
integer :: loc !tells which zone a point is in.
integer :: nseg_tip !number of the segments that the tip will cross through during the course of the model.
integer :: n !a counter
double precision :: dip,sigmadip !The the dip being moved and its uncertainty.
double precision :: v0term,exp1,exp2
double precision,dimension(nsegs) :: sqrtm
!Variables added or changed to accommodate the change to treat fault bends differently from changes in P/S or phi.
double precision,dimension(nangles) :: ramp_angle_loc,sin_ramp_loc,cos_ramp_loc,tan_ramp_loc,R_loc
double precision,dimension(nangles-1) :: axis_slopes !Slopes of fold axes.
double precision,dimension(2,nangles-1) :: bendxy_loc
integer,dimension(nsegs) :: tip_loc


!Separate changes in fault angle from other changes separating segments. This is added in quickly after I realized there is a problem with the way things are currently done.
!It would probably be better eventually to do this kind of thing up top or to more completely overhaul how ramp segments of different kinds are dealt with.
!The ending loc is because these go with the loc variable that denotes which fault ramp angle a point is above.
i = 1
ramp_angle_loc(1) = ramp_angle(1)
R_loc(1) = R(1)
tip_loc(1) = 1
do n = 2,nsegs
    if (ramp_angle(n) /= ramp_angle(n-1)) then
        ramp_angle_loc(i+1) = ramp_angle(n)
        R_loc(i+1) = R(n)
        bendxy_loc(:,i) = bendxy(:,n-1)
        axis_slopes(i) = tan(axis_dips(n-1))
        i = i+1
    end if
    tip_loc(n) = i !This lets us translate a tip segment to the corresponding value of the "loc" variable.
end do
sin_ramp_loc = sin(ramp_angle_loc)
cos_ramp_loc = cos(ramp_angle_loc)
tan_ramp_loc = tan(ramp_angle_loc)

!axis_slopes = tan(axis_dips)
sin_ramp = sin(ramp_angle)
cos_ramp = cos(ramp_angle)
!tan_ramp = tan(ramp_angle)
half_v0 = v0/2; !v0/2
xexp = 1/s; !exponent in vx term
yexp = (1+s)/s; !exponent in vy term
v0term = half_v0/s !v0/(2s) constant term that goes in front in dip equation
exp1 = (1+s)/(2*s) !exponent in 1st term for change in dip
exp2 = (1-s)/(2*s) !exponent in 2nd term for change in dip
sqrtm = sqrt(m) !square root of m
if (v0 >= 0) then !determine the sign of the slip
    slip_sign = 1
else
    slip_sign = -1
end if
nseg_tip = abs(tip_seg_end-tip_seg_start)+1 !Number of segments tip will pass through.
do n = tip_seg_start,tip_seg_end !Loop through each segment that the tip goes through
    do k = 1,ndips
        pt = pts(:,k)
        if (n == tip_seg_start) then
            tip = tip_start
        else
            tip = bendxy(:,n+(slip_sign-1)/2)
        end if
        slip = 0 !slip so far
        loc = 0 !location: 0 = unknown, 1 to nsegs = fault segments, nsegs+1 = trishear, nsegs+2 = footwall
        dip = dips(k)
        sigmadip = sigma_dip_array(k)
        rem_slip = slipseg(n)
        do while (rem_slip*slip_sign>0)
            !print*,loc
            if (loc == 0) then !need to find out location
                !I also need to deal with what happens if two fold axes intersect.
                if ((ramp_dir*(pt(1)-tip(1))>ramp_dir*(pt(2)-tip(2))/tan(ramp_angle(n)+ramp_dir*phi(n))) .and.&
                        ((pt(2)-tip(2))>(pt(1)-tip(1))*tan(ramp_angle(n)-ramp_dir*phi(n)))) then !Trishear zone
                    loc = nangles+1
                else
                    do i = tip_loc(n),(nangles-1) !Segment i is always the segment above bend i.
                        if (pt(1)*ramp_dir>bendxy_loc(1,i)*ramp_dir .and. pt(2)-bendxy_loc(2,i)<(pt(1)-bendxy_loc(1,i))&
                                *tan_ramp_loc(i)) then
                            loc = nangles+2 !Footwall, below segment i
                            exit
                        else if (ramp_dir*(pt(1)-bendxy_loc(1,i))>ramp_dir*(pt(2)-bendxy_loc(2,i))/axis_slopes(i)) then!In the hangingwall
                            loc = i !Above segment i
                            exit
                        end if
                    end do
                    if (loc ==0) then !If we still haven't assigned it to anywhere, the only remaining place is above or below the lowest ramp.
                        if ((nangles>1 .and. pt(2)-bendxy_loc(2,nangles-1)>(pt(1)-bendxy_loc(1,nangles-1))*tan_ramp_loc(nangles)) &
                            .or. (nangles==1 .and. pt(2)-tip(2)>(pt(1)-tip(1))*tan_ramp_loc(nangles))) then !If there's only one fault segment, then there are no bends, so no bendxy_loc
                            loc = nangles !Hanging wall, lowest segment.
                        else
                            loc = nangles+2 !Footwall, lowest segment.
                        end if
                    end if
                end if
            else if (loc == nangles+1) then !Trishear
                ptze(:,1) = pt
                ptze = xy_to_ze2D(ptze,tip(1),tip(2),min(ramp_angle(n),pi-ramp_angle(n)),ramp_dir); !Rotate into trishear coordinates.
                dip = (ramp_dir*dip+min(ramp_angle(n),pi-ramp_angle(n))) !Convert dips to ze coord. system
                x = ptze(1,1)
                y = ptze(2,1)
!                Refl = reshape((/ramp_dir,0,0,1/),(/2,2/))
!                Rot = reshape((/ramp_dir*cos_ramp(n),sin_ramp(n),-sin_ramp(n),ramp_dir*cos_ramp(n)/),&
!                    (/2,2/))
                !Note: If we bring in the points non-reflected, then we need to reflect here if necessary.
!                C = matmul(matmul(Refl,C),Refl) !Reflect C. (This will only have an effect if ramp_dir == -1; otherwise Refl is an identity matrix.)
!                C = matmul(matmul(transpose(Rot),C),Rot) !Change C into the rotated (zeta,eta) coordinate system.
                slip = 0
                if (s==1) then
                    !print*,k
                    !print*,x,y,dip,sigmadip,m(n),PoverS(n),R(n)*rem_slip,slip_sign,slip
                    call analytic_trishear_dip(x,y,dip,sigmadip,m(n),PoverS(n),R(n)*rem_slip,slip_sign,slip)
                    slip = slip/R(n) !Convert slip in the trishear segment to the slip on the lowest segment that we're using to keep track of the slip accrued.
                else
                    do while (slip == 0 .or. (y<=m(n)*x .and. y>=-m(n)*x))
                        if (x<0) then
                            print*,'Error: x<0, x = ',x,' slip = ',slip
                            print*,'y = ',y,'m*x = ',m*x
                            print*,'tip segment = ',n,' tip = ',tip
                            print*,'Data type = dip'
                        end if
                        if (y >= 0) then
                            signy = 1
                        else
                            signy = -1
                        end if
                        !First slip
                        slip = slip+increment
                        sigmadip = sigmadip*abs(1.+signy*(R(n)*increment/(s*x))*((abs(y)/(m(n)*x))**xexp)&
                            *(0.5*((x/y)-(y/x))*sin(2*dip)+cos(2*dip)))
!                        sigmadip = sqrt((sigmadip**2)*(1.+signy*(R(n)*increment/(s*x))*((abs(y)/(m(n)*x))**xexp)&
!                            *(0.5*((x/y)-(y/x))*sin(2*dip)+cos(2*dip)))**2+(50**2)*(signy*(R(n)*increment/(2*(s*x)**2))*&
!                            ((abs(y)/(m(n)*x))**xexp)*((1+2*s)*(y/x)*(cos(dip)**2)+(2+2*s)*sin(dip)*cos(dip)+&
!                            (x/y)*(sin(dip))**2))**2+(50**2)*(-signy*(R(n)*increment/(2*(s*x)**2))*((abs(y)/(m(n)*x))**xexp)&
!                            *((1+s)*(cos(dip))**2+2*(x/y)*sin(dip)*cos(dip)+(1-s)*((x/y)**2)*(sin(dip))**2))**2)
                        dip = dip+(R(n)*v0term/x)*(signy*sqrtm(n)*((abs(y)/(m(n)*x))**(exp1))*cos(dip)+(1/sqrtm(n))&
                            *((abs(y)/(m(n)*x))**(exp2))*sin(dip))**2
                        vx = R(n)*half_v0*((signy*(abs(y)/(m(n)*x))**xexp)+1)
                        vy = R(n)*half_v0*(m(n)/(1+s))*(((abs(y)/(m(n)*x))**yexp)-1)
                        x = x+vx-R(n)*PoverS(n)*increment
                        y = y+vy
                        if (abs(slip)>abs(rem_slip)) then
                            exit
                        end if
                    end do
                end if
                if (y==m(n)*x) then
                    loc = tip_loc(n) !In the hanging wall, tip segment.
                else if (y==-m(n)*x) then
                    loc = nangles+2 !In the footwall
                else
                    loc = 0 !Unknown
                end if
                ptze(:,1) = (/x+R(n)*PoverS(n)*slip,y/) !The +PoverS(n)*slip is because we are transforming back using the same tip position we used to transform into trishear coordinates.
                ptze = ze_to_xy2D(ptze,tip(1),tip(2),min(ramp_angle(n),pi-ramp_angle(n)),ramp_dir)!Rotate out of trishear coordinates.
                pt = ptze(:,1)
                dip = ramp_dir*(dip - min(ramp_angle(n),pi-ramp_angle(n))) !Convert dip back into section coordinates.
!                C = matmul(matmul(Rot,C),transpose(Rot)) !Transform C into the zeta, eta coordinate system after propagation and rotation.
!                C = matmul(matmul(Refl,C),Refl) !Reflect C. (This will only have an effect if ramp_dir == -1; otherwise Refl is an identity matrix.)
            else if(loc==nangles+2) then !In the footwall
                vb = v0*R(n)*PoverS(n)*(/cos_ramp(n),sin_ramp(n)/) !Velocity of trishear boundary.
                slip_tri = calc_slip(pt,tip,v0,(/0.0d0,0.0d0/),vb,tan(ramp_angle(n)-ramp_dir*phi(n)))
                if((n==nsegs .or. pt(1)*ramp_dir>bendxy(1,n)*ramp_dir) .and. slip_sign*slip_tri<slip_sign*rem_slip &
                    .and. slip_sign*slip_tri>0) then !Goes into trishear segment.
                    slip = slip_tri
                    loc = nangles+1
                else !Stays in FW.
                    slip = rem_slip
                end if
            else if (loc==tip_loc(n)) then !In the tip segment.
                v = v0*R(n)*(/cos_ramp(n),sin_ramp(n)/)
                vb = v0*R(n)*PoverS(n)*(/cos_ramp(n),sin_ramp(n)/)
                slip_tri = calc_slip(pt,tip,v0,v,vb,tan(ramp_angle(n)+ramp_dir*phi(n)))
                if (((slip_sign==-1 .and. PoverS(n)>1) .or. (slip_sign==1 .and. PoverS(n)<1)) .and. &
                    slip_sign*slip_tri<slip_sign*rem_slip .and. slip_sign*slip_tri>0) then !Goes into the trishear segment.
                    slip = slip_tri
                    loc = nangles+1
                else if(slip_sign==-1 .and. loc<nangles) then
                    slip_next = calc_slip(pt,bendxy_loc(:,loc),v0,v,(/0.0d0,0.0d0/),axis_slopes(loc))
                    if(slip_next*slip_sign<rem_slip*slip_sign  .and. slip_sign*slip_next>0) then
                        slip = slip_next !Goes into next segment.
                        vnext = v0*R_loc(loc+1)*(/cos_ramp_loc(loc+1),sin_ramp_loc(loc+1)/)
                        call calc_dip_change(dip,sigmadip,v,vnext,(/0.0d0,0.0d0/),axis_slopes(loc))
                        loc = loc+1
                    else !Stays in this fault segment.
                        slip = rem_slip
                    end if
                else
                    slip = rem_slip
                end if
                pt = pt+v*(slip/v0)
            else !In one of the fault segments
                v = v0*R_loc(loc)*(/cos_ramp_loc(loc),sin_ramp_loc(loc)/)
                if (loc == nangles .and. slip_sign == -1) then !In the last segment and moving backwards. Can't go anywhere else.
                    slip = rem_slip
                else
                    if(slip_sign==1) then !Forwards slip, can slip into segment above.
                        slip_next = calc_slip(pt,bendxy_loc(:,loc-1),v0,v,(/0.0d0,0.0d0/),axis_slopes(loc-1))
                    else    !Backwards slip, can slip into segment below.
                        slip_next = calc_slip(pt,bendxy_loc(:,loc),v0,v,(/0.0d0,0.0d0/),axis_slopes(loc))
                    end if
                    if(slip_next*slip_sign<rem_slip*slip_sign  .and. slip_sign*slip_next>0) then
                        slip = slip_next
                        vnext = v0*R_loc(loc-slip_sign)*(/cos_ramp_loc(loc-slip_sign),sin_ramp_loc(loc-slip_sign)/)
                        call calc_dip_change(dip,sigmadip,v,vnext,(/0.0d0,0.0d0/),axis_slopes(loc-(1+slip_sign)/2))
                        loc = loc-slip_sign
                    else
                        slip = rem_slip
                    end if
                end if
                pt = pt+v*(slip/v0)
            end if
            rem_slip = rem_slip - slip !remaining slip
            tip(1) = tip(1)+slip*R(n)*PoverS(n)*cos_ramp(n)
            tip(2) = tip(2)+slip*R(n)*PoverS(n)*sin_ramp(n)
            slip = 0 !Otherwise, when we go to loc = 0, we reuse the same slip twice.
        end do
        pts(:,k) = pt
        dips(k) = dip
        if (propagate == 1) then !It might be good to not calculate the propagation in the first place if we're not using propagated errors, but this will work for now.
            sigma_dip_array(k) = sigmadip
        end if
    end do
end do
end subroutine trishear_func_dip_multi_bend

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
double precision,dimension(2),intent(in) :: vb !Velocity at which ptb (and therefore the boundary) is moving
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

end module trishear_multibend