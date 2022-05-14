program Setup_Trishear
!A program to allow the user to easily set up an options file for running the Trishear program.

implicit none
character (len=50) :: save_file !options file that everything will be saved to
character (len=50) :: file_name !generic file name
logical :: ex !tells whether the file named exists
integer :: save_id = 4
integer :: ndatatypes
integer,dimension(:),allocatable :: datalist
integer :: goodinput !Tells whether the input is an acceptable value
integer :: i !A counter
integer,dimension(5) :: DatTypes !Types of data
integer :: DatType_beds=1,DatType_dips=2,DatType_terr=3,DatType_faultpts=4,DatType_beds_restored = 5 !Tells which number (and position in DatTypes) corresponds to which data type.
integer :: FaultType
integer :: nparams
integer :: propagate !Tells whether or not to propagate errors.
integer :: decoly_at_tipy
integer :: bend_at_tipy !Tells whether decollement depth and tipy must be the same (1) or not (0)
integer :: tipy_can_be_below_bend !Tells whether tipy (initial) can be below the bend
integer :: bend_start !The bend at which the tip must start if bend_at_tipy == 1
integer :: tipy_above_bend !Tells whether the initial tipy is forced to be above one of the bends, which is specified.
integer :: bend_tip_above !The number of the bend that the tip must be above if tipy_above_bend == 1.
integer :: nbends !Tells the total number of fault bends / points of change in ramp angle, P/S, or phi.
integer :: nsegs !Tells the total number of fault segments. This will be nbends+1
integer :: nangles,nPoverS,nphi !Numbers of different ramp angle, PoverS, and phi values. The sum of these three is nsegs.
integer :: backlimb_fold_type !(1) Fault parallel flow, (2) Fault bend folding, (3) Inclined Shear
integer :: phi_le_ramp !Tells whether or not phi is required to be less than the ramp angle in any segment that the tip passes through. (For multibend fault model.)
integer :: let_trishear_intersect_axes !Tells whether the trishear zone is allowed to intersect the fold axes above some specified elevation. (1 = yes, 0 = no) (For multibend fold model.)
double precision :: intersect_elevation !If trishear zone is allowed to intersect fold axes above some elevation, this is the elevation. (For multibend fault model.)
!double precision :: shear_angle !For inclined shear. 0 deg is horizontal. 90 deg is vertical. Only antithetic shear is allowed.
integer :: final_tip_in_upper_segment !Requires the final position of the fault tip to be in the uppermost fault segment, so that there are no imaginary, unused fault segments above it.
integer :: method
integer(kind = 8) :: nmodels !number of models to run
integer :: init_mod_rand !For MCMC, tells if the initial model is (1) random or (2) specified by the user.
integer :: nparams_total !For APT: The total number of parameters. When APT options is called, some parameters such as restored bed depths may not have been added into the nparams variable yet. This corrects for that omisison.
double precision, dimension(:),allocatable :: init_model !For MCMC, if initial model is not random, contains the initial parameters
integer :: begin_adapt !For AM algorithm, model number at which to begin adapting the covariance matrix
integer nlevels
integer :: save_at_interval,nsave !For APT algorithm, tells whether to save only at some interval and what that interval should be.
!integer :: max_iter !For RML algorithm, maximum number of iterations
!double precision,dimension(:),allocatable :: tol !For RML algorithm, parameter tolerances
!integer :: initrand !For RML algorithm, tells whether the initial simplex should be random
integer :: TipToSolve
integer :: constrain_nonparam_tip !Tells whether or not to constrain the tip position (initial or final) that isn't fit for as a model parameter.
integer :: constraint_type !Tells which type of constraint to use if the non-parameter tip position is being constrained.
double precision :: np_tip_xmin,np_tip_xmax,np_tip_ymin,np_tip_ymax !Limits for the tip position that isn't being fit for as a model parameter to use if constrain_nonparam_tip == 1.
double precision :: constraint_dip,constraint_intercept !dip, slope, and intercept of a line used as a constraint on the initial fault tip.
integer :: tip_relative_to_line !Tells the position of the tip relative to the line if constraint_type is a line. 1 = above (or left for vertical lines) and 2 = below (or right for vertical lines).
integer :: FitType !Type of fit to perform to data
double precision :: RegDip !Regional Dip
double precision,dimension(:),allocatable :: intercepts !y intercepts of known lines to fit beds to.
integer :: nbeds
integer :: ResultType !Tells what type of result the run should produce
double precision :: sigma_dip,sigma_fault !Bed, dip, and fault uncertainties
double precision,dimension(:),allocatable :: sigma_bed,sigma_bed_x,sigma_bed_y !Separate x and y uncertainties in the position of a point for beds.
double precision :: sigma_bed_restored,sigma_bed_restored_x,sigma_bed_restored_y !Uncertainties in points expected to be on restored beds.
integer :: sigma_xy_diff !Tells if sigma x and sigma y are different for beds.
integer :: diff_sigma_bed !Tells whether there is a different uncertainty for each bed.
double precision,dimension(:),allocatable :: sigma_terr !Uncertainty in terrace points. This can be different for each terrace.
double precision,dimension(:),allocatable :: sigma_orig_elev !Uncertainty in the restored inner edge elevation(s) of terraces (i.e. uncertainty in paleosealevel).
double precision :: sigma_restored_dip !Uncertainty in the restored dip values, that is not propagated through.
integer :: use_sigma_restored_dip !Tells whether or not to use sigma_restored_dip (0/1).
integer :: state_for_Lc !Tells whether distance for beds correlated along length should be calculated in the deformed state (1) or the restored state (2).
integer :: fit_Lc !Tells whether to fit for Lc as a model parameter.
double precision :: Lc !correlation length
integer :: nterr !Number of terraces
double precision,dimension(:),allocatable :: terr_dip_deg !Terrace surface dip
double precision,dimension(:),allocatable :: terr_orig_elev !Terrace inner edge original elevations
integer :: terr_age_order !Tells if terraces are in order by age: 0 = no order; 1 = youngest to oldest; 2 = oldest to youngest
integer :: beds_age_order !Tells if beds are in order by age: 0 = no order; 1 = youngest to oldest; 2 = oldest to youngest
integer :: fit_dips !Tells whether we should not fit for beds as a parameter (0), fit for a single restored dip for all data (1), or fit for each bed's restored dip individually (2).
integer :: n_bed_groups !Tells the number of groups of beds with the same dip, if fit_dips==3.
integer :: bed_group_start !Tells on which bed number each group starts.
integer :: fit_beds !Tells whether we should fit for bed intercepts as a parameter (0=no, 1=yes).
integer :: n_bed_segs !Number of segments per bed if fitting to multi-segment restored beds.
integer :: fault_shape !Tells the shape of the fault to be approximated by segments. 1 = circle, 2 = ellipse.
integer :: fit_sigma !Tells whether we should fit for data uncertainties as model parameters.
integer :: has_growth !Tells if the bed contains growth strata.
integer :: growth_slip_parameterization
character (len=50) :: errs_file_name
integer :: n_sigma_bed_groups, sigma_bed_group_start
integer :: ellipse_parameterization !How to parameterize the ellipse. (1) y_base, a, b; (2) y_base, y_str, e
integer :: limit_fault_x_pts !Tells whether or not to put limits on the x positions of fault points.
double precision,dimension(2) :: fault_x_limits !Minimum and maximum allowed value for fault points.
integer :: fault_concave_up !Tells whether the fault is required to be concave up (1) or can have convex bends (0).
integer :: beds_order_type !Tells how the order of beds is determined.
double precision,dimension(2) :: bed_domain_x !x coordinates of the domain in which to test for crossing.
integer :: nknots !Number of knots (control points) for the spline. This includes the fault tip.
integer :: constrain_backlimb_syncline !Tells whether to put constraints on x position of backlimb syncline (0 or 1).
double precision :: y_backlimb_syncline !Elevation at which to put the constraints on the x position of the backlimb syncline.
double precision,dimension(2) :: x_backlimb_syncline_limits !Minimum and maximum allowed x positions of backlimb syncline at y = y_backlimb_syncline
integer :: limit_bed_dips=0 !If FitType is 3, this tells whether or not to put min and max limits on the bed dips (0 or 1).
double precision,dimension(:),allocatable :: min_bed_dips,max_bed_dips !If limit_bed_dips is 1, these give the minimum and maximum allowed dip for each bed.

!Get the file to save options to.
print*, 'Enter name of a file to save the options to: '
read(*,*) save_file
open(unit=save_id, file = save_file)

!Get Data Type
print*,'Enter number of different data types: '
read(*,*) ndatatypes
write(save_id,*) ndatatypes
allocate(datalist(ndatatypes))
do 
    print*,'Data Are (list all that apply): '
    print*,'(1) beds '
    print*,'(2) dips '
    print*,'(3) points on a fault'
    print*,'(4) marine terraces'
    print*,'(5) points on restored state beds'
    read(*,*) datalist
    DatTypes = 0 !Initialize them all to zero
    goodinput = 1 !Initialize to 1
    do i=1,ndatatypes
        select case(datalist(i))
            case(1)
            DatTypes(DatType_beds) = 1
            case(2)
            DatTypes(DatType_dips) = 1
            case(3)
            DatTypes(DatType_faultpts) = 1
            case(4)
            DatTypes(DatType_terr) = 1
            case(5)
            DatTypes(DatType_beds_restored) = 1
            case default
            goodinput = 0 !Not good input. Needs to be redone.
        end select
    end do
    if (goodinput == 1) then !If all the inputs were good, we're finished here.
        exit
    else
        print*,'Invalid Command'
    end if
end do
write(save_id,*) datalist

!Read Beds
if (DatTypes(DatType_beds) == 1) then
    do !Do statement for reading bed file name
        print*,'Input bedding file name (with extension)'
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    write(save_id,*) file_name
end if

!Find out about growth strata parameterization if neceessary.
if (DatTypes(DatType_beds) == 1) then
    do
        print*,'Do the beds include growth strata? (0/1)'
        read(*,*) has_growth
        if (has_growth == 1) then
            do
                print*,'How is slip on growth strata parameterized?'
                print*,'    (1) Slip is total for each growth layer.'
                print*,'    (2) Slip is additional since previous layer.'
                print*,'        (In order in data file from top to bottom. Put youngest on top if using this option).'
                read(*,*) growth_slip_parameterization
                if (growth_slip_parameterization==1 .or. growth_slip_parameterization==2) then
                    write(save_id,*) 'growth_slip_parameterization'
                    write(save_id,*) growth_slip_parameterization
                    exit
                else
                    print*,'Error: Unrecognized growth strata slip parameterization.'
                end if
            end do
            exit
        else if (has_growth == 0) then
            exit
        else
            print*,'Error: Unrecognized response. Answer with 0 = False or 1 = True.'
        end if
    end do
end if

!Read Dips
if (DatTypes(DatType_dips) == 1) then
    do !Do statement for reading bed file name
        print*,'Input dip file name (with extension)'
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    write(save_id,*) file_name
    !Find out if height or depth
end if

!Read Terraces
if (DatTypes(DatType_terr) == 1) then
    do !Do statement for reading bed file name
        print*,'Input terrace file name (with extension)'
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    write(save_id,*) file_name
end if

!Read Fault Points
if (DatTypes(DatType_faultpts) == 1) then
    do !Do statement for reading fault points file name
        print*,'Input fault file name (with extension)'
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    write(save_id,*) file_name
end if

!Read Restored Beds
if (DatTypes(DatType_beds_restored) == 1) then
    do !Do statement for reading fault points file name
        print*,'Input restored bedding file name (with extension)'
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    write(save_id,*) file_name
end if

!Get Tip Type
do
    print*,'Tip position is:'
    print*,'(1) Initial'
    print*,'(2) Final'
    read(*,*) TipToSolve
    if (TipToSolve == 1 .or. TipToSolve == 2) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
write(save_id,*) TipToSolve
do
    if (TipToSolve==1) then
        print*,'Do you want to place constraints on the final tip position too? (0/1)'
    else
        print*,'Do you want to place constraints on the initial tip position too? (0/1)'
    end if
    read(*,*) constrain_nonparam_tip
    if (constrain_nonparam_tip==1) then
        print*,'Warning: Constraining non-parameter tip position will only work if you are using a Straight, Multi-bend, &
            Listric (approximate), or Spline fault type.'
        do
            print*,'Do you want to:'
            print*,'(1) Constrain the tip position within a box?'
            print*,'(2) Constrain the tip position relative to a line?'
            read(*,*) constraint_type
            if (constraint_type == 1) then
                print*,'Enter tip position limits as: xmin,xmax,ymin,ymax.'
                read(*,*) np_tip_xmin,np_tip_xmax,np_tip_ymin,np_tip_ymax
                write(save_id,*) 'constrain_non-parameter_tip'
                write(save_id,*) constraint_type
                write(save_id,*) np_tip_xmin,np_tip_xmax,np_tip_ymin,np_tip_ymax
                exit
            else if (constraint_type == 2) then
                print*,'Enter line as: line dip (in degrees, positive down to left), y-intercept.'
                read(*,*) constraint_dip,constraint_intercept
                write(save_id,*) 'constrain_non-parameter_tip'
                write(save_id,*) constraint_type
                write(save_id,*) constraint_dip,constraint_intercept
                if (constraint_dip==90) then !Vertical line
                    print*,'Must the tip be (1) left or (2) right of the line?'
                else if (constraint_dip>0) then
                    print*,'Must the tip be (1) above / left of the line or (2) below / right of the line?'
                else if (constraint_dip<0) then
                    print*,'Must the tip be (1) above / right of the line or (2) below / left of the line?'
                else !Horizontal line
                    print*,'Must the tip be (1) above or (2) below the line?'
                end if
                read(*,*) tip_relative_to_line
                write(save_id,*) tip_relative_to_line
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        exit
    else if (constrain_nonparam_tip==0) then
        exit
    else
        print*,'Invalid Command'
    end if
end do

!Get Fault Type
do
    print*,'Choose Fault Type'
    print*,'(1) Straight Fault Ramp'
    print*,'(2) Ramp from Horizontal Detachment'
    print*,'(3) Fault with a bend in it'
    print*,'(4) Listric Fault (fault parallel flow)'
    print*,'(5) Parallel Fault Propagation Fold'
    print*,'(6) Multi-bend fault'
    print*,'(7) Listric Fault (circular or elliptic, approximate)'
    print*,'(8) Spline Fault (Requires tip position is final)'
    read(*,*) FaultType
    if (1<=FaultType .and. FaultType<=7) then
        exit
    else if (FaultType==8 .and. TipToSolve==2) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
write(save_id,*) FaultType

!Options based on fault type
select case(FaultType)
    case(1) !Straight Fault
        nparams = 7
    case(2) !Horizontal Detachment
        do
            print*,'Must detachment depth be at initial fault tip depth? (0/1)'
            read(*,*) decoly_at_tipy
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
        write(save_id,*) decoly_at_tipy
    case(3) !Bend in the fault
        do
            if (TipToSolve==1) then
                print*,'Must depth of fault bend be at initial fault tip depth? (0/1)'
                read(*,*) bend_at_tipy
            else
                bend_at_tipy = 0
            end if
            tipy_can_be_below_bend = 0
            if (bend_at_tipy == 0) then
                nparams = 9
                do
                    print*,'Can initial tipy be below fault bend? (0/1)'
                    read(*,*) tipy_can_be_below_bend
                    if (tipy_can_be_below_bend == 0 .or. tipy_can_be_below_bend == 1) then
                        exit
                    else
                        print*,'Invalid Command'
                    end if
                end do
                exit
            else if (bend_at_tipy == 1) then
                nparams = 8
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        if (TipToSolve == 1) then
            write(save_id,*) bend_at_tipy
        end if
        if(bend_at_tipy == 0) then
            write(save_id,*) tipy_can_be_below_bend
        end if
    case(4) !Listric fault
        do
            if (TipToSolve==1) then !If solving for the initial tip
                print*,'Must detachment depth be at initial fault tip depth? (0/1)'
                read(*,*) decoly_at_tipy
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
        if (TipToSolve==1) then
            write(save_id,*) decoly_at_tipy
        end if
    case(5) !Parallel fault propagation fold.
        do
            print*,'Must detachment depth be at initial fault tip depth? (0/1)'
            print*,'Choosing 0 may not work properly.'
            read(*,*) bend_at_tipy
            if (bend_at_tipy == 0) then
                nparams = 6
                exit
            else if(bend_at_tipy == 1) then
                nparams = 5
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        write(save_id,*) bend_at_tipy
    case(6,7,8) !(6) Multi-bend fault or (7) Elliptical listric fault or (8) spline fault approximated by a multi-bend fault.
        if (FaultType == 7) then
            do
                print*,'Shape: '
                print*,'(1) Circle'
                print*,'(2) Ellipse'
                read(*,*) fault_shape
                if (fault_shape == 1) then !Circle
                    write(save_id,*) 'circle'
                    exit
                else if (fault_shape == 2) then !Ellipse
                    write(save_id,*) 'ellipse'
                    exit
                else
                    print*,'Invalid Command'
                end if
            end do
            if (fault_shape==2) then
                do
                    print*,'How should the ellipse geometry be parameterized?'
                    print*,'(Note: Base and top elevation are where the becomes straight.)'
                    print*,'    (1) base elevation, horizontal semiaxis, vertical semiaxis'
                    print*,'    (2) base elevation, top elevation, eccentricity'
                    read(*,*) ellipse_parameterization
                    if (ellipse_parameterization==1 .or. ellipse_parameterization==2) then
                        write(save_id,*) 'ellipse_parameterization'
                        write(save_id,*) ellipse_parameterization
                        if (ellipse_parameterization==2) then
                            print*,'For the purpose of this program, ellipticity should be positive if the horizontal semiaxis is&
                                 longer'
                            print*,'and negative if the vertical semiaxis is longer.'
                        end if
                        exit
                    else
                        print*,'Error: Unrecognized ellipse parameterization.'
                    end if
                end do
            end if
            print*,'Note: number of fault segments is the number of segments the listric fault will be divided into.'
            print*,'This includes the horizontal detachment and the upper straight segment.'
        end if
        if (TipToSolve==1) then
            !See if the fault should begin propagating from the detachment or above it.
            !This is only an option if fitting for the initial fault tip position, since in that case tipy will also be the detachment depth so it's easier that way.
            do
                print*,'Must detachment depth be at initial fault tip depth? (0/1)'
                read(*,*) decoly_at_tipy
                if (decoly_at_tipy==0) then
                    exit
                else if(decoly_at_tipy==1) then
                    write(save_id,*) 'decoly_at_tipy'
                    exit            
                else
                    print*,'Error: Unrecognized option.'
                end if
            end do
        end if
        if (FaultType == 8) then
            print*,'Enter number of spline knots (control points).'
            print*,'This includes the bottom end point, but not the fault tip.'
            do
                read(*,*) nknots
                if (nknots>2) then
                    write(save_id,*) nknots
                    exit
                else
                    print*,'Error: You must have at least 2 knots'
                end if
            end do
            !Otherwise, this uses the same options as trishear_multibend.
            print*,'Note: number of fault segments is the number of segments the ellipse or circle will be divided into.'
            print*,'This includes the horizontal detachment and the upper straight segment.'
        end if
        print*,'Enter number of fault segments with different ramp angles'
        read(*,*) nangles
        write(save_id,*) nangles
        print*,'Enter number of different phi values'
        read(*,*) nphi
        write(save_id,*) nphi
        print*,'Enter number of different P/S values'
        read(*,*) nPoverS
        write(save_id,*) nPoverS
        nsegs = nangles+nPoverS+nphi-2
        nbends = nsegs-1
        do
            print*,'Choose type of backlimb deformation:'
            print*,'(1) Fault-Parallel Flow (Fold axes bisect fault bends)'
            print*,'(2) Fault Bend Folding (Preserves line length)'
            print*,'(3) Inclined Shear'
            read(*,*) backlimb_fold_type
            if (backlimb_fold_type==1) then
                write(save_id,*) backlimb_fold_type
                exit
            else if(backlimb_fold_type==2) then
                print*,'Warning: Restored state beds are assumed to be flat for calculating fold axis orientation'
                write(save_id,*) backlimb_fold_type
                exit
            else if(backlimb_fold_type==3) then
                !print*,'Enter shear angle in degrees (range: [0,90], 90 deg = vertical):'
                print*,'Shear angle should be in degrees (range: [0,90], 90 deg = vertical):'
                print*,'(Only antithetic shear is allowed)'
                !read(*,*) shear_angle
                write(save_id,*) backlimb_fold_type
                !write(save_id,*) shear_angle
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do
            print*,'Must fault tip stay above a specified fault bend at all times? (0=no, 1=yes)'
            read(*,*) tipy_above_bend
            if (tipy_above_bend == 0) then
                write(save_id,*) tipy_above_bend
                exit
            else if (tipy_above_bend == 1) then
                do
                    print*,'Enter bend number that tip must stay above: (1 is the highest bend.)'
                    read(*,*) bend_tip_above
                    if (bend_tip_above >=1 .and. bend_tip_above <= nbends) then
                        write(save_id,*) tipy_above_bend
                        write(save_id,*) bend_tip_above
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
        if (TipToSolve == 1 .and. FaultType == 6) then !Right now, we're only allowing the initial fault tip position to be at one of the bends if the user is searching for initial tip position.
            do
                print*,'Must initial fault tip position be at one of the bends? (0/1)'
                read(*,*) bend_at_tipy
                if (bend_at_tipy == 0) then
                    write(save_id,*) bend_at_tipy
                    nparams = nbends+nsegs
                    exit
                else if (bend_at_tipy == 1) then
                    nparams = nbends+nsegs-1
                    do
                        print*,'Enter bend number that tip must start at: (1 is the highest bend.)'
                        read(*,*) bend_start
                        if (bend_start >= 1 .and. bend_start <= nbends) then
                            write(save_id,*) bend_at_tipy
                            write(save_id,*) bend_start
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
            nparams = nbends+nsegs
        end if
        do
            print*,'Must phi be less than the ramp angle in any fault segment the tip passes through? (0=no, 1=yes)'
            print*,'(This prevents lowering of the footwall due to a trishear zone extending downwards)'
            read(*,*) phi_le_ramp
            if (phi_le_ramp == 0 .or. phi_le_ramp == 1) then
            write(save_id,*) phi_le_ramp
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do
            print*,'Can the trishear zone boundary intersect the backlimb fold axes above some elevation? (0 = no, 1 = yes)'
            read(*,*) let_trishear_intersect_axes
            if (let_trishear_intersect_axes == 0) then
                exit
            else if(let_trishear_intersect_axes == 1) then
                write(save_id,*) 'let_trishear_intersect_axes'
                print*,'Enter elevation above which trishear zone and fold axes can intersect: (must be above all data)'
                read(*,*) intersect_elevation
                write(save_id,*) intersect_elevation
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do
            print*,'Must the final fault tip position be in the uppermost fault segment? (0 = no, 1= yes)'
            print*,'(If not, some fault segments may go unused.)'
            read(*,*) final_tip_in_upper_segment
            if (final_tip_in_upper_segment == 1) then
                write(save_id,*) 'final_tip_in_upper_segment'
            end if
            if (final_tip_in_upper_segment==0 .or. final_tip_in_upper_segment==1) then
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do
            print*,'Do you want to limit the minimum and maximum x values of fault bends? (0/1)'
            read(*,*) limit_fault_x_pts
            if (limit_fault_x_pts==1) then
                write(save_id,*) 'limit_fault_x'
                print*,'Enter minimum and maximum allowed x values for fault bends.'
                read(*,*) fault_x_limits
                write(save_id,*) fault_x_limits
                exit
            else if (limit_fault_x_pts==0) then
                exit
            else
                print*,'Error: Invalid Response'
            end if
        end do
        do
            print*,'Must the fault be concave upward at all points? (0 = no, 1= yes (default))'
            read(*,*) fault_concave_up
            if (fault_concave_up==1) then
                write(save_id,*) 'fault_concave_up'
                exit
            else if (fault_concave_up==0) then
                write(save_id,*) 'fault_concave_up_not_required'
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do 
            print*,'Do you want to put constraints on the x position of the backlimb syncline? (0/1)'
            read(*,*) constrain_backlimb_syncline
            if (constrain_backlimb_syncline == 1) then
                write(save_id,*) 'constrain_backlimb_syncline'
                print*,'Enter elevation (y coordinate) at which to measure the backlimb syncline position.'
                read(*,*) y_backlimb_syncline
                write(save_id,*) y_backlimb_syncline
                print*,'Enter minimum and maximum allowed x position of backlimb syncline at y.'
                read(*,*) x_backlimb_syncline_limits
                write(save_id,*) x_backlimb_syncline_limits
                exit
            else if (constrain_backlimb_syncline == 0) then
                exit
            else
                print*,'Invalid Command'
            end if
        end do
end select
!Get parameters file name
do
    print*,'Input Parameter File Name (with extension)'
    read(*,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
write(save_id,*) file_name

!Get Method
do
    print*,'Choose Method: '
    print*,'(1) Grid Search'
    print*,'(2) Grid Monte Carlo'
    print*,'(3) Monte Carlo from a normal distribution'
    print*,'(4) Metropolis-Hastings Algorithm'
    print*,'(5) Adaptive Metropolis'
    print*,'(6) Robust Adaptive Metropolis'
    print*,'(7) Adaptive Parallel Tempering'
    !print*,'(8) Reverse Maximum Likelihood'
    read(*,*) method
    if (method>0 .and. method<8) then
        write(save_id,*) method
        exit
    else
        print*,'Invalid Command'
    end if
end do

!Get method-specific options
select case(method)
case(1) !Grid
    !No special options needed.
case(2) !Grid MC
    print*,'Enter number of models to run: '
    read(*,*) nmodels
    write(save_id,*) nmodels
case(3) !normal distribution MC
    print*,'Enter number of models to run: '
    read(*,*) nmodels
    do
        print*,'Enter filename to read proposal distribution from: '
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    write(save_id,*) nmodels
    write(save_id,*) file_name
case(4) !MCMC
    print*,'Enter number of models to run: '
    read(*,*) nmodels
    write(save_id,*) nmodels
    do
        print*,'Should initial model be (1) random or (2) specified?'
        read(*,*) init_mod_rand
        if (init_mod_rand == 1) then
            write(save_id,*) init_mod_rand
            exit
        else if (init_mod_rand == 2) then
            print*,'Enter initial values for all ',nparams,'  parameters: '
            allocate(init_model(nparams))
            read(*,*) init_model
            write(save_id,*) init_mod_rand
            write(save_id,*) init_model
            exit
        else
            print*,'Invalid Command'
        end if
    end do
case(5) !AM
    print*,'Enter number of models to run: '
    read(*,*) nmodels
    write(save_id,*) nmodels
    do
        print*,'Should initial model be (1) random or (2) specified?'
        read(*,*) init_mod_rand
        if (init_mod_rand == 1) then
            write(save_id,*) init_mod_rand
            exit
        else if (init_mod_rand == 2) then
            print*,'Enter initial values for all ',nparams,'  parameters: '
            allocate(init_model(nparams))
            read(*,*) init_model
            write(save_id,*) init_mod_rand
            write(save_id,*) init_model
            exit
        else
            print*,'Invalid Command'
        end if
    end do
    print*,'Enter model number after which to begin using adapted covariance matrix'
    read(*,*) begin_adapt
    write(save_id,*) begin_adapt
case(6) !RAM
    print*,'Enter number of models to run: '
    read(*,*) nmodels
    write(save_id,*) nmodels
    do
        print*,'Should initial model be (1) random or (2) specified?'
        read(*,*) init_mod_rand
        if (init_mod_rand == 1) then
            write(save_id,*) init_mod_rand
            exit
        else if (init_mod_rand == 2) then
            print*,'Enter initial values for all ',nparams,'  parameters: '
            allocate(init_model(nparams))
            read(*,*) init_model
            write(save_id,*) init_mod_rand
            write(save_id,*) init_model
            exit
        else
            print*,'Invalid Command'
        end if
    end do
case(7) !APT
    print*,'Enter number of models to run: '
    read(*,*) nmodels
    write(save_id,*) nmodels
    do
        print*,'Do you want to save only every n models? (0/1)'
        read(*,*) save_at_interval
        if (save_at_interval == 1) then
            write(save_id,*) 'save_at_interval'
            print*,'Enter interval at which to save model results:'
            read(*,*) nsave
            write(save_id,*) nsave
            exit
        else if (save_at_interval == 0) then
            exit
        else
            print*,'Error: Invalid Response'
        end if
    end do
    do
        print*,'Should initial model be (1) random or (2) specified?'
        read(*,*) init_mod_rand
        if (init_mod_rand == 1) then
            write(save_id,*) init_mod_rand
            exit
        else if (init_mod_rand == 2) then
            print*,'Enter number of model parameters: (currently expecting at least ',nparams,')'
            read(*,*) nparams_total
            write(save_id,*) 'nparams_total'
            write(save_id,*) nparams_total
            print*,'Enter initial values for all ',nparams_total,'  parameters: '
            allocate(init_model(nparams))
            read(*,*) init_model
            write(save_id,*) init_mod_rand
            write(save_id,*) init_model
            exit
        else
            print*,'Invalid Command'
        end if
    end do
    do
        print*,'Enter number of temperature levels to use: '
        read(*,*) nlevels
        if (nlevels < 2) then
            print*,'Error: Must have 2 or more levels for APT algorithm.'
        else
            exit
        end if
    end do
    write(save_id,*) nlevels
case default
    print*,'Error: Unrecognized Method'
end select

!Find out the type of fit to do
do
    print*,'Objective Function: '
    print*,'(1) Fit to flat line.'
    print*,'(2) Fit to known dip.'
    print*,'(3) Fit to best fit line / mean dip.'
    if (DatTypes(1) == 1) then !only useful if beds are included
        print*,'(4) Fit to known line'
    end if
    if (DatTypes(1)==1 .or. DatTypes(2)==1) then
        print*,'(5) Fit for restored dips and/or bed elevations as model parameters'
    end if
    if (DatTypes(DatType_beds)==1) then
        print*,'(6) Multi-segment restored beds'
    end if
    read(*,*) FitType
    if (FitType == 1) then
        print*,'Note: This method fits to an average (unweighted) of the restored bed points'
        write(save_id,*) FitType
        exit
    else if (FitType == 2) then
        print*,'Note: This method fits to an average (unweighted) of the restored bed intercepts corresponding to the given slope'
        print*,'Enter regional dip in degrees. (Negative is down to the left, positive right)'
        read(*,*) RegDip
        write(save_id,*) FitType
        write(save_id,*) RegDip
        exit
    else if (FitType == 3) then
        print*,'Note: This method uses an unweighted least squares fit'
        if (DatTypes(1) == 1) then !Data are beds
            print*,'Warning: If there is more than one bed, each bed will be fit separately for dip.'
        end if
        if (DatTypes(1) == 2 .and. DatTypes(2) == 1) then !Data are both beds and dips
            print*,'Warning: Bed and dip data will be fit separately, not necessarily to the same dip.'
        end if
        write(save_id,*) FitType
        print*,'Do you want to put limits on restored bed dips? (0=no, 1=yes)'
        read(*,*) limit_bed_dips
        if (limit_bed_dips==1) then
            write(save_id,*) 'limit_bed_dips'
            print*,'Enter number of beds: '
            read(*,*) nbeds
            allocate(min_bed_dips(nbeds),max_bed_dips(nbeds))
            print*,'Note: Bed dips are negative down to left and positive down to right.'
            print*,'Enter minimum dips for all ',nbeds,' beds:'
            read(*,*) min_bed_dips
            write(save_id,*) min_bed_dips
            print*,'Enter maximum dips for all ',nbeds,' beds:'
            read(*,*) max_bed_dips
            write(save_id,*) max_bed_dips
        end if
        exit
    else if (FitType == 4) then
        if (sum(DatTypes)>1) then
            print*,'Note: Only beds will be fit to specified lines.'
            print*,'Dips, terraces, and points on the fault will not be affected.'
        end if
        print*,'Enter regional dip in degrees. (Negative is down to the left, positive right)'
        read(*,*) RegDip
        print*,'Enter number of beds: '
        read(*,*) nbeds
        allocate(intercepts(nbeds))
        print*,'Enter y intercepts of all ',nbeds,' beds in order'
        read(*,*) intercepts
        write(save_id,*) FitType
        write(save_id,*) RegDip
        write(save_id,*) intercepts
        exit
    else if (FitType == 5) then
        do
            write(save_id,*) FitType
            print*,'Fit for dips? (0/1)'
            read(*,*) fit_dips
            write(save_id,*) fit_dips
            if (fit_dips==0) then
                print*,'Enter regional dip in degrees. (Negative is down to the left, positive right)'
                read(*,*) RegDip
                write(save_id,*) RegDip
                exit
            else if (fit_dips==1) then
                if (DatTypes(DatType_beds)==1) then
                    do
                        print*,'For bed / contact data:'
                        print*,'(1) Use same restored dip for all beds (and dip data).' !Later I might want to allow using a different one for dip data, but I'm not sure when you'd want to do that.
                        print*,'(2) Fit dip separately for each bed.'
                        print*,'(3) Fit dips for groups of beds.'
                        read(*,*) fit_dips
                        write(save_id,*) fit_dips
                        if (fit_dips==1) then
                            exit
                        else if(fit_dips==2) then
                            if (DatTypes(DatType_dips)==1) then
                                print*,'Restored dip for dip data will be fit as a single additional parameter'
                            end if
                            exit
                        else if(fit_dips==3) then
                            print*,'Enter number of groups: '
                            read(*,*) n_bed_groups
                            write(save_id,*) n_bed_groups
                            print*,'Beds are numbered by their order in the data file.'
                            print*,'Note: All beds in a group must be together in the data file.'
                            do i = 1,n_bed_groups
                                print*,'Enter number of first bed in group ',i
                                read(*,*) bed_group_start
                                write(save_id,*) bed_group_start
                            end do
                            if (DatTypes(DatType_dips)==1) then
                                print*,'Restored dip for dip data will be fit as a single additional parameter'
                            end if
                            exit
                        else
                            print*,'Invalid Command'
                        end if
                    end do
                end if
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        print*,'Fit for y intercepts of beds (0/1)?'
        read(*,*) fit_beds
        write(save_id,*) fit_beds
        do
            if(fit_beds==0) then
                exit
            else if(fit_beds==1) then
                exit
            else
                print*,'Invalid Command.'
            end if
        end do
        write(save_id,*) file_name !Since the program will ask for the parameters file name again.
        exit
    else if (FitType == 6) then !Multi-segment restored beds.
        write(save_id,*) FitType
        print*,'Note: Did data will all be fit to the same restored dip'
        print*,'Restored dip for dip data will be fit as a single additional parameter'
        print*,'Enter number of segments per bed: ' !Currently this must be the same for all beds.
        read(*,*) n_bed_segs
        write(save_id,*) n_bed_segs
        !For each bed, parameters will be dips of all n_bed_segs segments in order from left to right, x-coordinates of all n_bed_segs-1 bends from left to right, y-coordinate of the first bend
        write(save_id,*) file_name !Since the program will ask for the parameters file name again.
        exit
    else
        print*,'Invalid Command'
    end if
end do
!Find out if beds are in any particular order
print*,'Are beds in order by age in the data file (as read from top to bottom)?'
print*,'(0) No order'
print*,'(1) Youngest to oldest'
print*,'(2) Oldest to youngest'
read(*,*) beds_age_order
if (beds_age_order == 1) then
    write(save_id,*) 'beds_yng_to_old'
else if (beds_age_order == 2) then
    write(save_id,*) 'beds_old_to_yng'
    beds_age_order = 2
end if
if (beds_age_order==1 .or. beds_age_order==2) then
    !print*,'Warning: The order of restored beds will be measured at x = 0.'
    !print*,'It will not take into account angular unconformities, and will project lower beds above them.'
    print*,'How do you want to determine order?'
    print*,'(1) Use order at x = 0.'
    print*,'(2) Prevent beds from crossing only in the domain of x coordinates spanned by both restored beds.'
    print*,'(3) Prevent beds from crossing anywhere within a specified domain.'
    read(*,*) beds_order_type
    write(save_id,*) 'beds_order_type'
    write(save_id,*) beds_order_type
    if (beds_order_type == 3) then
        print*,'Enter minimum and maximum x coordinates of the domain to consider.'
        read(*,*) bed_domain_x
        write(save_id,*) bed_domain_x
    end if
end if
if (DatTypes(3) == 1) then
    print*,'Enter number of terraces'
    read(*,*) nterr
    print*,'Enter original terrace surface dip at time of formation in degrees for all ',nterr,&
        ' terraces (Negative is down to the left, positive right): '
    allocate(terr_dip_deg(nterr))
    read(*,*) terr_dip_deg
    write(save_id,*) terr_dip_deg
    allocate(terr_orig_elev(nterr))
    print*,'Enter terrace inner edge elevations at time of formation for all ',nterr,' terraces.'
    read(*,*) terr_orig_elev
    write(save_id,*) terr_orig_elev
    print*,'Are terraces in order by age?'
    print*,'(0) No order'
    print*,'(1) Youngest to oldest'
    print*,'(2) Oldest to youngest'
    read(*,*) terr_age_order
    if (terr_age_order==1) then
        write(save_id,*) 'terr_yng_to_old'
    else if (terr_age_order==2) then
        write(save_id,*) 'terr_old_to_yng'
    end if
end if
!Get the type of results to output
do
    print*,'Results To Calculate: '
    print*,'(1) RMS Error'
    print*,'(2) Probability (unnormalized) for uncorrelated data'
    print*,'(3) Chi-square statistic (Cardozo, 2005)'
    if (DatTypes(1) == 1) then !Only for beds
        print*,'(4) Probability (unnormalized) for data correlated along a bed'
        print*,'(5) Uncorrelated and correlated probability along a bed (combination of 1 and 4)'
    end if
    read(*,*) ResultType
    if (ResultType == 1 .or. ResultType == 3) then
        if (method >= 4) then !Markov Chain or RML, requires probability
            print*,'Error: For a Metropolis Algorithm you must choose probability'
        else if (sum(DatTypes)>1 .and. ResultType == 1) then
            print*,'Error: If fitting more than one data type, you must choose probability.'
        else if (ResultType == 3 .and. DatTypes(4) == 1) then
            print*,'Chi squared is not supported for fault points.'
        else
            write(save_id,*) ResultType
            exit
        end if
    else if (ResultType == 2 .or. ResultType == 4 .or. ResultType == 5) then
        if (ResultType == 5) then
            print*,'Warning: Second (correlated) error will not be propagated.'
            print*,'It is recommended not to use error propagation with this option, &
                as it will only appply to the uncorrelated error.'
        end if
        write(save_id,*) ResultType
        do !Find out if uncertainties should be fit for as a model parameter.
            print*,'Should data uncertainties be fit for as a model parameter? (0/1)'
            read(*,*) fit_sigma
            if (fit_sigma==1) then
                write(save_id,*) 'fit_sigma'
                exit
            else if (fit_sigma == 0) then
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do
            print*,'Do you want to propagate errors? (0/1): '
            read(*,*) propagate
            if (propagate == 1 .or. propagate == 0) then
                write(save_id,*) propagate
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        if (DatTypes(DatType_beds) == 1) then
            print*,'Are uncertainties different for different beds?'
            print*,'    (0) All bed uncertainties are the same.'
            print*,'    (1) All beds have independent uncertainties.'
            print*,'    (2) Group of beds have different uncertainties.'
            do                
                read(*,*) diff_sigma_bed
                if (diff_sigma_bed==1) then
                    print*,'Enter number of beds: '
                    read(*,*) nbeds
                    write(save_id,*) 'diff_sigma_bed'
                    allocate(sigma_bed(nbeds),sigma_bed_x(nbeds),sigma_bed_y(nbeds))
                    exit
                else if (diff_sigma_bed==2) then
                    print*,'Enter number of groups: '
                    write(save_id,*) 'diff_sigma_bed_groups'
                    read(*,*) n_sigma_bed_groups
                    write(save_id,*) n_sigma_bed_groups
                    print*,'Beds are numbered by their order in the data file.'
                    print*,'Note: All beds in a group must be together in the data file.'
                    do i = 1,n_sigma_bed_groups
                        print*,'Enter number of first bed in group ',i
                        read(*,*) sigma_bed_group_start
                        write(save_id,*) sigma_bed_group_start
                    end do
                    allocate(sigma_bed(n_sigma_bed_groups),sigma_bed_x(n_sigma_bed_groups),sigma_bed_y(n_sigma_bed_groups))
                    exit
                else if (diff_sigma_bed==0) then
                    allocate(sigma_bed(1),sigma_bed_x(1),sigma_bed_y(1))
                    exit
                else
                    print*,'Invalid Command'
                end if
            end do
            print*,'Are x and y uncertainties different for bed data? (0/1)'
            do
                read(*,*) sigma_xy_diff
                if (sigma_xy_diff == 0) then
                    if (fit_sigma == 0) then
                        if (ResultType==5) print*,'Enter uncorrelated uncertainty at first prompts. &
                            Prompts for correlated uncertainty will come later.'
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Uncertainty in Bed Data for all ',nbeds,' beds: '
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Uncertainty in Bed Data for all ',n_sigma_bed_groups,' bed groups: '
                        else
                            print*,'Enter Uncertainty in Bed Data: '
                        end if
                        read(*,*) sigma_bed
                        write(save_id,*) sigma_bed
                    end if
                    exit
                else if(sigma_xy_diff == 1) then
                    write(save_id,*) 'sigma_xy_diff'
                    if (fit_sigma == 0) then
                        if (ResultType==5) print*,'Enter uncorrelated uncertainty at first prompts. &
                            Prompts for correlated uncertainty will come later.'
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Uncertainty in Bed Data (x position) for all ',nbeds,' beds: '
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Uncertainty in Bed Data (x position) for all ',n_sigma_bed_groups,' bed groups: '
                        else
                            print*,'Enter Uncertainty in Bed Data (x position): '
                        end if
                        read(*,*) sigma_bed_x
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Uncertainty in Bed Data (y position) for all ',nbeds,' beds: '
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Uncertainty in Bed Data (y position) for all ',n_sigma_bed_groups,' bed groups: '
                        else
                            print*,'Enter Uncertainty in Bed Data (y position): '
                        end if
                        read(*,*) sigma_bed_y
                        write(save_id,*) sigma_xy_diff
                        write(save_id,*) sigma_bed_x
                        write(save_id,*) sigma_bed_y
                    end if
                    exit
                else
                    print*,'Invalid Command'
                end if
            end do
            if (ResultType == 4 .or. ResultType == 5) then
                print*,'Note: Only point data (beds) will be affected by correlation length.'
                print*,'    Dips and other data will still be considered uncorrelated.'
                print*,'Do you want to calculate correlation matrix based on (1) deformed state or (2) restored state bed length?'
                do
                    read(*,*) state_for_Lc
                    if (state_for_Lc==1) then
                        write(save_id,*) 'Lc_deformed_state'
                        exit
                    elseif (state_for_Lc==2) then
                        write(save_id,*) 'Lc_restored_state'
                        exit
                    else
                        print*,'Invalid Command'
                    end if
                end do
                print*,'Do you want to fit for the correlation length as a model parameter? (0/1)'
                do
                    read(*,*) fit_lc
                    if (fit_Lc == 1) then
                        write(save_id,*) 'fit_Lc'
                    end if
                    if  (fit_Lc==0 .or. fit_Lc==1) then
                        exit
                    else
                        print*,'Invalid Command'
                    end if
                end do
                if (fit_Lc == 0) then
                    print*,'Enter correlation length for bed data: '
                    read(*,*) Lc
                    write(save_id,*) Lc
                end if
            end if
            if (ResultType==5) then !A second sigma or set of sigmas for beds are needed.
                if (sigma_xy_diff == 0) then
                    if (fit_sigma == 0) then
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Second (Correlated)Uncertainty in Bed Data for all ',nbeds,' beds: '
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data for all ',n_sigma_bed_groups,&
                                ' bed groups: '
                        else
                            print*,'Enter Second (Correlated)Uncertainty in Bed Data: '
                        end if
                        read(*,*) sigma_bed
                        write(save_id,*) sigma_bed
                    end if
                    exit
                else if(sigma_xy_diff == 1) then
                    if (fit_sigma == 0) then
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data (x position) for all ',nbeds,' beds: '
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data (x position) for all ',&
                                n_sigma_bed_groups,' bed groups: '
                        else
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data (x position): '
                        end if
                        read(*,*) sigma_bed_x
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data (y position) for all ',nbeds,' beds: '
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data (y position) for all ',&
                                n_sigma_bed_groups,' bed groups: '
                        else
                            print*,'Enter Second (Correlated) Uncertainty in Bed Data (y position): '
                        end if
                        read(*,*) sigma_bed_y
                        write(save_id,*) sigma_bed_x
                        write(save_id,*) sigma_bed_y
                    end if
                    exit
                else
                    print*,'Invalid Command'
                end if
            end if
        end if
        if (DatTypes(DatType_dips) == 1) then
            if (fit_sigma == 0) then
                print*,'Enter Uncertainty in Dip Data: '
                read(*,*) sigma_dip
                write(save_id,*) sigma_dip
            end if
            do
                print*,'Do you want to include a separate uncertainty in the restored state dips &
                    (which will not be propagated through restoration)? (0/1) '
                read(*,*) use_sigma_restored_dip
                if (use_sigma_restored_dip == 1) then
                    write(save_id,*) 'use_sigma_restored_dip'
                    write(save_id,*) use_sigma_restored_dip
                    if (fit_sigma == 0) then
                        print*,'Enter Uncertainty in Restored State Dip Data: '
                        read(*,*) sigma_restored_dip
                        write(save_id,*) sigma_restored_dip
                    end if
                    exit
                else if (use_sigma_restored_dip == 0) then
                    exit
                else
                    print*,'Invalid Command'
                end if        
            end do    
        end if        
        if (DatTypes(DatType_terr) == 1) then
            if (fit_sigma == 0) then
                print*,'Enter Uncertainties in Terrace Data for all ',nterr,' terraces: '
                allocate(sigma_terr(nterr))
                read(*,*) sigma_terr
                write(save_id,*) sigma_terr
                print*,'Enter Uncertainty in all ',nterr,' Terrace Inner Edge Original Elevations'
                allocate(sigma_orig_elev(nterr))
                read(*,*) sigma_orig_elev
                write(save_id,*) sigma_orig_elev
            end if
        end if
        if (DatTypes(DatType_faultpts) == 1) then
            if (fit_sigma == 0) then
                print*,'Enter Uncertainty in Fault Point Data: '
                read(*,*) sigma_fault
                write(save_id,*) sigma_fault
            end if
        end if
        if (DatTypes(DatType_beds_restored) == 1) then
            !Note: Correlation along the length of a bed is not applied to restored bed points.
            !Also, we are assuming that sigma_xy_diff is the same as for beds (when not restored).
            if (fit_sigma == 0) then
                if (sigma_xy_diff == 0) then
                    print*,'Enter Uncertainty in Restored Bed Data: '
                    read(*,*) sigma_bed_restored
                    write(save_id,*) sigma_bed_restored
                    exit
                else if(sigma_xy_diff == 1) then
                    print*,'Enter Uncertainty in Restored Bed Data (x position): '
                    read(*,*) sigma_bed_restored_x
                    print*,'Enter Uncertainty in Restored Bed Data (y position): '
                    read(*,*) sigma_bed_restored_y
                    write(save_id,*) sigma_bed_restored_x
                    write(save_id,*) sigma_bed_restored_y
                    exit
                else
                    print*,'Invalid Command'
                end if
            end if
        end if
        exit
    else if(ResultType == 3) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
if (fit_sigma==1) then
    write(save_id,*) file_name !Since the program will ask for the parameters file name again.
end if

!Get errors file name
print*,'Input file name to save results to.'
read(*,*) errs_file_name !File name to store the errors in
write(save_id,*) errs_file_name

!Put a blank line at the end
write(save_id,*)

close(save_id)

end program Setup_Trishear