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

module user_input
!Assignments and subprograms related to reading files or getting user input.

implicit none

contains

subroutine GetOptionsID
use options, only: options_id
implicit none
integer :: choice !1 or 2 to specify which input option is chosen
character (len=50) :: file_name !file name where the options come from if choice == 2
logical :: ex !tells whether the file named exists
do
print*, 'Choose source of run options (not parameter space):'
print*, '(1) Input manually'
print*, '(2) Read from file'
read(*,*) choice
if (choice == 1) then
    options_id = 5
    exit
elseif (choice == 2) then
    options_id = 4
    do !Do statement for reading options id file name
        print*,'Input file name (with extension)'
        read(*,*) file_name
        inquire(file=file_name, exist=ex) !Check that the file exists
        if (ex .eqv. .true.) then
            exit
        else
            print*,'Invalid File Name'
        end if
    end do
    open(unit=4, file = file_name)
    exit
end if
end do
end subroutine GetOptionsID

subroutine GetDatType ! Find out data type
use options, only: options_id,DatTypes,ndatatypes,DatType_beds,DatType_dips,DatType_terr,DatType_faultpts,DatType_beds_restored
implicit none
integer,dimension(:),allocatable :: datalist
integer :: i,goodinput
print*,'Enter number of different data types: '
read(options_id,*) ndatatypes
allocate(datalist(ndatatypes))
do 
    print*,'Data Are (list all that apply): '
    print*,'(1) beds '
    print*,'(2) dips '
    print*,'(3) points on a fault'
    print*,'(4) marine terraces'
    print*,'(5) points on restored state beds'
    read(options_id,*) datalist
    DatTypes = 0 !Initialize them all to zero
    goodinput = 1 !Initialize to 1
    do i=1,ndatatypes
        select case(datalist(i)) !I could replace all this with DatTypes(datalist(i)) = 1, except I need to know if a value is not good.
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
end subroutine GetDatType

subroutine ReadBeds !Read in a bedding data file
use options, only: options_id
use data_module, only: bed,beds,nbeds,ngrowth
use growth_strata, only: GetGrowthParameterization
implicit none
character (len=50) :: file_name
logical :: ex !tells whether the file named exists
double precision :: x,y !Temporary variables for the x and y coordinates of each point
character (len=20) :: strx,stry !x and y as strings
character (len=20) :: bedname !the name of a bed
logical :: isgrowth !Tells whether the bed is growth strata or not.
double precision,dimension(:,:),allocatable :: bedxy    !combined x and y coordinates of the bed being read
type(bed) :: newbed !the bed data type containing all the data about the bed being read
integer :: IOstatus !status of reading
integer,dimension(:),allocatable :: bed_lengths
integer :: n,i !counters of beds and points within bed respectively
!Open file
do !Do statement for reading bed file name
    print*,'Input bedding file name (with extension)'
    read(options_id,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
open(unit=1, file = file_name)
!Read in the bed information
nbeds = 0 !Initialize with no beds
ngrowth = 0 !Initialize with no growth strata.
allocate(bed_lengths(0)) !Initialize empty
do
    read(1,*,IOSTAT = IOstatus) strx,stry
    if (IOstatus >0) then
        print*,'Problem Encountered in file.'
    else if (IOstatus < 0) then
        !End of file reached
        exit
    else
        if (strx == 'bed') then
            !New bed
            nbeds = nbeds+1
            bed_lengths = (/bed_lengths,0/)
        else if(strx == 'growth') then
            !New growth strata bed
            nbeds = nbeds+1
            ngrowth = ngrowth+1
            bed_lengths = (/bed_lengths,0/)
        else
            !Still in the same bed
            bed_lengths(nbeds) = bed_lengths(nbeds)+1
        end if
    end if
end do
rewind(1)
allocate(beds(nbeds))
print*,'Reading ',nbeds,' beds'
print*,'Name:        number of points:' !Column headers for what will be printed out.
do n = 1,nbeds
    allocate(bedxy(2,bed_lengths(n)))
    do i = 1,bed_lengths(n)+1
        if (i == 1) then
            read(1,*) strx,bedname
            if (strx == 'growth') then
                isgrowth = .true.
            else
                isgrowth = .false.
            end if
        else
            read(1,*) x,y
            bedxy(1,i-1) = x
            bedxy(2,i-1) = y
        end if
    end do
    newbed = bed(bedxy,bed_lengths(n),bedname,isgrowth)
    beds(n) = newbed
    deallocate(bedxy)
    print*,bedname,bed_lengths(n)
end do
close(1)
if (ngrowth > 0.0) then
    call GetGrowthParameterization
end if

end subroutine ReadBeds

subroutine ReadDips
use constants
use data_module
use options, only: options_id
implicit none
character (len=50) :: file_name
logical :: ex !tells whether the file named exists
integer :: lines !number of lines
integer :: IOstatus !status of reading
integer :: i !Counter
double precision :: dipx,dipy,dip !Temporary dip position and values
!Open file
do !Do statement for reading dip file name
    print*,'Input dip file name (with extension)'
    read(options_id,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
open(unit=1, file = file_name)
!Find out number of lines
lines = 0
do
    read(1,*,IOSTAT = IOstatus)
    if (IOstatus >0) then
        print*,'Problem Encountered in file on line',lines+1
    else if (IOstatus < 0) then
        exit
    else
        lines = lines+1
    end if
end do
!Then allocate arrays and read file
allocate(dip_pos(2,lines),dips(lines))
rewind(1)
do i = 1,lines
    read(1,*) dipx,dipy,dip !read in x position, y position, and the dip
    dip_pos(1,i) = dipx
    dip_pos(2,i) = dipy
    dips(i) = dip*pi/180.
end do
ndips = size(dips,1) !number of dip measurements
print*,'Read ',ndips,' dips'
close(1)
end subroutine ReadDips

subroutine ReadTerraces !Read in a terrace data file
use data_module
use options, only: options_id
implicit none
character (len=50) :: file_name
logical :: ex !tells whether the file named exists
double precision :: x,y !Temporary variables for the x and y coordinates of each point
character (len=20) :: strx,stry !x and y as strings
character (len=20) :: terrname !the name of a terrace
double precision,dimension(:,:),allocatable :: terrxy    !combined x and y coordinates of the terrace being read
type(bed) :: newterr !the bed data type containing all the data about the terrace being read
integer :: IOstatus !status of reading
integer,dimension(:),allocatable :: terr_lengths
integer :: n,i !counters of terraces and points within a terrace respectively
!Open file
do !Do statement for reading bed file name
    print*,'Input terrace file name (with extension)'
    read(options_id,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
open(unit=1, file = file_name)
print*,'Note: First point in terrace file is considered to be inner edge elevation.'
!Read in the bed information
nterr = 0 !Initialize with 0
allocate(terr_lengths(0)) !Initialize empty
do
    read(1,*,IOSTAT = IOstatus) strx,stry
    if (IOstatus >0) then
        print*,'Problem Encountered in file.'
    else if (IOstatus < 0) then
        !End of file reached
        exit
    else
        if (strx == 'terrace') then
            !New terrace
            nterr = nterr+1
            terr_lengths = (/terr_lengths,0/)
        else
            !Still in the same terrace
            terr_lengths(nterr) = terr_lengths(nterr)+1
        end if
    end if
end do
rewind(1)
allocate(terraces(nterr))
do n = 1,nterr
    allocate(terrxy(2,terr_lengths(n)))
    do i = 1,terr_lengths(n)+1
        if (i == 1) then
            read(1,*) strx,terrname
        else
            read(1,*) x,y
            terrxy(1,i-1) = x
            terrxy(2,i-1) = y
        end if
    end do
    newterr = bed(terrxy,terr_lengths(n),terrname,.false.)
    terraces(n) = newterr
    deallocate(terrxy)
end do
close(1)
print*,'Read ',nterr,' terraces'
end subroutine ReadTerraces

subroutine ReadFaultPts !Read in a file containing points on the fault.
use data_module
use options, only: options_id
implicit none
character (len=50) :: file_name
logical :: ex !tells whether the file named exists
double precision :: x,y !Temporary variables for the x and y coordinates of each point
character (len=20) :: strx,stry !x and y as strings
integer :: IOstatus !status of reading
integer :: i !counters of terraces and points within a terrace respectively
!Open file
do !Do statement for reading bed file name
    print*,'Input fault file name (with extension)'
    read(options_id,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
open(unit=1, file = file_name)
!Read in the fault information
do
    read(1,*,IOSTAT = IOstatus) strx,stry
    if (IOstatus >0) then
        print*,'Problem Encountered in file.'
    else if (IOstatus < 0) then
        !End of file reached
        exit
    else
        nfaultpts = nfaultpts+1
    end if
end do
rewind(1)
allocate(faultpts(2,nfaultpts))
do i = 1,nfaultpts
    read(1,*) x,y
    faultpts(1,i) = x
    faultpts(2,i) = y
end do
close(1)
print*,'Read ',nfaultpts,' points on the fault'
end subroutine ReadFaultPts

subroutine ReadRestoredBedPts
use options, only: options_id
use data_module, only: bed,restored_beds,n_restored_beds
implicit none
character (len=50) :: file_name
logical :: ex !tells whether the file named exists
double precision :: x,y !Temporary variables for the x and y coordinates of each point
character (len=20) :: strx,stry !x and y as strings
character (len=20) :: bedname !the name of a bed
logical :: isgrowth !Tells whether the bed is growth strata or not.
double precision,dimension(:,:),allocatable :: bedxy    !combined x and y coordinates of the bed being read
type(bed) :: newbed !the bed data type containing all the data about the bed being read
integer :: IOstatus !status of reading
integer,dimension(:),allocatable :: bed_lengths
integer :: n,i !counters of beds and points within bed respectively
!Open file
do !Do statement for reading bed file name
    print*,'Input restored bedding file name (with extension)'
    read(options_id,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
open(unit=1, file = file_name)
!Read in the bed information
n_restored_beds = 0 !Initialize with no beds
allocate(bed_lengths(0)) !Initialize empty
do
    read(1,*,IOSTAT = IOstatus) strx,stry
    if (IOstatus >0) then
        print*,'Problem Encountered in file.'
    else if (IOstatus < 0) then
        !End of file reached
        exit
    else
        if (strx == 'bed' .or. strx == 'growth') then !For this purpose, it doesn't matter if it's a growth or pregrowth bed, but either title can be used.
            !New bed
            n_restored_beds = n_restored_beds+1
            bed_lengths = (/bed_lengths,0/)
        else
            !Still in the same bed
            bed_lengths(n_restored_beds) = bed_lengths(n_restored_beds)+1
        end if
    end if
end do
rewind(1)
allocate(restored_beds(n_restored_beds))
print*,'Reading ',n_restored_beds,' restored beds'
print*,'Name:        number of points:' !Column headers for what will be printed out.
do n = 1,n_restored_beds
    allocate(bedxy(2,bed_lengths(n)))
    do i = 1,bed_lengths(n)+1
        if (i == 1) then
        read(1,*) strx,bedname
            if (strx == 'growth') then
                isgrowth = .true.
            else
                isgrowth = .false.
            end if
        else
            read(1,*) x,y
            bedxy(1,i-1) = x
            bedxy(2,i-1) = y
        end if
    end do
    newbed = bed(bedxy,bed_lengths(n),bedname,isgrowth)
    restored_beds(n) = newbed
    deallocate(bedxy)
    print*,bedname,bed_lengths(n)
end do
close(1)
end subroutine ReadRestoredBedPts

subroutine GetFaultType
use options, only: options_id,FaultType
use options,only:TipToSolve
implicit none
do
    print*,'Choose Fault Type'
    print*,'(1) Straight Fault Ramp'
    print*,'(2) Ramp from Horizontal Detachment'
    print*,'(3) Fault with a bend in it'
    print*,'(4) Listric Fault (fault parallel flow)'
    print*,'(5) Parallel Fault Propagation Fold'
    print*,'(6) Multi-bend Fault'
    print*,'(7) Listric Fault (circular or elliptic, approximate)'
    print*,'(8) Spline Fault (Requires tip position is final)'
    read(options_id,*) FaultType
    if (1<=FaultType .and. FaultType<=7) then
        exit
    else if (FaultType==8 .and. TipToSolve==2) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
end subroutine GetFaultType

subroutine ReadParams
use parameters,only: nparams,mins,maxs,steps,slip_sense,increment,v0
use options
implicit none
character (len=50) :: file_name
logical :: ex !tells whether the file named exists
integer :: i !a counter
real :: slip_sense_real !This is a terrible, cobbled together solution for the fact that if we have to call read params twice, the number that first gets read for this variable may not be an integer.
!Get the file name and open it
do !Do statement for reading parameter file name
    print*,'Input Parameter File Name (with extension)'
    read(options_id,*) file_name
    inquire(file=file_name, exist=ex) !Check that the file exists
    if (ex .eqv. .true.) then
        exit
    else
        print*,'Invalid File Name'
    end if
end do
open(unit=params_id,file=file_name)
!Read the parameters
if (allocated(mins) .and. allocated(maxs) .and. allocated(steps)) then !Occurs when we call ReadParams again when fitting for bed dip / restored depth.
    deallocate(mins,maxs,steps)
end if
allocate(mins(nparams),maxs(nparams),steps(nparams))
do i = 1,nparams
    read(params_id,*) mins(i),maxs(i),steps(i)
end do
if (maxs(5) == 90) then !i.e. if phi = 90 deg
    print*,'Error: phi = 90 degrees, not allowed, using 89 degrees'
end if
read(params_id,*) slip_sense_real
slip_sense = -nint(slip_sense_real) !Since this is inversion, make it the opposite of what it was
read(params_id,*) increment
close(params_id) !Close the file
increment = increment*slip_sense !Make it increment in the right direction
v0 = increment ! phrases it as a velocity
end subroutine ReadParams

subroutine GetMethod
use options, only: options_id,method
implicit none
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
    read(options_id,*) method
    if (method>0 .and. method<8) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
end subroutine GetMethod


subroutine GetFitType !Find out the type of fit to do
use options, only: FitType,DatTypes,options_id,DatType_beds,DatType_dips,fit_dips,fit_beds,beds_age_order,n_bed_segs, &
    n_bed_groups,bed_group_start
use parameters, only: RegDip,RegDipCos,RegDipSin,RegSlope,intercepts,nparams,nparam_start_restored_fit,nrestored_fit_params
use data_module, only: nbeds,beds
use beds_order, only: GetBedsOrderType
use options, only: options_id,limit_bed_dips,min_bed_slopes,max_bed_slopes
use constants
implicit none
character (len = 25) :: str !A string to tell whether beds are in a particular age order.
integer :: i !A counter
double precision,dimension(:),allocatable :: min_bed_dips,max_bed_dips !If limit_bed_dips is 1, these give the minimum and maximum allowed dip for each bed.
do
    print*,'Objective Function for beds and/or dips:'
    print*,'(1) Fit to flat line.'
    print*,'(2) Fit to known dip.'
    print*,'(3) Fit to best fit line / mean dip.'
    if (DatTypes(DatType_beds) == 1) then !only useful if beds are included
        print*,'(4) Fit to known line'
    end if
    if (DatTypes(DatType_beds)==1 .or. DatTypes(DatType_dips)==1) then
        print*,'(5) Fit for restored dips and/or bed elevations as model parameters'
    end if
    if (DatTypes(DatType_beds)==1) then
        print*,'(6) Multi-segment restored beds'
    end if
    read(options_id,*) FitType
    if (FitType == 1) then
        print*,'Note: This method fits to an average (unweighted) of the restored bed points'
        RegDip = 0
        RegSlope = 0
        RegDipSin = 0
        RegDipCos = 1
        exit
    else if (FitType == 2) then
        print*,'Note: This method fits to an average (unweighted) of the restored bed intercepts corresponding to the given slope'
        call GetRegDip
        exit
    else if (FitType == 3) then
        print*,'Note: This method uses an unweighted least squares fit'
        if (DatTypes(DatType_beds) == 1 .and. DatTypes(DatType_dips) == 1) then !Data are both beds and dips
            print*,'Warning: Bed and dip data will be fit separately, not necessarily to the same dip.'
        end if
        print*,'Do you want to put limits on restored bed dips? (0=no, 1=yes)'
        if (options_id == 5) then !Manual input
            read(options_id,*) limit_bed_dips
        else
            read(options_id,*) str
            if (str == 'limit_bed_dips') then                    
                limit_bed_dips = 1
            else
                backspace(options_id)
                limit_bed_dips = 0
            end if                
        end if
        if (limit_bed_dips==1) then
            allocate(min_bed_dips(nbeds),max_bed_dips(nbeds),min_bed_slopes(nbeds),max_bed_slopes(nbeds))
            print*,'Note: Bed dips are negative down to left and positive down to right.'
            print*,'Enter minimum dips for all ',nbeds,' beds:'
            read(options_id,*) min_bed_dips
            min_bed_slopes = tan(min_bed_dips*pi/180)
            print*,'Enter maximum dips for all ',nbeds,' beds:'
            read(options_id,*) max_bed_dips
            max_bed_slopes = tan(max_bed_dips*pi/180)
        end if
        exit
    else if (FitType == 4) then
        if (sum(DatTypes) > 1) then
            print*,'Note: Only beds will be fit to specified lines.'
            print*,'Dips, terraces, and points on the fault will not be affected.'
        end if
        call GetRegDip
        allocate(intercepts(nbeds))
        print*,'Enter y intercepts of all ',nbeds,' beds in order'
        read(options_id,*) intercepts
        exit
    else if (FitType == 5) then
        nrestored_fit_params = 0 !Number of parameters to be fit for. it will increase.
        nparam_start_restored_fit = nparams+1 !First parameter number for the parameters used for the restored bed / dip values.
        do
            print*,'Fit for dips? (0/1)'
            read(options_id,*) fit_dips
            if (fit_dips==0) then
                call GetRegDip
                exit
            else if (fit_dips==1) then
                if (DatTypes(DatType_beds)==1) then
                    do
                        print*,'For bed / contact data:'
                        print*,'(1) Use same restored dip for all beds (and dip data).' !Later I might want to allow using a different one for dip data, but I'm not sure when you'd want to do that.
                        print*,'(2) Fit dip separately for each bed.'
                        print*,'(3) Fit dips for groups of beds.'
                        read(options_id,*) fit_dips
                        if (fit_dips==1) then
                            nrestored_fit_params = nrestored_fit_params+1
                            exit
                        else if(fit_dips==2) then
                            nrestored_fit_params = nrestored_fit_params+nbeds
                            if (DatTypes(DatType_dips)==1) then
                                print*,'Restored dip for dip data will be fit as a single additional parameter' !Maybe I should make it possible to include dips in one of the groups.
                                nrestored_fit_params = nrestored_fit_params+1
                            end if
                            exit
                        else if(fit_dips==3) then
                            print*,'Enter number of groups: '
                            read(options_id,*) n_bed_groups
                            nrestored_fit_params = nrestored_fit_params+n_bed_groups
                            print*,'Beds are numbered by their order in the data file as follows:'
                            do i = 1,nbeds
                                print*,i,': ',beds(i)%ident
                            end do
                            print*,'Note: All beds in a group must be together in the data file.'
                            allocate(bed_group_start(n_bed_groups))
                            do i = 1,n_bed_groups
                                print*,'Enter number of first bed in group ',i
                                read(options_id,*) bed_group_start(i)
                            end do
                            if (DatTypes(DatType_dips)==1) then
                                print*,'Restored dip for dip data will be fit as a single additional parameter'
                                nrestored_fit_params = nrestored_fit_params+1
                            end if
                            exit
                        else
                            print*,'Invalid Command'
                        end if
                    end do
                else
                    nrestored_fit_params = nrestored_fit_params+1
                end if
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        print*,'Fit for y intercepts of beds (0/1)?'
        read(options_id,*) fit_beds
        do
            if(fit_beds==0) then
                exit
            else if(fit_beds==1) then
                nrestored_fit_params = nrestored_fit_params+nbeds
                exit
            else
                print*,'Invalid Command.'
            end if
        end do
        if (fit_dips==0 .and. fit_beds==0) then !I'm not sure why any one would do this, but if you choose these options it becomes the same as FitType = 2, so you might as well treat it as such.
            FitType=2
        end if
        nparams = nparams+nrestored_fit_params
        call ReadParams !ReadParams is called before GetFitType. We've added extra parameters, so we need to reread them.
        exit
    else if (FitType == 6) then !Multi-segment restored beds.
        print*,'Note: Did data will all be fit to the same restored dip'
        print*,'Restored dip for dip data will be fit as a single additional parameter'
        nparam_start_restored_fit = nparams+1 !First parameter number for the parameters used for the restored bed / dip values.
        print*,'Enter number of segments per bed: ' !Currently this must be the same for all beds.
        read(options_id,*) n_bed_segs
        !For each bed, parameters will be dips of all n_bed_segs segments in order from left to right, x-coordinates of all n_bed_segs-1 bends from left to right, y-coordinate of the first bend
        !In the parameters file, these are all together in order for each bed, then the next bed, and so on. The restored dip for dip data comes after all of these.
        nrestored_fit_params = 2*n_bed_segs*nbeds
        !print*,'Must restored beds be parallel (0/1)?'
        !!If yes (1), then the order of parameters will be n dips for all n_bed_segs segments, x-coordinates of all n_bed_segs-1 bends on the first bed, axis dip of all n_bed_segs-1 fold axes, y-coordinates of the first bend for all beds.
        if (DatTypes(DatType_dips)==1) then
            print*,'Restored dip for dip data will be fit as a single additional parameter'
            nrestored_fit_params = nrestored_fit_params+1
        end if
        nparams = nparams+nrestored_fit_params
        call ReadParams !ReadParams is called before GetFitType. We've added extra parameters, so we need to reread them.
        exit
    else
        print*,'Invalid Command'
    end if
end do
print*,'Are beds in order by age in the data file (as read from top to bottom)?'
print*,'(0) No order'
print*,'(1) Youngest to oldest'
print*,'(2) Oldest to youngest'
if (options_id == 5) then !Manual input
    read(options_id,*) beds_age_order
else
    read(options_id,*) str
    if (str == 'beds_yng_to_old') then                    
        beds_age_order = 1
    else if (str == 'beds_old_to_yng') then
        beds_age_order = 2
    else
        backspace(options_id)
        beds_age_order = 0
    end if                
end if
if (beds_age_order==1 .or. beds_age_order==2) then
    !print*,'Warning: The order of restored beds will be measured at x = 0.'
    !print*,'It will not take into account angular unconformities, and will project lower beds above them.'
    call GetBedsOrderType
end if
end subroutine GetFitType

subroutine GetTerrFitInfo !Get information necessary for fitting terrace data.
use parameters, only: terr_dip,terr_slope,terr_dip_deg,terr_orig_elev
!use parameters, only: terr_slope,terr_orig_elev
use data_module, only: nterr
use options, only: options_id,terr_age_order
use constants
implicit none
!double precision,dimension(:),allocatable :: terr_dip_deg,terr_dip !Dip at which terrace surfaces are formed.
character (len = 25) :: str !A string to tell whether terraces are in a particular age order.
print*,'Enter original terrace surface dip at time of formation in degrees for all ',nterr,&
    ' terraces (Negative is down to the left, positive right): '
allocate(terr_dip_deg(nterr),terr_dip(nterr),terr_slope(nterr))
read(options_id,*) terr_dip_deg
terr_dip = terr_dip_deg*(pi/180.)
terr_slope = tan(terr_dip) !Note that sign convention is opposite of regular usage.
allocate(terr_orig_elev(nterr))
print*,'Enter terrace inner edge elevations at time of formation for all ',nterr,' terraces.'
read(options_id,*) terr_orig_elev
print*,'Are terraces in order by age?'
print*,'(0) No order'
print*,'(1) Youngest to oldest'
print*,'(2) Oldest to youngest'
if (options_id == 5) then !Manual input
    read(options_id,*) terr_age_order
else
    read(options_id,*) str
    if (str == 'terr_yng_to_old') then                    
        terr_age_order = 1
    else if (str == 'terr_old_to_yng') then
        terr_age_order = 2
    else
        backspace(options_id)
        terr_age_order = 0
    end if                
end if
end subroutine GetTerrFitInfo

subroutine GetRegDip !Find out the regional dip
use parameters, only: RegDip,RegDipDeg,RegDipCos,RegDipSin,RegSlope
use options, only: options_id
use constants
implicit none
print*,'Enter regional dip in degrees. (Negative is down to the left, positive right)'
read(options_id,*) RegDipDeg
RegDip = RegDipDeg*pi/180
RegDipSin = sin(RegDip)
RegDipCos = cos(RegDip)
RegSlope = tan(RegDip) !Express it as a slope instead of an angle Note: sign convention is backwards of regular usage
end subroutine GetRegDip

subroutine GetResultType !Find out what kind of result to output
use options, only: ResultType, method,options_id,DatTypes,propagate,DatType_beds,DatType_faultpts,DatType_dips,DatType_terr, &
    DatType_beds_restored,use_sigma_restored_dip,sigma_xy_diff,diff_sigma_bed,state_for_Lc,fit_sigma,n_sigma_bed_groups, &
    sigma_bed_group_start,fit_lc
use parameters, only:  nparams
use data_uncertainties, only: Lc_fixed,sigma_bed_x_fixed,sigma_bed_y_fixed,sigma_dip_fixed,sigma_restored_dip_fixed,&
        sigma_terr_fixed,sigma_orig_elev_fixed,sigma_fault_fixed,sigma_bed_restored_x_fixed,sigma_bed_restored_y_fixed,&
        nparam_start_sigma_bed,nparam_start_sigma_dip,nparam_start_sigma_terr,nparam_start_sigma_fault,&
        nparam_start_sigma_bed_restored,nparam_Lc,sigma2_bed,sigma2_bed_x_fixed,sigma2_bed_y_fixed,&
        nparam_start_sigma2_bed
use data_module, only:nterr,nbeds,beds,spher_var_R
use err_and_prob, only: make_spherical_variogram_R
implicit none
character (len = 25) :: str !A string to tell whether to use an error in the restored state dips.
double precision :: sigma_bed,sigma_bed_restored !A variable into which to read a single bed before assigning it to all the elements of a sigma_bed... array.
integer :: i
do
    print*,'Results To Calculate: '
    print*,'(1) RMS Error'
    print*,'(2) Probability (unnormalized) for uncorrelated data'
    print*,'(3) Chi-square statistic (Cardozo, 2005)'
    if (DatTypes(DatType_beds) == 1) then !Only for beds
        print*,'(4) Probability (unnormalized) for data correlated along a bed'
        print*,'(5) Uncorrelated and correlated probability along a bed (combination of 1 and 4)'
    end if
    read(options_id,*) ResultType
    if (ResultType == 1 .or. ResultType == 3) then
        if (method >= 4) then !Markov Chain or RML, requires probability
            print*,'Error: For a Metropolis Algorithm you must choose probability'
        else if (sum(DatTypes)>1 .and. ResultType == 1) then
            print*,'Error: If fitting more than one data type, you must choose probability or Chi-square.'
        else if (ResultType == 3 .and. DatTypes(DatType_faultpts) == 1) then
            print*,'Chi squared is not supported for fault points.'
        else
            exit
        end if
    else if (ResultType == 2 .or. ResultType == 4 .or. ResultType == 5) then
        if (ResultType == 5) then
            print*,'Warning: Second (correlated) error will not be propagated.'
            print*,'It is recommended not to use error propagation with this option, &
                as it will only appply to the uncorrelated error.'
        end if
        do !Find out if uncertainties should be fit for as a model parameter.
            print*,'Should data uncertainties be fit for as a model parameter? (0/1)'
            if (options_id == 5) then !Manual input
                read(options_id,*) fit_sigma
            else
                read(options_id,*) str
                if (str == 'fit_sigma') then                    
                    fit_sigma = 1
                else
                    backspace(options_id)
                    fit_sigma = 0
                end if                
            end if
            if (fit_sigma==1 .or. fit_sigma == 0) then
                exit
            else
                print*,'Invalid Command'
            end if
        end do
        do
            print*,'Do you want to propagate errors? (0/1): '
            read(options_id,*) propagate
            if (propagate == 1 .or. propagate == 0) then
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
                if (options_id == 5) then !Manual input
                    read(options_id,*) diff_sigma_bed
                else
                    read(options_id,*) str
                    if (str == 'diff_sigma_bed') then                    
                        diff_sigma_bed = 1
                    else if (str == 'diff_sigma_bed_groups') then
                        diff_sigma_bed = 2
                        print*,'Enter number of groups: '
                        read(options_id,*) n_sigma_bed_groups
                        print*,'Beds are numbered by their order in the data file as follows:'
                        do i = 1,nbeds
                            print*,i,': ',beds(i)%ident
                        end do
                        print*,'Note: All beds in a group must be together in the data file.'
                        allocate(sigma_bed_group_start(n_sigma_bed_groups))
                        do i = 1,n_sigma_bed_groups
                            print*,'Enter number of first bed in group ',i
                            read(options_id,*) sigma_bed_group_start(i)
                        end do
                    else
                        backspace(options_id)
                        diff_sigma_bed = 0
                    end if                
                end if
                if (diff_sigma_bed==1 .or. diff_sigma_bed==0 .or. diff_sigma_bed==2) then
                    allocate(sigma_bed_x_fixed(nbeds),sigma_bed_y_fixed(nbeds))
                    exit
                else
                    print*,'Invalid Command'
                end if
            end do
            print*,'Are x and y uncertainties different for bed data? (0/1)'
            do
                if (options_id == 5) then !Manual input
                    read(options_id,*) sigma_xy_diff
                else
                    read(options_id,*) str
                    if (str == 'sigma_xy_diff') then                    
                        read(options_id,*) sigma_xy_diff
                    else
                        backspace(options_id)
                        sigma_xy_diff = 0
                    end if                
                end if
                if (fit_sigma == 0) then
                    if (ResultType==5) print*,'Enter uncorrelated uncertainty at first prompts. &
                        Prompts for correlated uncertainty will come later.'
                    if (sigma_xy_diff == 0) then
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Uncertainty in Bed Data for all ',nbeds,' beds: '
                            read(options_id,*) sigma_bed_x_fixed
                            sigma_bed_y_fixed = sigma_bed_x_fixed
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Uncertainty in Bed Data for all ',n_sigma_bed_groups,' bed groups: '
                            read(options_id,*) sigma_bed_x_fixed
                            sigma_bed_y_fixed = sigma_bed_x_fixed
                        else
                            print*,'Enter Uncertainty in Bed Data: '
                            read(options_id,*) sigma_bed
                            sigma_bed_x_fixed = sigma_bed
                            sigma_bed_y_fixed = sigma_bed
                        end if
                        exit
                    else if(sigma_xy_diff == 1) then
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Uncertainty in Bed Data (x position) for all ',nbeds,' beds: '
                            read(options_id,*) sigma_bed_x_fixed
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Uncertainty in Bed Data (x position) for all ',n_sigma_bed_groups,' bed groups: '
                            read(options_id,*) sigma_bed_x_fixed
                        else
                            print*,'Enter Uncertainty in Bed Data (x position): '
                            read(options_id,*) sigma_bed
                            sigma_bed_x_fixed = sigma_bed
                        end if
                        if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                            print*,'Enter Uncertainty in Bed Data (y position) for all ',nbeds,' beds: '
                            read(options_id,*) sigma_bed_y_fixed
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            print*,'Enter Uncertainty in Bed Data (y position) for all ',n_sigma_bed_groups,' bed groups: '
                            read(options_id,*) sigma_bed_y_fixed
                        else
                            print*,'Enter Uncertainty in Bed Data (y position): '
                            read(options_id,*) sigma_bed
                            sigma_bed_y_fixed = sigma_bed
                        end if
                        exit
                    else
                        print*,'Invalid Command'
                    end if
                else !Fitting for sigma.
                    nparam_start_sigma_bed = nparams+1
                    if (sigma_xy_diff ==0) then
                        if (diff_sigma_bed == 1) then !Different uncertainties for each bed.
                            nparams = nparams+nbeds
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            nparams = nparams+n_sigma_bed_groups
                        else
                            nparams = nparams+1
                        end if
                        exit
                    else if(sigma_xy_diff == 1) then
                        if (diff_sigma_bed == 1) then  !Different uncertainties for each bed.
                            nparams = nparams+nbeds*2
                        else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                            nparams = nparams+n_sigma_bed_groups*2
                        else
                            nparams = nparams+2
                        end if
                    else
                        print*,'Invalid Command'
                    end if
                end if
            end do
            if (ResultType == 4 .or. ResultType == 5) then
                print*,'Note: Only point data (beds) will be affected by correlation length.'
                print*,'    Dips and other data will still be considered uncorrelated.'
                print*,'Do you want to calculate correlation matrix based on (1) deformed state or (2) restored state bed length?'
                do
                    if (options_id == 5) then !Manual input
                        read(options_id,*) state_for_Lc
                    else
                        read(options_id,*) str
                        if (str == 'Lc_deformed_state') then                    
                            state_for_Lc = 1
                        else if(str == 'Lc_restored_state') then
                            state_for_Lc = 2
                        else
                            backspace(options_id)
                            state_for_Lc = 2
                        end if                
                    end if
                    if  (state_for_Lc==1 .or. state_for_Lc==2) then
                        exit
                    else
                        print*,'Invalid Command'
                    end if
                end do
                print*,'Do you want to fit for the correlation length as a model parameter? (0/1)'
                do
                    if (options_id == 5) then  !Manual input
                        read(options_id,*) fit_Lc
                    else
                        read(options_id,*) str
                        if (str=='fit_Lc') then
                            fit_Lc = 1
                        else
                            backspace(options_id)
                            fit_Lc = 0
                        end if
                    end if
                    if  (fit_Lc==0 .or. fit_Lc==1) then
                        exit
                    else
                        print*,'Invalid Command'
                    end if
                end do
                if (fit_Lc == 0) then
                    print*,'Enter correlation length for bed data: '
                    read(options_id,*) Lc_fixed
                    print*,Lc_fixed
                    if (state_for_Lc==1) then !Lc is in the deformed state.
                        call make_spherical_variogram_R(beds,Lc_fixed,spher_var_R)
                    end if
                else
                    nparams = nparams+1
                    nparam_lc = nparams
                end if
            end if
            if (ResultType==5) then !A second sigma or set of sigmas for beds are needed.
                if (fit_sigma == 0) then
                        if (sigma_xy_diff == 0) then
                            if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data for all ',nbeds,' beds: '
                                read(options_id,*) sigma2_bed_x_fixed
                                sigma_bed_y_fixed = sigma2_bed_x_fixed
                            else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data for all ',n_sigma_bed_groups,&
                                    ' bed groups: '
                                read(options_id,*) sigma2_bed_x_fixed
                                sigma_bed_y_fixed = sigma2_bed_x_fixed
                            else
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data: '
                                read(options_id,*) sigma2_bed
                                sigma_bed_x_fixed = sigma2_bed
                                sigma_bed_y_fixed = sigma2_bed
                            end if
                            exit
                        else if(sigma_xy_diff == 1) then
                            if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data (x position) for all ',nbeds,' beds: '
                                read(options_id,*) sigma2_bed_x_fixed
                            else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data (x position) for all ',&
                                    n_sigma_bed_groups,' bed groups: '
                                read(options_id,*) sigma2_bed_x_fixed
                            else
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data (x position): '
                                read(options_id,*) sigma2_bed
                                sigma_bed_x_fixed = sigma2_bed
                            end if
                            if (diff_sigma_bed==1) then !Different uncertainties for each bed.
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data (y position) for all ',nbeds,' beds: '
                                read(options_id,*) sigma2_bed_y_fixed
                            else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data (y position) for all ',&
                                    n_sigma_bed_groups,' bed groups: '
                                read(options_id,*) sigma2_bed_y_fixed
                            else
                                print*,'Enter Second (Correlated) Uncertainty in Bed Data (y position): '
                                read(options_id,*) sigma2_bed
                                sigma_bed_y_fixed = sigma2_bed
                            end if
                            exit
                        else
                            print*,'Invalid Command'
                        end if
                    else !Fitting for sigma.
                        nparam_start_sigma2_bed = nparams+1
                        if (sigma_xy_diff ==0) then
                            if (diff_sigma_bed == 1) then !Different uncertainties for each bed.
                                nparams = nparams+nbeds
                            else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                                nparams = nparams+n_sigma_bed_groups
                            else
                                nparams = nparams+1
                            end if
                            exit
                        else if(sigma_xy_diff == 1) then
                            if (diff_sigma_bed == 1) then  !Different uncertainties for each bed.
                                nparams = nparams+nbeds*2
                            else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                                nparams = nparams+n_sigma_bed_groups*2
                            else
                                nparams = nparams+2
                            end if
                        else
                            print*,'Invalid Command'
                        end if
                    end if
            end if
        end if
        if (DatTypes(DatType_dips) == 1) then
            if (fit_sigma == 0) then
                print*,'Enter Uncertainty in Dip Data: '
                read(options_id,*) sigma_dip_fixed
            else
                nparam_start_sigma_dip = nparams+1
                nparams = nparams+1
            end if
            print*,'Do you want to include a separate uncertainty in the restored state dips &
                (which will not be propagated through restoration)? (0/1)'
            do
                if (options_id == 5) then !Manual input
                    read(options_id,*) use_sigma_restored_dip
                else
                    read(options_id,*) str
                    if (str == 'use_sigma_restored_dip') then                    
                        read(options_id,*) use_sigma_restored_dip
                    else
                        backspace(options_id)
                        use_sigma_restored_dip = 0
                    end if                
                end if
                if (use_sigma_restored_dip == 0) then
                    sigma_restored_dip_fixed = 0.;
                    exit
                else if(use_sigma_restored_dip == 1) then
                    if (fit_sigma == 0) then
                        print*,'Enter Uncertainty in Restored State Dip Data: '
                        read(options_id,*) sigma_restored_dip_fixed
                    else
                        nparams = nparams+1
                    end if
                    exit
                else
                    print*,'Invalid Command'
                end if
            end do
        end if
        if (DatTypes(DatType_terr) == 1) then
            allocate(sigma_terr_fixed(nterr))
            allocate(sigma_orig_elev_fixed(nterr))
            if (fit_sigma == 0) then
                print*,'Enter Uncertainties in Terrace Data for all ',nterr,' terraces: '
                read(options_id,*) sigma_terr_fixed
                print*,'Enter Uncertainty in all ',nterr,' Terrace Inner Edge Original Elevations'
                read(options_id,*) sigma_orig_elev_fixed
            else
                nparam_start_sigma_terr = nparams+1
                nparams = nparams+2*nterr
            end if
        end if
        if (DatTypes(DatType_faultpts) == 1) then
            if (fit_sigma == 0) then
                print*,'Enter Uncertainty in Fault Point Data: '
                read(options_id,*) sigma_fault_fixed
            else
                nparam_start_sigma_fault = nparams+1
                nparams = nparams+1
            end if
        end if
        if (DatTypes(DatType_beds_restored) == 1) then
            !Note: Correlation along the length of a bed is not applied to restored bed points.
            !Also, we are assuming that sigma_xy_diff is the same as for beds (when not restored).
            if (sigma_xy_diff == 0) then
                if (fit_sigma == 0) then
                    print*,'Enter Uncertainty in Restored Bed Data: '
                    read(options_id,*) sigma_bed_restored
                    sigma_bed_restored_x_fixed = sigma_bed_restored
                    sigma_bed_restored_y_fixed = sigma_bed_restored
                else
                    nparam_start_sigma_bed_restored = nparams+1
                    nparams = nparams+1
                end if
            else if(sigma_xy_diff == 1) then
                if (fit_sigma == 0) then
                    print*,'Enter Uncertainty in Restored Bed Data (x position): '
                    read(options_id,*) sigma_bed_restored_x_fixed
                    print*,'Enter Uncertainty in Restored Bed Data (y position): '
                    read(options_id,*) sigma_bed_restored_y_fixed
                else
                    nparam_start_sigma_bed_restored = nparams+1
                    nparams = nparams+2
                end if
            else
                print*,'Invalid Command'
            end if
        end if
        exit
    else
        print*,'Invalid Command'
    end if
end do
if (fit_sigma==1) then
    call ReadParams !ReadParams is called before GetResultType. We've added extra parameters, so we need to reread them.
end if
end subroutine GetResultType

subroutine GetErrsFile !Get a name for the errors file and create it
use options, only: errs_id,options_id
implicit none
character (len=50) :: errs_file_name !file name in which to store errors
print*,'Input file name to save results to.'
read(options_id,*) errs_file_name !File name to store the errors in
open(unit=errs_id, file = errs_file_name)
end subroutine GetErrsFile

end module user_input