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

module fault_tip
!Contains subroutines for dealing with the fault tip position, such as calculating the initial or final fault tip position from the other 
!and checking whether the fault tip is in an acceptable location.

!Possible changes to make / Features to add:
!1) Do the calculation of tip_seg directly in find_tip_position function, rather than take it in as an argument.
!2) Allow constraints on the unknown position (initial or final) for all fault models.
!3) Use the fault_tip module for all fault models, not just trishear_multi_bend.
!4) Move the TipToSolve option from the options module into this module.
!5) Maybe change the "known" / "unknown" terminology for the tip position not being fit for to something else such as "fit" and "nonfit" or "param" and "nonparam".

integer :: constrain_nonparam_tip !Tells whether or not to constrain the tip position (initial or final) that isn't fit for as a model parameter.
integer :: constraint_type !Tells which type of constraint to use if the non-parameter tip position is being constrained.
double precision :: np_tip_xmin,np_tip_xmax,np_tip_ymin,np_tip_ymax !Limits for the tip position that isn't being fit for as a model parameter to use if constrain_nonparam_tip == 1.
double precision :: constraint_dip,constraint_slope,constraint_intercept !dip, slope, and intercept of a line used as a constraint on the initial fault tip.
integer :: tip_relative_to_line !Tells the position of the tip relative to the line if constraint_type is a line. 1 = above (or left for vertical lines) and 2 = below (or right for vertical lines).

contains

subroutine GetTipOptions !Find out whether to solve for initial or final tip position. (Initial or final in terms of forward motion.)
use options, only: TipToSolve,options_id
use constants, only: pi
implicit none
character (len = 27) :: str !A string to tell whether to constrain the tip position that isn't being fit for.
do
    print*,'Tip position is:'
    print*,'(1) Initial'
    print*,'(2) Final'
    read(options_id,*) TipToSolve
    if (TipToSolve == 1 .or. TipToSolve == 2) then
        exit
    else
        print*,'Invalid Command'
    end if
end do
do
    if (TipToSolve==1) then
        print*,'Do you want to place constraints on the final tip position too? (0/1)'
    else
        print*,'Do you want to place constraints on the initial tip position too? (0/1)'
    end if
    if (options_id == 5) then !Manual input
        read(options_id,*) constrain_nonparam_tip
    else !Reading from a file
        read(options_id,*) str
        if (str == 'constrain_non-parameter_tip') then                    
            constrain_nonparam_tip = 1
        else
            backspace(options_id)
            constrain_nonparam_tip = 0
        end if   
    end if
    if (constrain_nonparam_tip==1) then
        print*,'Warning: Constraining non-parameter tip position will only work if you are using a Straight, Multi-bend, &
            Listric (approximate), or Spline fault type.'
        do
            print*,'Do you want to:'
            print*,'(1) Constrain the tip position within a box?'
            print*,'(2) Constrain the tip position relative to a line?'
            read(options_id,*) constraint_type
            if (constraint_type == 1) then
                print*,'Enter tip position limits as: xmin,xmax,ymin,ymax.'
                read(options_id,*) np_tip_xmin,np_tip_xmax,np_tip_ymin,np_tip_ymax
                exit
            else if (constraint_type == 2) then
                print*,'Enter line as: line dip (in degrees, positive down to left), y-intercept.'
                read(options_id,*) constraint_dip,constraint_intercept
                constraint_slope = tan(constraint_dip*pi/180)
                if (constraint_dip==90) then !Vertical line
                    print*,'Must the tip be (1) left or (2) right of the line?'
                else if (constraint_dip>0) then
                    print*,'Must the tip be (1) above / left of the line or (2) below / right of the line?'
                else if (constraint_dip<0) then
                    print*,'Must the tip be (1) above / right of the line or (2) below / left of the line?'
                else !Horizontal line
                    print*,'Must the tip be (1) above or (2) below the line?'
                end if
                read(options_id,*) tip_relative_to_line
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
end subroutine GetTipOptions

subroutine find_tip_position(tipxknown,tipyknown,tipxinit,tipyinit,tipxfinal,tipyfinal,tip_seg_known,tip_seg_init,tip_seg_final,&
    nsegs,slipseg,bendxy,total_slip,PoverS,R,ramp_angle)
!Calculate either the initial fault tip position from the final position or the final fault tip position from the initial position.
use parameters, only:slip_sense
use options, only: TipToSolve
implicit none
double precision,intent(in) :: tipxknown,tipyknown !The known fault tip position, either initial or final depending on what TipToSolve is.
double precision,intent(out) :: tipxinit,tipyinit,tipxfinal,tipyfinal !Initial and final fault tip positions
integer,intent(in) :: tip_seg_known !Fault segment that the known fault tip position in in.
integer,intent(out) :: tip_seg_init,tip_seg_final !Fault segments that the fault tip is in at its initial and final positions.
integer :: nsegs !Number of fault segments.
double precision,dimension(1:),intent(out) :: slipseg !The amoung of slip with the tip in each segment. (Absolute value.) In the order that the tip moves through them.
double precision,dimension(1:,1:),intent(in) :: bendxy !x,y coordinates of the fault bends. Numbered from highest = 1 to lowest = nbends.
double precision,intent(in) :: total_slip !Total slip on the fault (lowest segment).
double precision,dimension(1:),intent(in) :: PoverS !Propagation to slip ratio.
double precision,dimension(1:),intent(in) :: R !Ratio of slip on a given segment relative to slip on the lowest segment.
double precision,dimension(1:) :: ramp_angle !Fault dips for all segments.
double precision :: slip_acc !The accumulated slip during the process of finding the final tip
double precision :: dist !The distance from the current tip position to the next bend and from (tipx,tipy) to the current bend
integer i
!Determine which segment the tip (tipx,tipy) is in. Note: this can be initial or final tip depending which one tipx,tipy represent.
if (TipToSolve == 1) then !If the initial fault tip position is known.
    !Note: If we allowed tipxfinal to be lower than tipxinit, this wouldn't be so simple. Same for if TipToSolve = 2 below
    tip_seg_init = tip_seg_known
    tipxinit = tipxknown
    tipyinit = tipyknown
    tipxfinal = tipxinit !Initialize these. We will add to them.
    tipyfinal = tipyinit
    slip_acc = 0 !Accumulated slip as we move the tip up the fault toward its final position.
    slipseg = 0 !Slip in each segment
    do i = tip_seg_init,1,-1
        if(i==1) then
            dist = abs(total_slip)*PoverS(i)*R(i)
        else if (i==tip_seg_init) then
            dist = sqrt((tipxinit-bendxy(1,i-1))**2 + (tipyinit-bendxy(2,i-1))**2)
        else !This might need some special case for if i == nsegs. It should be okay, though, because if i==nsegs, then tip_seg_init == nsegs, so the above case will be triggered instead.
            dist = sqrt((bendxy(1,i)-bendxy(1,i-1))**2 + (bendxy(2,i)-bendxy(2,i-1))**2) !Slip in the next segment
        end if
        if(i==1 .or. (abs(total_slip)-abs(slip_acc)<=dist/(R(i)*PoverS(i)))) then !Ends in current segment.
            slipseg(i) = total_slip-slip_acc
            tipxfinal = tipxfinal-slipseg(i)*PoverS(i)*R(i)*cos(ramp_angle(i)) !Note: If normal faults are allowed, there might be a sign issue. Right now, slipseg is negative, so - negative = positive.
            tipyfinal = tipyfinal-slipseg(i)*PoverS(i)*R(i)*sin(ramp_angle(i))
            tip_seg_final = i
            exit
        else
            slipseg(i) = slip_sense*dist/(R(i)*PoverS(i)) !Slip (as measured at the lowest segment) needed to move the tip through this segment.
            slip_acc = slip_acc+slipseg(i)
            tipxfinal = bendxy(1,i-1) !Move the tip positions up to the beginning of the next fault segment.
            tipyfinal = bendxy(2,i-1)
        end if
    end do
else !If the final fault tip position is known.
    tip_seg_final = tip_seg_known
    tipxfinal = tipxknown
    tipyfinal = tipyknown
    tipxinit = tipxfinal !Initialize these. We will move it down from here.
    tipyinit = tipyfinal
    slip_acc = 0 !Accumulated slip as we move the tip down the fault toward its initial position.
    slipseg = 0 !Slip in each segment
    do i = tip_seg_final,nsegs
        if (i==tip_seg_known) then
            dist = sqrt((tipxfinal-bendxy(1,i))**2 + (tipyfinal-bendxy(2,i))**2)
        else if(i/=nsegs) then
            dist = sqrt((bendxy(1,i)-bendxy(1,i-1))**2 + (bendxy(2,i)-bendxy(2,i-1))**2) !Slip in the next segment
        end if
        if(i==nsegs .or. (abs(total_slip)-abs(slip_acc)<=dist/(R(i)*PoverS(i)))) then !Ends in current segment.
            slipseg(i) = total_slip-slip_acc
            tipxinit = tipxinit+slipseg(i)*PoverS(i)*R(i)*cos(ramp_angle(i))
            tipyinit = tipyinit+slipseg(i)*PoverS(i)*R(i)*sin(ramp_angle(i))
            tip_seg_init = i
            exit
        else
            slipseg(i) = slip_sense*dist/(R(i)*PoverS(i)) !Slip (as measured at the lowest segment) needed to move the tip through this segment.
            slip_acc = slip_acc+slipseg(i)
            tipxinit = bendxy(1,i) !Move the tip positions down to the beginning of the next fault segment. Note: If normal faults are allowed, this might be different.
            tipyinit = bendxy(2,i)
        end if
    end do
end if
end subroutine find_tip_position

function check_fault_tip(tipxinit,tipyinit,tipxfinal,tipyfinal)
!If there are restrictions on the non-parameter tip, this checks to see if the tip is within the allowed range for that tip.
!Returns a logical variable (true if the tip is okay, false if it is not).
use options, only: TipToSolve
implicit none
double precision,intent(in) :: tipxinit,tipyinit,tipxfinal,tipyfinal !Initial and final fault tip positions
logical :: check_fault_tip
if (constrain_nonparam_tip==1) then !There is some redundancy in this code that could perhaps be reduced. For instance by making the arguments for this be parameter tip vs. non-parameter tip.
    if (constraint_type==1) then !Box
        if (TipToSolve==1) then !Initial tip is a parameter, final is not.
            check_fault_tip = tipxfinal>=np_tip_xmin .and. tipxfinal<=np_tip_xmax .and. tipyfinal>=np_tip_ymin &
                .and. tipyfinal<=np_tip_ymax
        else if(TipToSolve==2) then !Final tip is a parameter, initial is not.
            check_fault_tip = tipxinit>=np_tip_xmin .and. tipxinit<=np_tip_xmax .and. tipyinit>=np_tip_ymin &
                .and. tipyinit<=np_tip_ymax
        else    !This shouldn't happen
            print*,'Error: TipToSolve value not recognized'
        end if
    else !Line
        if (TipToSolve==1) then !Initial tip is a parameter, final is not.
            if (constraint_dip/=90) then !Not vertical
                if (tip_relative_to_line==1) then
                    check_fault_tip = tipyfinal>constraint_slope*tipxfinal+constraint_intercept
                else
                    check_fault_tip = tipyfinal<constraint_slope*tipxfinal+constraint_intercept
                end if
            else !Vertical line
                if (tip_relative_to_line==1) then !Left
                    check_fault_tip = tipxfinal<(tipyfinal-constraint_intercept)/constraint_slope
                else !Right
                    check_fault_tip = tipxfinal>(tipyfinal-constraint_intercept)/constraint_slope
                end if
            end if
        else if(TipToSolve==2) then !Final tip is a parameter, initial is not.
            if (constraint_dip/=90) then !Not vertical
                if (tip_relative_to_line==1) then
                    check_fault_tip = tipyinit>constraint_slope*tipxinit+constraint_intercept
                else
                    check_fault_tip = tipyinit<constraint_slope*tipxinit+constraint_intercept
                end if
            else !Vertical line
                if (tip_relative_to_line==1) then !Left
                    check_fault_tip = tipxinit<(tipyinit-constraint_intercept)/constraint_slope
                else !Right
                    check_fault_tip = tipxinit>(tipyinit-constraint_intercept)/constraint_slope
                end if
            end if
        else    !This shouldn't happen
            print*,'Error: TipToSolve value not recognized'
        end if
    end if
else
    check_fault_tip = .true.
end if
end function check_fault_tip

end module fault_tip