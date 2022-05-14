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

module growth_strata
!These are some functions for interpreting slip to restore growth strata, so I don't have to continually repeat this code.
!In future, I want to create a beds module and incorporate this into it.

integer :: growth_slip_parameterization !1: Total slip for each bed. 2: Additional slip since previous bed.

contains

subroutine GetGrowthParameterization
!Find out how growth strata slip should be parameterized.
use options, only: options_id
implicit none
character (len = 30) :: str
do
    print*,'How is slip on growth strata parameterized?'
    print*,'    (1) Slip is total for each growth layer.'
    print*,'    (2) Slip is additional since previous layer.'
    print*,'        (In order in data file from top to bottom. Put youngest on top if using this option).'
    if (options_id == 5) then !Manual input
        read(options_id,*) growth_slip_parameterization
    else !Reading from a file.
        read(options_id,*) str
        if (str == 'growth_slip_parameterization') then
            read(options_id,*) growth_slip_parameterization
        else
            growth_slip_parameterization = 1
            backspace(options_id)
        end if
    end if
    if (growth_slip_parameterization==1 .or. growth_slip_parameterization==2) then
        exit
    else
        print*,'Error: Unrecognized growth strata slip parameterization.'
    end if
end do
end subroutine GetGrowthParameterization

subroutine GetGrowthSlip(growth_slip_params,total_slip_mag,slip_sense,growth_slip,goodrun)
    use parameters, only: nparam_start_growth
    use data_module, only: ngrowth
    use options, only: beds_age_order
    implicit none
    double precision, dimension(:), intent(in) :: growth_slip_params !The parameters that refer to slip needed to restore growth strata.
    double precision, dimension(:), intent(out) :: growth_slip !Slip associated with each growth bed.
    double precision, intent(in) :: total_slip_mag !Magnitude (absolute value) of total slip on the fault.
    integer, intent(in) :: slip_sense !Sense of slip on the fault.
    logical, intent(out) :: goodrun !Tells if these parameters are acceptable.
    integer :: i !A counter
    if (growth_slip_parameterization==1) then
        growth_slip = growth_slip_params
    else
        growth_slip(1) = growth_slip_params(1)
        do i=2,ngrowth
            growth_slip(i) = growth_slip(i-1)+growth_slip_params(i)
        end do
    end if
    if (maxval(growth_slip) > total_slip_mag) then !Note: terr_slip at this point should be positive.
        goodrun = .false.
    end if
    if (beds_age_order==1 .and. ngrowth>1) then !Youngest to oldest
        if (minval(growth_slip(2:ngrowth)-growth_slip(1:(ngrowth-1)))<0) then
            goodrun = .false.
        end if
    else if(beds_age_order==2 .and. ngrowth>1) then !Oldest to youngest
        if (minval(growth_slip(1:(ngrowth-1))-growth_slip(2:ngrowth))<0) then
            goodrun = .false.
        end if
    end if
    growth_slip = growth_slip*slip_sense
end subroutine GetGrowthSlip

end module growth_strata