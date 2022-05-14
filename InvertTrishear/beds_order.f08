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

module beds_order

integer :: beds_order_type !Tells how the order of beds is determined.
double precision,dimension(2) :: bed_domain_x !x coordinates of the domain in which to test for crossing.

contains

!Ask how to determine the order. This could be included in user_input or can be a separate function to get called.
subroutine GetBedsOrderType
use options, only: options_id
implicit none
character (len = 25) :: str !A string to tell how the order of the beds is determined.
print*,'How do you want to determine order?'
print*,'(1) Use order at x = 0.'
print*,'(2) Prevent beds from crossing only in the domain of x coordinates spanned by both restored beds.'
print*,'(3) Prevent beds from crossing anywhere within a specified domain.'
if (options_id == 5) then !Manual input
    read(options_id,*) beds_order_type
else
    read(options_id,*) str
    if (str=='beds_order_type') then
        read(options_id,*) beds_order_type
    else
        backspace(options_id)
        beds_order_type = 1
    end if
end if
if (beds_order_type == 3) then
    print*,'Enter minimum and maximum x coordinates of the domain to consider.'
    read(options_id,*) bed_domain_x
end if
end subroutine GetBedsOrderType

function check_beds_order(n,b,slope,b_last,slope_last,bed_xlims,bed_xlims_last) result(right_order)
use options, only: beds_age_order
implicit none
integer :: n !The bed number
double precision,intent(in) :: b,b_last !Current and previous bed intercepts.
double precision,intent(in) :: slope,slope_last !Current and previous bed slopes. Note that the sign convention on this is backwards as in the rest of the program.
double precision,dimension(2),intent(in) :: bed_xlims,bed_xlims_last !Current and previous minimum and maximum x coordinates of the restored bed.
double precision,dimension(2) :: xlims !The x limits to consider in determining if the beds are in order.
logical :: right_order
if (beds_order_type==2) then !Get only the x values spanned by both beds.
    xlims(1) = max(bed_xlims(1),bed_xlims_last(1))
    xlims(2) = min(bed_xlims(2),bed_xlims_last(2))
else if (beds_order_type==3) then
    xlims = bed_domain_x
end if
right_order = .true.
if (beds_age_order==1 .and. n>1) then !Youngest to oldest
    if (beds_order_type==1 .and. b>b_last) then !Younger to older means this bed should be older and deeper than the last.
        right_order = .false.
    else if (beds_order_type /= 1) then
        if ((-slope*xlims(1)+b>-slope_last*xlims(1)+b_last) .or. (-slope*xlims(2)+b>-slope_last*xlims(2)+b_last)) then
            right_order = .false.
        end if
    end if
else if(beds_age_order==2 .and. n>1) then !Oldest to youngest
    if (beds_order_type==1 .and. b<b_last) then !Older to younger means this bed should be younger and higher than the last.
        right_order = .false.
    else if (beds_order_type /= 1) then
        if ((-slope*xlims(1)+b<-slope_last*xlims(1)+b_last) .or. (-slope*xlims(2)+b>-slope_last*xlims(2)+b_last)) then
            right_order = .false.
        end if
    end if
end if
end function

end module beds_order