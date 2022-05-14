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

module options
!Contains variables specifying user options.
integer :: options_id !Identifier of the source of run options, 5 for keyboard or 4 for a file
integer :: params_id = 2 !Identifier of the file used to read parameter values.
integer :: errs_id = 3 !Identifier of the file used to store error values
integer :: method !Tells what method of solution to use (grid search or Monte Carlo)
integer :: ndatatypes !Number of different data types to be used.
integer :: ndatatypes_poss = 5 !Total possible number of data types
integer,dimension(5) :: DatTypes !Tells whether each type of data is present (0/1). 1st spot =beds, 2nd=dips, 3rd=terraces, 4th = points on the fault
integer :: DatType_beds=1,DatType_dips=2,DatType_terr=3,DatType_faultpts=4,DatType_beds_restored = 5 !Tells which number (and position in DatTypes) corresponds to which data type.
integer :: FaultType !Tells the type of fault. 1=straight, 2=ramp from detachment
integer :: FitType !Type of fit to perform to data
integer :: ResultType !Tells what type of result the run should produce
integer :: TipToSolve !Tip to solve for. 1 = initial, 2 = final
integer :: propagate = 0 !Tells whether or not to propagate errors.
integer :: sigma_xy_diff !Tells whether uncertainty in bed positions is different in x and y.
integer :: diff_sigma_bed !Tells whether there is a different uncertainty for each bed.
integer :: n_sigma_bed_groups !If bed uncertainties are grouped, number of groups there are.
integer,dimension(:),allocatable :: sigma_bed_group_start !Tells on which bed number each bed uncertainty group starts.
integer :: use_sigma_restored_dip !Tells whether or not to use sigma_restored_dips (0/1).
integer :: fit_sigma !Tells whether or not we should fit for sigma as a model parameter.
integer :: state_for_Lc !Tells whether distance for beds correlated along length should be calculated in the deformed state (1) or the restored state (2).
integer :: fit_Lc !Tells whether or not to fit for correlation length as a model parameter. (0 = don't fit, 1 = do fit).
integer :: fit_dips !Tells whether we should not fit for beds as a parameter (0), fit for a single restored dip for all data (1), fit for each bed's restored dip individually (2), fit for groups of beds with the same dip (3).
integer :: n_bed_groups !Tells the number of groups of beds with the same dip, if fit_dips==3.
integer,dimension(:),allocatable :: bed_group_start !Tells on which bed number each group starts.
integer :: fit_beds !Tells whether we should fit for bed intercepts as a parameter (0=no, 1=yes).
integer :: n_bed_segs !Number of segments per bed if fitting to multi-segment restored beds.
integer :: terr_age_order !Tells if terraces are in order by age: 0 = no order; 1 = youngest to oldest; 2 = oldest to youngest
integer :: beds_age_order !Tells if beds are in order by age: 0 = no order; 1 = youngest to oldest; 2 = oldest to youngest
integer :: limit_bed_dips=0 !If FitType is 3, this tells whether or not to put min and max limits on the bed dips (0 or 1).
double precision,dimension(:),allocatable :: min_bed_slopes,max_bed_slopes !If limit_bed_dips is 1, these give the minimum and maximum allowed slope for each bed.
end module options