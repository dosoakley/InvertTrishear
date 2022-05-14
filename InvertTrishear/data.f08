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

module data_module
!The data used in a run, stays the same through all of the models
double precision,dimension(:,:),allocatable :: dip_pos !positions of dips
double precision,dimension(:),allocatable :: dips !values of dips
integer :: ndips !number of dips
type bed
    double precision, dimension(:,:),allocatable :: pts !The points that make up the bed
    integer :: npts !The number of points in the bed
    character (len=20) :: ident !Bed identifier
    logical :: growth !Tells if a bed is in growth strata or not.
end type bed
type(bed),dimension(:),allocatable :: beds !Array of all the beds
integer :: nbeds=0 !number of beds
type(bed),dimension(:),allocatable :: spher_var_R !R matrix for spherical variogram for beds correlated along length.
type(bed),dimension(:),allocatable :: growth_strata !Array of growth strata
integer :: ngrowth=0 !number of growth strata beds.
type(bed),dimension(:),allocatable :: terraces !Array of all the terraces
integer :: nterr=0 !Nuimber of terraces
double precision,dimension(:,:),allocatable :: faultpts !Points on the fault.
integer :: nfaultpts=0 !Number of fault points.
type(bed),dimension(:),allocatable :: restored_beds !Points expected to be on a bed when it is restored.
integer :: n_restored_beds=0 !Number of beds with restored points
end module data_module