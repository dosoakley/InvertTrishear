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

module parameters
!Parameters of the run as a whole, stay the same for each model.
double precision,dimension(:),allocatable :: mins,maxs,steps !minimum, maximum, and step values for parameters.
integer :: nparams !number of model parameters
integer :: slip_sense !tells whether slip should be thrust (1) or reverse (-1)
double precision :: increment !increment of slip
double precision :: v0 !slip velocity along fault
integer(kind = 8) :: nmodels !number of models to run
!Parameters telling the expected dip or intercepts of  data:
double precision :: RegDip,RegDipDeg,RegSlope,RegDipSin,RegDipCos !For regional dip. Note: These will probably be removed and transferred to a different module if I separate out different fit methods.
double precision,dimension(:),allocatable :: intercepts !y intercepts of known lines to fit beds to. Note: As for regional dip, will probably be moved soon.
!Parameters relating to growth strata
integer :: nparam_start_growth !Parameter number of the first parameter with growth_slip values.
!Parameters relating to the restored geometry of marine terraces:
integer :: nparam_start_terr !Parameter number of the first parameter with terr_slip values.
double precision,dimension(:),allocatable :: terr_dip,terr_dip_deg,terr_slope !Dip at which terrace surfaces are formed.
double precision,dimension(:),allocatable :: terr_orig_elev !Original terrace inner edge elevations
!Parameters to be used if FitType == 5 (fit for restored dips and/or bed intercepts as parameters):
integer :: nparam_start_restored_fit !Parameter number of the first parameter for the restored lines to fit to.
integer :: nrestored_fit_params=0 !Number of restored geometry parameters to be fit.
end module parameters