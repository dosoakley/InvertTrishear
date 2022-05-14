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

module choose_fault

contains

subroutine choose_fault_model(output,params)
use options, only: FaultType
use parameters, only: nparams
use trishear_straight, only: run_model_str
use trishear_decol, only: run_model_decol
use trishear_bend, only: run_model_bend
use trishear_listric, only: run_model_listric
use Parallel_FPF, only: run_model_pFPF
use trishear_multibend, only: run_model_multibend
use trishear_elliptic_approx, only: run_model_elliptic
use trishear_spline, only: run_model_spline
implicit none
double precision, dimension(nparams),intent(inout) :: params
double precision, intent(out) :: output
select case(FaultType)
case(1)
    call run_model_str(output,params)
case(2)
    call run_model_decol(output,params)
case(3)
    call run_model_bend(output,params)
case(4)
    call run_model_listric(output,params)
case(5)
    call run_model_pFPF(output,params)
case(6)
    call run_model_multibend(output,params)
case(7)
    call run_model_elliptic(output,params)
case(8)
    call run_model_spline(output,params)
end select
end subroutine choose_fault_model

end module