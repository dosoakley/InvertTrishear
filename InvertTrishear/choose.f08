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

module choose
!Includes subroutines for choosing which algorithms to use for:
!1) Solution algorithm (grid search, monte carlo, etc)
!2) Fault type

contains

subroutine choose_fault_options
use options, only: FaultType
use trishear_straight, only: trishear_str_options
use trishear_decol, only: trishear_decol_options
use trishear_bend, only: trishear_bend_options
use trishear_listric, only: trishear_listric_options
use Parallel_FPF, only:parallel_FPF_options
use trishear_multibend, only: trishear_multibend_options
use trishear_elliptic_approx, only: trishear_elliptic_options
use trishear_spline, only: trishear_spline_options
implicit none
select case(FaultType)
case(1)
    call trishear_str_options
case(2)
    call trishear_decol_options
case(3)
    call trishear_bend_options
case(4)
    call trishear_listric_options
case (5)
    call parallel_FPF_options
case(6)
    call trishear_multibend_options
case(7)
    call trishear_elliptic_options
case(8)
    call trishear_spline_options
case default
    print*,'Error: Unrecognized Fault Type'
end select
end subroutine choose_fault_options

subroutine choose_method_options !Choose the options subroutine related to the chosen method
use options, only: method
use GridSearch_Module, only: GridSearch_options
use MC_Grid_Module, only: MCGrid_options
use MCMC_module, only: MCMC_options
use AM_Module, only: AM_options
use RAM_Module, only: RAM_options
use APT_Module, only: APT_options
use MC_normal_module, only: MCnormal_options
!use RML_module, only: RML_options
implicit none
select case(method)
case(1)
    call GridSearch_options
case(2)
    call MCGrid_options
case(3)
    call MCnormal_options
case(4)
    call MCMC_options
case(5)
    call AM_options
case(6)
    call RAM_options
case(7)
    call APT_options
!case(8)
!    call RML_options
case default
    print*,'Error: Unrecognized Method'
end select
end subroutine choose_method_options

subroutine choose_method !Choose the solution method
use options, only: method
use GridSearch_Module,only: GridSearch
use MC_Grid_Module, only: MonteCarlo_Grid
use MCMC_module, only: MCMC
use AM_Module, only: AM
use RAM_Module, only: RAM
use APT_Module, only: APT
use MC_normal_module, only: MC_normal
!use RML_module, only: RML
implicit none
select case(method)
case(1)
    call GridSearch
case(2)
    call MonteCarlo_Grid
case(3)
    call MC_normal
case(4)
    call MCMC
case(5)
    call AM
case(6)
    call RAM
case(7)
    call APT
!case(8)
!    call RML
end select
end subroutine choose_method

end module choose