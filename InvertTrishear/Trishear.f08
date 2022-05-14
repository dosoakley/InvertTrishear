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

program Trishear
!Program for trishear inversion of bed and dip data for fault-propagation folds over thrust faults.
!If you want to invert a normal fault, or run a thrust forward, you should use slip_sense = -1 and search for the initial tip, but tell the program it is the final tip.
!Trishear is calculated using the linear velocity field equations of Zehnder and Allmendinger (2000).

use user_input
use fault_tip, only: GetTipOptions
use choose, only: choose_fault_options,choose_method_options,choose_method
use options, only: DatTypes,errs_id,options_id,DatType_beds,DatType_dips,DatType_faultpts,DatType_terr,DatType_beds_restored
use parameters, only: nmodels,nparams,nparam_start_growth,nparam_start_terr
use data_module, only: ngrowth,nterr

implicit none

integer :: tstart, tend, rate !start and end times, and count rate.

call GetOptionsID
call GetDatType
if (DatTypes(DatType_beds) == 1) then
    call ReadBeds
end if
if (DatTypes(DatType_dips) == 1) then
    call ReadDips
end if
if (DatTypes(DatType_terr) == 1) then
    call ReadTerraces
end if
if (DatTypes(DatType_faultpts) == 1) then
    call ReadFaultPts
end if
if (DatTypes(DatType_beds_restored) == 1) then
    call ReadRestoredBedPts
end if
call GetTipOptions
call GetFaultType
call choose_fault_options
!The fault options should set the initial number of parameters. After that, add more if necessary.
if (DatTypes(DatType_beds) == 1) then !If there are growth strata, the slip to restore each of these is a parameter.
    nparam_start_growth = nparams+1
    nparams = nparams+ngrowth
end if
if (DatTypes(DatType_terr) == 1) then !If there are terraces, the slip to restore each terrace adds a new parameter.
    nparam_start_terr = nparams+1
    nparams = nparams+nterr
end if
call ReadParams
call GetMethod
call choose_method_options
if (DatTypes(DatType_beds)==1 .or. DatTypes(DatType_dips)==1) then
    call GetFitType
end if
if (DatTypes(DatType_terr)==1) then
    call GetTerrFitInfo
end if
call GetResultType
call GetErrsFile
if (options_id /= 5) then
    close(options_id) !Close the file that run options were read from.
end if

call system_clock(tstart,rate) !Start timing for the main part of the program.
call choose_method
call system_clock(tend,rate)
print*,'Ran ',nmodels,' models in ',(tend-tstart)/rate,' seconds.'
close(errs_id) !Close the errors file
print*,'Finished'
read(*,*)
end program Trishear