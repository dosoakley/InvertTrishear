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

module MC_Grid_Module

real :: nmodels_real !nmodels as a real to check it's the same as the integer value
integer(kind = 8) :: nmodels_poss !number of possible models (for Monte Carlo only)

contains

subroutine MCGrid_options
use parameters, only: maxs,mins,steps,nmodels
use options, only: options_id
implicit none
nmodels_real = product((maxs-mins)/steps+1)
nmodels_poss = nint(nmodels_real,8)
if (abs(nmodels_poss - nmodels_real) < 1e-5) then !Using < 1e-5 rather than = 0 to avoid small errors due to the computer's imprecision
    print*,'Enter number of models to run: '
    read(options_id,*) nmodels
    print*,nmodels,'models will be run.'
else
    print*,'Error: (max-min) is not evenly divisible by step for all parameters.' !Then what do I do?
end if
end subroutine MCGrid_options

subroutine MonteCarlo_Grid
use parameters, only: nmodels,nparams,mins,maxs,steps
use options, only: errs_id
use math, only: init_random_seed,Makekback,ModelParamsFromNum
use choose_fault, only: choose_fault_model
!$ use omp_lib, only: omp_get_num_threads
implicit none
double precision,dimension(nparams) :: params !Values of all the parameters put into one vector.
integer,dimension(nparams) :: kback !Vector of cumulative product backward of n_vals
real :: randn ! a random number
integer(kind = 8) :: model = 0 !model number for counting how many models have been done
integer(kind = 8) :: modelnum = 0 !model number for determining the model parameters, can be 1 to nmodels_poss
double precision :: RMS !RMS error for a model
integer :: nthreads !The number of parallel threads being run.
!Make the kback vector and initialize the random number generator
call Makekback(kback,nparams,mins,maxs,steps)
!$omp parallel default(none) private(model,params,RMS,randn,modelnum) &
!$omp shared(nmodels,nparams,kback,mins,steps,errs_id,nthreads,nmodels_poss)
call init_random_seed()
!$omp single
!$ nthreads = omp_get_num_threads()
!$ print*,'Starting ',nthreads,' threads'
!$omp end single
!$omp do schedule(guided,100)
!Loop through all the runs
do model = 1,nmodels
    !Randomly choose model parameters
    call random_number(randn)
    modelnum = int(randn*nmodels_poss)+1 !Model number for choosing parameters.
    call ModelParamsFromNum(modelnum,nparams,params,kback,mins,steps) !Determine parameter values
    call choose_fault_model(RMS,params)
    if (mod(model,10000)==0) then !Print model number every 10,000 models
        print*, 'ran model ',model,' of ',nmodels
    end if
    write(errs_id,*) RMS,params!Write the error to the file along with the parameters.
end do
!$omp end do
!$omp end parallel
end subroutine MonteCarlo_Grid

end module MC_Grid_Module