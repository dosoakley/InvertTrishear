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

module GridSearch_Module

double precision :: nmodels_real !nmodels as a real to check it's the same as the integer value

contains

subroutine GridSearch_options
use parameters, only: maxs,mins,steps,nmodels,nparams
implicit none
double precision,dimension(nparams) :: a
a = (maxs-mins)/steps+1
!nmodels_real = product((maxs-mins)/steps+1) !I used to do it this way, which is simpler, as it's all one step, but I started getting an error.
nmodels_real = product(a)
nmodels = nint(nmodels_real,8)
if (abs(nmodels - nmodels_real) < 1e-5) then !Using < 1e-5 rather than = 0 to avoid small errors due to the computer's imprecision
    print*,nmodels,'models will be run.'
else
    print*,'Error: (max-min) is not evenly divisible by step for all parameters.' !Then what do I do?
end if
end subroutine GridSearch_options

subroutine GridSearch
use parameters, only: nparams,mins,maxs,steps,nmodels
use options, only: errs_id
use choose_fault, only: choose_fault_model
use math, only: Makekback,ModelParamsFromNum
!$ use omp_lib, only: omp_get_num_threads,omp_get_thread_num
implicit none
double precision,dimension(nparams) :: params !Values of all the parameters put into one vector.
integer,dimension(nparams) :: kback !Vector of cumulative product backward of n_vals
integer(kind=8) :: model !model number
integer(kind=8) :: model_last = 0 !The last model run in the previous loop through buffer_size models
double precision :: RMS !RMS error for a model
integer :: nwrite !Number of models to write from this block.
integer :: buffer_size !Size of the write buffer. Once this is filled. It gets written to file.
double precision,dimension(:),allocatable :: write_buffer !Holds values to be written to file
integer :: nthreads !The number of parallel threads being run.
integer :: nblocks,block !nblocks = Number of blocks in which the results are calculated and then written all together.
!Make kback vector for use with ModelParamsFromNum
call Makekback(kback,nparams,mins,maxs,steps)
!Just iterate through all the models
buffer_size = 100000 !Seems like a good number.
allocate(write_buffer(buffer_size))
nblocks = ceiling(nmodels_real/buffer_size)
do block = 1,nblocks
    !$omp parallel default(none) private(model,params,RMS) &
    !$omp shared(nmodels,nparams,kback,mins,steps,errs_id,model_last,buffer_size,write_buffer,nthreads,block)
    !$omp single
    if (block == 1) then
        !$ nthreads = omp_get_num_threads()
        !$ print*,'Starting ',nthreads,' threads'
    end if
    !$omp end single
    !$omp do schedule(guided,100)
    do model = model_last+1,min(model_last+buffer_size,nmodels)
        call ModelParamsFromNum(model,nparams,params,kback,mins,steps)
        call choose_fault_model(RMS,params)
        write_buffer(model-model_last) = RMS
    end do
    !$omp end do
    !$omp end parallel
    nwrite = min(buffer_size,nmodels-model_last)
    model_last = min(model_last+buffer_size,nmodels)
    print*,'ran model',model_last,' of ',nmodels
    do model = 1,nwrite
        write(errs_id,*) write_buffer(model)!Write the error to the file.
    end do
end do
end subroutine GridSearch

end module GridSearch_Module