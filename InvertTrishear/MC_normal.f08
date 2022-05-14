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

module MC_normal_module
!Module for Monte Carlo simulations that draw samples from a multivariate normal distribution.

double precision, dimension(:),allocatable :: means !means of the multivariate normal distribution
double precision,dimension(:,:),allocatable :: C !Covariance of the multivariate normal distribution
integer :: dist_id = 8 !id for the file that holds the parameters for the normal distribution

contains

subroutine MCnormal_options
use parameters, only: nmodels,nparams
use options, only: options_id
implicit none
character (len=50) :: file_name
integer :: i !a counter
print*,'Enter number of models to run: '
read(options_id,*) nmodels
print*,nmodels,'models will be run.'
print*,'Enter filename to read proposal distribution from: '
read(options_id,*) file_name
allocate(means(nparams))
allocate(C(nparams,nparams))
open(unit=dist_id,file=file_name)
read(dist_id,*) means
read(dist_id,*) !Skip a line.
do i = 1,nparams
    read(dist_id,*) C(i,:) !Read each line of the covariance matrix
end do
close(dist_id)
end subroutine MCnormal_options

subroutine MC_normal
use parameters, only: nmodels,nparams
use options, only: errs_id
use math, only: init_random_seed,randnormal_multi,chol
use choose_fault, only: choose_fault_model
!$ use omp_lib, only: omp_get_num_threads
implicit none
double precision,dimension(nparams,nparams) :: L !Cholesky decomposition of C
double precision,dimension(nparams) :: params !Values of all the parameters put into one vector.
integer(kind = 8) :: model = 0 !model number for counting how many models have been done
double precision :: error !error for a model
integer :: nthreads !The number of parallel threads being run.
L = chol(nparams,C)
!$omp parallel default(none) private(model,params,error) &
!$omp shared(nmodels,nparams,means,C,L,errs_id,nthreads)
call init_random_seed() !Initialize the seed for the pseudorandom number generator.
!$omp single
    !$ nthreads = omp_get_num_threads()
    !$ print*,'Starting ',nthreads,' threads'
!$omp end single
!$omp do schedule(guided,1000)
do model = 1,nmodels
    params = randnormal_multi(nparams,means,L) !Note: This does Cholesky decomposition of C every time. It would be better to do just once.
    call choose_fault_model(error,params)
    if (mod(model,10000)==0) then !Print model number every 10,000 models
        print*, 'ran model ',model,' of ',nmodels
    end if
    write(errs_id,*) error,params!Write the error to the file along with the parameters.
end do
!$omp end do
!$omp end parallel
end subroutine MC_normal

end module MC_normal_module