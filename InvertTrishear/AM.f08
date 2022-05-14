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
module AM_Module

integer :: init_mod_rand !For MCMC, tells if the initial model is (1) random or (2) specified by the user.
double precision, dimension(:),allocatable :: init_model !For MCMC, if initial model is not random, contains the initial parameters
integer :: begin_adapt !For AM algorithm, model number at which to begin adapting the covariance matrix

contains

subroutine AM_options
use parameters, only: nparams,nmodels
use options, only: options_id
print*,'Enter number of models to run: '
read(options_id,*) nmodels
do
    print*,'Should initial model be (1) random or (2) specified?'
    read(options_id,*) init_mod_rand
    if (init_mod_rand == 1) then
        exit
    else if (init_mod_rand == 2) then
        print*,'Enter initial values for all ',nparams,'  parameters: '
        allocate(init_model(nparams))
        read(options_id,*) init_model
        exit
    else
        print*,'Invalid Command'
    end if
end do
print*,'Enter model number after which to begin using adapted covariance matrix'
read(options_id,*) begin_adapt
end subroutine AM_options

subroutine AM !Adaptive Metropolis (Haario et al. 2001)
use parameters, only: nparams,mins,maxs,steps,nmodels
use math, only: randnormal_multi,init_random_seed,chol
use options, only: errs_id
use choose_fault, only: choose_fault_model
implicit none
double precision,dimension(:,:),allocatable :: C,C0,Ident !Covariance Matrix, Initial Covariance Matrix and Identity Matrix
double precision,dimension(:,:),allocatable :: L !Cholesky decomposition of C or C0
double precision,dimension(:,:),allocatable :: mean,old_mean !matrix of mean values and their values from the previous step
double precision :: sd,eps !scaling parameter and epsilon for calculating new covariance matrix
double precision,dimension(nparams) :: params,new_params !Current and proposed values for all the 8 parameters.
integer,dimension(nparams) :: tosolve !Tells which of the variables(model parameters) are actually to be solved for and which are known.
integer :: ntosolve = 0 !Tells how many parameters are to be solved for
double precision,dimension(nparams) :: randnums ! an array of random numbers
integer(kind = 8) :: model = 0 !model number for counting how many models have been done
double precision,dimension(:),allocatable :: x,xnew !Current and proposed model positions
double precision,dimension(:,:),allocatable :: x2d !rank 2 version of x
double precision,dimension(:,:),allocatable :: xlims !Minimum and maximum limits of the parameters being solved for.
logical :: outofbounds !Specifies if a proposed model is out of the bounds of the model space.
integer :: i,j,k ! Counters
double precision :: lnp,lnpnew !Current and proposed model probabilities
double precision :: pacc !acceptance probability
double precision :: randn !random number used to decide whether to accept or not
real :: acc_count=0 !Counts how many acceptances have been made
!Determine the variables that are actually variable
do i = 1,nparams
    if (maxs(i) /= mins(i)) then
        tosolve(i) = 1
        ntosolve = ntosolve + 1
    else
        tosolve(i) = 0
    end if
end do
allocate(x(ntosolve))
allocate(xnew(ntosolve))
allocate(C(ntosolve,ntosolve))
allocate(C0(ntosolve,ntosolve))
allocate(L(ntosolve,ntosolve))
allocate(Ident(ntosolve,ntosolve))
allocate(mean(ntosolve,1))
allocate(old_mean(ntosolve,1))
allocate(xlims(ntosolve,2))
allocate(x2d(ntosolve,2))
C = 0. !Initialize all the covariances to 0
C0 = 0.
sd = (2.4**2)/ntosolve !Optimal choice according to Haario et al. (2001) citing work by Gelman et al. (1996)
eps = 1e-9 !A very small number
Ident = 0
mean = 0
!Get specified initial model or choose a random initial model
call init_random_seed()
if (init_mod_rand == 1) then
    call random_number(randnums)
    params = mins+randnums*(maxs-mins)
else
    params = init_model
end if
print*,'Initial Model: ',params
call choose_fault_model(lnp,params)
!Make a vector of the x values and the standard deviations for the jumping distances
j=1
do i = 1,nparams
    if (tosolve(i) == 1) then
        x(j) = params(i)
        C0(j,j) = steps(i)**2
        Ident(j,j) = 1
        xlims(j,1) = mins(i)
        xlims(j,2) = maxs(i)
        j = j+1
    end if
end do
mean(:,1) = x !Initial mean is just the initial parameters.
!Loop through the full number of models
do model = 1,nmodels
    if (model > begin_adapt) then !If sufficient models have been run for it to be time to begin using the adapted covariance matrix
        L = chol(ntosolve,C)
    else
        L = chol(ntosolve,C0)
    end if
    xnew = randnormal_multi(ntosolve,x,L) !Propose a new model
    j = 1
    do i = 1,nparams !Put in the new values
        if (tosolve(i) == 1) then
            new_params(i) = xnew(j)
            j = j+1
        else
            new_params(i) = params(i)
        end if
    end do
    outofbounds = .false.
    do k = 1,ntosolve
        if (xnew(k)<xlims(k,1) .or. xnew(k)>xlims(k,2)) then
            lnpnew = -huge(0.) !As close to p = 0 and ln(p) = -Inf as I can reasonably make it.
            outofbounds = .true.
            exit
        end if
    end do
    if (outofbounds .eqv. .false.) then
        call choose_fault_model(lnpnew,new_params)
        !Decide whether to accept or not. If accepting, then move to the new point.
        pacc = min(1.,exp(lnpnew-lnp))
        call random_number(randn)
        if (randn<pacc) then
            acc_count = acc_count+1 !Add to the acceptance count
            x = xnew
            lnp = lnpnew
            params = new_params
        end if
    end if
    !Calculate the new mean and covariance
    old_mean = mean
    mean(:,1) = ((model-1)*old_mean(:,1)+x)/model
    x2d(:,1) = x !Just params in a rank 2 matrix so that it can be transposed.
    !I'm not sure if this is entirely right. We may need to initialize C in some manner that is not done here.
    C = ((model-1.)/model)*C+(sd/model)*(model*matmul(old_mean,transpose(old_mean))-(model+1)*matmul(mean,transpose(mean)) &
        +matmul(x2d,transpose(x2d))+eps*Ident)
    !Either way, write the results to the file
    write(errs_id,*) lnp,params
    if (mod(model,10000)==0) then !Print model number every 1000 models
        print*, 'ran model ',model,' of ',nmodels
    end if
end do
print*,'Final Covariance matrix:'
do i = 1,ntosolve
    print*,real(C(i,:))
end do
print*,'Acceptance Ratio = ',acc_count/nmodels
end subroutine AM

end module AM_Module