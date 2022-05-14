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

module APT_Module
!Adaptive parallel tempering metropolis algorithm.
!This version uses OpenMP to run the different temperature levels in parallel.

integer :: init_mod_rand !For MCMC, tells if the initial model is (1) random or (2) specified by the user.
double precision, dimension(:),allocatable :: init_model !For MCMC, if initial model is not random, contains the initial parameters
integer :: nlevels !For APT algorithm, number of temperature levels to use.
integer :: nsave !If greater than 1, only every nsave model results will be saved. This is a way to reduce
integer :: nparams_total !The total number of parameters. When APT options is called, some parameters such as restored bed depths may not have been added into the nparams variable yet. This corrects for that omisison.

contains

subroutine APT_options
use parameters, only: nparams,nmodels
use options, only: options_id
implicit none
character (len = 25) :: str
print*,'Enter number of models to run: '
read(options_id,*) nmodels
if (options_id == 5) then !Manual input
    print*,'Enter interval at which to save model results:'
    read(options_id,*) nsave
else
    read(options_id,*) str
    if (str == 'save_at_interval') then                    
        read(options_id,*) nsave
    else
        backspace(options_id)
        nsave=1
    end if                
end if
do
    print*,'Should initial model be (1) random or (2) specified?'
    read(options_id,*) init_mod_rand
    if (init_mod_rand == 1) then
        exit
    else if (init_mod_rand == 2) then
        if (options_id == 5) then !Manual input
            print*,'Enter number of model parameters: (currently expecting at least ',nparams,')'
            read(options_id,*) nparams_total
        else !Reading from a file
            read(options_id,*) str
            if (str == 'nparams_total') then
                read(options_id,*) nparams_total
            else
                backspace(options_id)
                nparams_total = nparams
            end if
        end if
        print*,'Enter initial values for all ',nparams_total,'  parameters: '
        allocate(init_model(nparams_total))
        read(options_id,*) init_model
        exit
    else
        print*,'Invalid Command'
    end if
end do
print*,'Enter number of temperature levels to use: '
do
    read(options_id,*) nlevels
    if (nlevels < 2) then
        print*,'Error: Must have 2 or more levels for APT algorithm.'
    else
        exit
    end if
end do
end subroutine APT_options

subroutine APT !Adaptive Parallel Tempering Algorithm of Miasojedow et al. (2013). Individual chains use the RAM algoirthm (Vihola 2011).
use parameters, only: nparams,mins,maxs,steps,nmodels
use math, only: randnormal,init_random_seed,chol
use options, only: errs_id
use choose_fault, only: choose_fault_model
!$ use omp_lib, only: omp_get_num_threads,omp_get_thread_num
implicit none
double precision,dimension(:,:,:),allocatable :: C,S !Covariance Matrix and Lower triangular decomposition of C. 3rd dimension is temp level.
double precision,dimension(:,:),allocatable :: Ident !Identity Matrix
double precision,dimension(nparams,nlevels) :: params !Current values for all 8 parameters for each level
double precision,dimension(nparams) :: new_params !Proposed values for all the 8 parameters.
integer,dimension(nparams) :: tosolve !Tells which of the variables(model parameters) are actually to be solved for and which are known.
integer :: ntosolve = 0 !Tells how many parameters are to be solved for
double precision,dimension(nparams) :: randnums ! an array of random numbers
double precision :: randn !a single random number
integer(kind = 8) :: model = 0 !model number for counting how many models have been done
double precision,dimension(:,:),allocatable :: x,xnew !Current and proposed model positions. First dimension is params. Second is temp levels.
double precision,dimension(:,:),allocatable :: U,jump !A random vector used in proposing new x, and distance to jump for proposal
double precision,dimension(:,:),allocatable :: xlims !Minimum and maximum limits of the parameters being solved for.
logical :: outofbounds !Specifies if a proposed model is out of the bounds of the model space.
integer :: i,j,k,l ! Counters
double precision,dimension(nlevels) :: lnp !Current model ln of probabilities
double precision :: lnpnew !Proposed model ln of probability
double precision :: pacc,pdes!acceptance probability and desired acceptance probability
double precision :: eta !step size for proposal adaptation (via RAM), decaying to zero
double precision :: gamma !step size for temperature adaptation, decaying to zero
double precision,dimension(nlevels) :: beta !Inverse temperatures.
double precision,dimension(nlevels) :: temp !Temperature = 1/beta
double precision,dimension(nlevels-1) :: rho !Vector of log differences between adjacent temperatures.
double precision,dimension(:),allocatable :: mu_stand,sigma_stand !vector of mu and sigma for standard normal distribution
real,dimension(nlevels) :: acc_count !Counts how many acceptances have been made
real :: acc_count_sw = 0 !Count of how many swap acceptances have been made.
double precision,dimension(nparams) :: params_l !parameters from level l for swapping
double precision,dimension(:),allocatable :: x_l !x values from level l for swapping
double precision :: lnp_l !ln of probability for level l for swapping
integer :: l_swap !The level to swap with the one above it.
!Determine the variables that are actually variable
do i = 1,nparams
    if (maxs(i) /= mins(i)) then
        tosolve(i) = 1
        ntosolve = ntosolve + 1
    else
        tosolve(i) = 0
    end if
end do
allocate(x(ntosolve,nlevels),xnew(ntosolve,nlevels))
allocate(x_l(ntosolve))
allocate(U(ntosolve,1),jump(ntosolve,1))
allocate(C(ntosolve,ntosolve,nlevels),S(ntosolve,ntosolve,nlevels))
allocate(Ident(ntosolve,ntosolve))
allocate(xlims(ntosolve,2))
allocate(mu_stand(ntosolve),sigma_stand(ntosolve))
C = 0 !Initialize all the covariances to 0
Ident = 0
pdes = 0.234 !Desired probability of acceptance
mu_stand = 0
sigma_stand = 1
acc_count = 0 !Initialize the acceptance count to 0
!Get specified initial model or choose a random initial model
call init_random_seed()
!Note: Gfortran documentation says the random number generator is thread safe for OMP, but it does generate all numbers in sequence from a single source, so it may not be ideal.
if (init_mod_rand == 1) then
    i = 0 !Counts how many random initial parameter sets we have tried.
    do
        call random_number(randnums) !If I put this inside the loop, I would get a different random number for each level, which might help the algorithm find good regions quickly.
        do l = 1,nlevels !Is there a way to do this in a single statement without looping through all the levels?
            params(:,l) = mins+randnums*(maxs-mins)
        end do
        call choose_fault_model(lnpnew,params(:,1)) !Just using params for the first level, since at this point they are the same for all levels.
        if (lnpnew /= -huge(0.)) exit
        i = i+1
        if (i == 1e4) then !This allows a maximum of 1e4 sets of parameter values to be tried, but this can be changed.
            print*,'Error: Can not find an allowed initial model.'
            EXIT
        end if
    end do
else
    do l = 1,nlevels
        params(:,l) = init_model
    end do
end if
print*,'Initial Model: ',params(:,1)
call choose_fault_model(lnpnew,params(:,1)) !Just using params for the first level, since at this point they are the same for all levels.
lnp = lnpnew !Since all have the same first model, they all have the same initial probability.
!Make a vector of the x values and the standard deviations for the jumping distances
j=1
do i = 1,nparams
    if (tosolve(i) == 1) then
        x(j,:) = params(i,:)
        C(j,j,:) = steps(i)**2
        Ident(j,j) = 1
        xlims(j,1) = mins(i)
        xlims(j,2) = maxs(i)
        j = j+1
    end if
end do
do l = 1,nlevels
    S(:,:,l) = chol(ntosolve,C(:,:,l)) !Initial S matrix
end do
rho = 1 !Initialize all values of rho to 1.
temp(1) = 1 !Temperature for the base level
do j = 1,nlevels-1  !Calculate temp for the higher levels
    temp(j+1) = temp(j) + exp(rho(j))
end do
beta = 1/temp !Calculate beta = inverse temperature
print*,'Initial beta values: ',beta
!Loop through the full number of models
!open(unit=10, file = 'debugging_output.txt')
do model = 1,nmodels
    eta = min(1.,ntosolve*model**(-2./3.)) !Update eta for use in the RAM updates to the covariance matrices. 
    !$omp parallel default(none) private(l,i,j,k,U,jump,new_params,outofbounds,pacc,lnpnew,randn) &
    !$omp shared(nlevels,nparams,ntosolve,mu_stand,sigma_stand,S,x,xnew,tosolve,params,xlims,beta,lnp,acc_count,eta,model,C, &
    !$omp Ident,pdes)
    !Run the RAM algorithm
    !$omp do 
        do l = 1,nlevels !loop through all the temperature levels.
        U(:,1) = randnormal(ntosolve,mu_stand,sigma_stand) !Random vector from standard normal distribution
        jump = matmul(S(:,:,l),U) !Distance to jump to new model.
        xnew(:,l) = x(:,l)+jump(:,1)!Propose a new model
        j = 1
        do i = 1,nparams !Put in the new values
            if (tosolve(i) == 1) then
                new_params(i) = xnew(j,l)
                j = j+1
            else
                new_params(i) = params(i,l)
            end if
        end do
        !write(10,*),l
        !write(errs_id,*),l,new_params
        outofbounds = .false.
        do k = 1,ntosolve !Check whether any of the parameters values are outside the model space.
            if (xnew(k,l)<xlims(k,1) .or. xnew(k,l)>xlims(k,2)) then
                lnpnew = -huge(0.) !As close to p = 0 and ln(p) = -Inf as I can reasonably make it.
                outofbounds = .true.
                pacc = 0 !No chance of accepting this, but need to know pacc for the adaptation step.
                exit
            end if
        end do
        !write(10,*),outofbounds
        if (outofbounds .eqv. .false.) then
            !write(10,*),'a'
            call choose_fault_model(lnpnew,new_params)
            !write(10,*),'b'
            !Decide whether to accept or not. If accepting, then move to the new point.
            pacc = min(1.,exp(beta(l)*(lnpnew-lnp(l))))
            call random_number(randn)
            if (randn<pacc) then
                acc_count(l) = acc_count(l)+1 !Add to the acceptance count
                x(:,l) = xnew(:,l)
                lnp(l) = lnpnew
                params(:,l) = new_params
            end if
        end if
        !write(10,*),'c'
        !Calculate the new covariance
        C(:,:,l) = matmul(matmul(S(:,:,l),(Ident+eta*(pacc-pdes)*matmul(U,transpose(U))/sum(U**2))),transpose(S(:,:,l)))
        S(:,:,l) = chol(ntosolve,C(:,:,l))
        if (isnan(S(ntosolve,ntosolve,l)) .eqv. .true.) then
           print*,'Error: S is NaN at model ',model
            print*,'level = ',l
        end if
    end do !End the loop through the levels
    !$omp end do
    !$omp end parallel
    !Swap states
    call random_number(randn)
    l_swap = int(randn*(nlevels-1))+1 !Choose a level to propose swapping with the one above it.
    pacc = min(1.,exp((beta(l_swap)-beta(l_swap+1))*(lnp(l_swap+1)-lnp(l_swap))))
    call random_number(randn)
    if (randn<pacc) then !Swap states and associated probabilities.
        x_l = x(:,l_swap)
        params_l = params(:,l_swap)
        lnp_l = lnp(l_swap)
        x(:,l_swap) = x(:,l_swap+1)
        params(:,l_swap) = params(:,l_swap+1)
        lnp(l_swap) = lnp(l_swap+1)
        x(:,l_swap+1) = x_l
        params(:,l_swap+1) = params_l
        lnp(l_swap+1) = lnp_l
        acc_count_sw = acc_count_sw + 1
    end if
    !Adapt the temperature schedule
    gamma = (model+1.)**(-0.6) !The same gamma used in the Miasojedow et al. Matlab script. Note: This is very similar to eta.
    !Note on calculating rho: Miasojedow et al. (2013) have a factor of L-1 (nlevels-1 in my code) before gamma in there in their Matlab code, and I can't see where it came from.
    !   It's not the Pi_rho term in Eqn. 14 of their paper, because that gets multiplied by everything including the old rho. However, from the discussion of gamma in their
    !   paper, it sounds like gamma can have any positive constant in front of it. So maybe (nlevels-1) is their choice for that constant. In my models, however, I find that
    !   either it makes little difference or including the (nlevels-1) factor results in much worse mixing, depending on the model. So I'm leaving it out.
    !rho(l_swap) = rho(l_swap) + (nlevels-1)*gamma*(pacc-pdes)
    rho(l_swap) = rho(l_swap) + gamma*(pacc-pdes)
    do j = 1,nlevels-1  !Calculate temp for the higher levels
        temp(j+1) = temp(j) + exp(rho(j))
    end do
    beta = 1./temp !Calculate beta = inverse temperature
    !Write the results to the file
    if (mod(model,nsave)==0) then
        write(errs_id,*) lnp(1),params(:,1)
    end if
    if (mod(model,10000)==0) then !Print model number every 10000 models
    !if (mod(model,1000)==0) then !Print model number every 10 models
        print*, 'ran model ',model,' of ',nmodels
    end if
end do
!C(:,:,1) = matmul(S(:,:,1),transpose(S(:,:,1)))
print*,'Final Covariance matrix for level 1:'
do i = 1,ntosolve
    print*,real(C(i,:,1))
end do
do l = 1,nlevels
    print*,'Acceptance Ratio for level ',l,' = ',acc_count(l)/nmodels
end do
print*,'Acceptance ratio for swaps = ',acc_count_sw/nmodels
print*,'Final Beta Values: ',beta
end subroutine APT

end module APT_Module