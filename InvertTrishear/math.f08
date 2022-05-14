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

module math

contains

function xy_to_ze2D(pts,tipx,tipy,ramp_acute,ramp_dir) !Convert from regular x,y coordinates to trishear zeta, eta coordinates
implicit none
double precision, dimension(1:,1:) :: pts
double precision, dimension(:,:),allocatable :: ptsze,ptszet
double precision :: tipx, tipy, ramp_cos, ramp_sin
double precision :: ramp_acute !acute ramp angle
double precision, dimension(size(pts,1),size(pts,2)) :: xy_to_ze2D
integer :: i
integer :: ramp_dir !+1 for a right verging fault, -1 for a left verging fault
allocate(ptsze(size(pts,1),size(pts,2)),ptszet(size(pts,1),size(pts,2)))
ramp_cos = cos(-ramp_acute) !It's quicker to just calculate these once
ramp_sin = sin(-ramp_acute)
do i = 1,size(pts,2)
!Translate
ptszet(1,i) = pts(1,i) - tipx
ptszet(2,i) = pts(2,i) - tipy
!Reflect
if (ramp_dir == -1) then
    ptszet(1,i) = -ptszet(1,i);
end if
!Rotate
ptsze(1,i) = ramp_cos*ptszet(1,i)-ramp_sin*ptszet(2,i)
ptsze(2,i) = ramp_sin*ptszet(1,i)+ramp_cos*ptszet(2,i)
end do
xy_to_ze2D = ptsze
end function

function ze_to_xy2D(pts,tipx,tipy,ramp_acute,ramp_dir) !Convert from trishear zeta, eta coordinates to regular x,y coordinates
double precision, dimension(1:,1:) :: pts
double precision, dimension(:,:),allocatable :: ptsxy,ptsxyr
double precision :: tipx, tipy, ramp_cos, ramp_sin
double precision :: ramp_acute !acute ramp angle
double precision, dimension(size(pts,1),size(pts,2)) :: ze_to_xy2D
integer :: i
integer :: ramp_dir !+1 for a right verging fault, -1 for a left verging fault
allocate(ptsxy(size(pts,1),size(pts,2)),ptsxyr(size(pts,1),size(pts,2)))
ramp_cos = cos(ramp_acute) !It's quicker to just calculate these once
ramp_sin = sin(ramp_acute)
do i = 1,size(pts,2)
!Rotate
ptsxyr(1,i) = ramp_cos*pts(1,i)-ramp_sin*pts(2,i)
ptsxyr(2,i) = ramp_sin*pts(1,i)+ramp_cos*pts(2,i)
!Reflect
if (ramp_dir == -1) then
    ptsxyr(1,i) = -ptsxyr(1,i);
end if
!Translate
ptsxy(1,i) = ptsxyr(1,i) + tipx
ptsxy(2,i) = ptsxyr(2,i) + tipy
end do
ze_to_xy2D = ptsxy
end function

function rot2D(pts,angle) !Rotates n points in a 2xn array through a given angle
double precision, dimension(1:,1:) :: pts
double precision :: angle, angle_cos, angle_sin
double precision, dimension(:,:),allocatable :: pts_rot !The rotated points
double precision,dimension(size(pts,1),size(pts,2)) :: rot2D
integer :: i
allocate(pts_rot(size(pts,1),size(pts,2)))
angle_cos = cos(angle)
angle_sin = sin(angle)
do i = 1,size(pts,2)
    pts_rot(1,i) = angle_cos*pts(1,i)-angle_sin*pts(2,i)
    pts_rot(2,i) = angle_sin*pts(1,i)+angle_cos*pts(2,i)
end do
rot2D = pts_rot
end function

subroutine ModelParamsFromNum(model,nparams,params,kback,mins,steps)
!Determines model parameters based on the model number
implicit none
integer(kind=8),intent(in) :: model !model number
double precision,dimension(nparams),intent(in) :: mins,steps !minimum, maximum, and step values for parameters.
integer,intent(in) :: nparams !number of model parameters
integer,dimension(nparams),intent(in) :: kback !Vector of cumulative product backward of n_vals
double precision,dimension(nparams),intent(out) :: params !model parameters
double precision, dimension(nparams) :: subsc !vector of subscript values
integer(kind=8) :: ndx !keeps track of remaining number
integer :: j ! a counter
integer(kind=8) :: vi,vj !internal variables in calculating subscripts from indices
!Convert the model number to a subscript
ndx = model !keeps track of remaining number
do j = 1,nparams
    vi = mod(ndx-1, kback(j)) + 1
    vj = (ndx - vi)/kback(j) !No +1 here because we want to start at 0 not 1
    subsc(j) = vj
    params(j) = mins(j)+steps(j)*subsc(j) !Calculate each of the parameter values
    ndx = vi
end do
end subroutine ModelParamsFromNum

subroutine Makekback(kback,nparams,mins,maxs,steps)
implicit none
integer,intent(in) :: nparams
double precision,dimension(nparams),intent(in) :: mins,maxs,steps !minimum, maximum, and step values for parameters.
integer,dimension(nparams),intent(out) :: kback !Vector of cumulative product backward of n_vals
integer,dimension(nparams) :: n_vals !Vector of all the numbers of values for the different parameters
integer :: i !a counter
!Determine the number of values tried for each parameter
do i = 1,nparams
    n_vals(i) = nint((maxs(i)-mins(i))/steps(i))+1
end do
!Make the kback vector
do i = 1,nparams-1
    kback(i) = product(n_vals(i+1:nparams))
end do
kback(nparams) = 1
end subroutine Makekback

function randnormal(n,mu,sigma) !Random uncorrelated numbers from a normal distribution, using Box-Muller transform
use constants
implicit none
integer,intent(in) :: n !number of samples to take
double precision,intent(in),dimension(n) :: mu,sigma !means and standard deviations of the distribution
double precision,dimension(2) :: rn !Random numbers
double precision :: r,theta !r and theta for Box-Muller transform
integer :: i !a counter for how many samples have been taken
double precision,dimension(n) :: randnormal !Output vector
i = 0
do
    if (i < n) then
        i = i+1
        call random_number(rn)
        r = (-2.*log(rn(1)))**0.5
        theta = 2.*pi*rn(2)
        randnormal(i) = r*cos(theta)*sigma(i)+mu(i)
        if (i<n) then
            i = i+1
            randnormal(i) = r*sin(theta)*sigma(i)+mu(i)
        end if
    else
        exit
    end if
end do
end function

function randnormal_multi(n,mu,L) !Random numbers from a multivariate normal distribution. Can only take one sample at a time.
use constants
implicit none
integer,intent(in) :: n !number of dimensions
double precision,intent(in),dimension(n) :: mu !vector of means for each dimension
double precision,intent(in),dimension(n,n) :: L !Cholesky decomposition of Covariance matrix
double precision,dimension(n) :: z !vector of normally distributed random numbers.
double precision,dimension(n) :: mu_stand,sigma_stand !vector of mu and sigma for standard normal distribution
double precision,dimension(n) :: randnormal_multi !Output vector
!Get random numbers from a standard normal distribution
mu_stand = 0
sigma_stand = 1
z = randnormal(n,mu_stand,sigma_stand)
!Convert to random numbers from the multivariate distribution desired
randnormal_multi = mu + matmul(L,z)
end function randnormal_multi

function chol(n,A) !Takes the lower triangular cholesky decomposition of a positive definite square matrix, using formulas as given in Kreyszig (2006).
integer,intent(in) :: n !matrix dimension
double precision,intent(in),dimension(n,n) :: A !Matrix to be decomposed
double precision,dimension(n,n) :: L !Lower triangular Cholesky decomposition of A
integer :: p,j,s !counters
double precision :: sum1,sum2 !Keep track of sums used in calculating terms for L
double precision,dimension(n,n) :: chol !The result of the function. Equals L.
L = 0 !Initialize all to zero
L(1,1) = sqrt(A(1,1)) !Calculate first element
do j = 2,n !Do the j,1 terms
        L(j,1) = A(j,1)/L(1,1)
end do      
do j = 2,n
    sum1 = 0
    do s = 1,j-1
        sum1 = sum1+L(j,s)**2
    end do
    L(j,j) = sqrt(A(j,j) -sum1)
    if (A(j,j)-sum1 <= 0) then
        print*,'Error: Not positive definite'
        print*,j
        print*,A(j,j)
        print*,sum1
    end if
    do p = j+1,n
        sum2 = 0
        do s = 1,j-1
            sum2 = sum2+L(j,s)*L(p,s)
        end do
        L(p,j) = (1/L(j,j))*(A(p,j)-sum2)
    end do
end do
chol = L
end function chol

function matinv(n,C) !Calculates the inverse (C^-1) of a matrix (C)
implicit none
integer, intent(in) :: n !matrix dimension
double precision,dimension(n,n),intent(in) :: C !The matrix to be inverted
double precision,dimension(n,n) :: Cinv !The inverse of C
double precision,dimension(n,n) :: L,U,Y !Lower and upper triangular matrices, and Y = U*cinv
integer :: i,j,k !Counters
double precision,dimension(n,n) :: matinv !The result of the function. Equals Cinv

!Perform LU decomposition
L = 0 !Initialize L and U to zero
U = 0
L(1,1) = 1
U(1,1) = C(1,1)
do j = 2,n
    U(1,j) = C(1,j)/L(1,1)
    L(j,1) = C(j,1)/U(1,1)
end do
do i = 2,n-1
    L(i,i) = 1
    U(i,i) = C(i,i)
    do k = 1,i-1 !This is really just a matrix multiplication, but I can't seem to do it with matmul without getting an error.
        U(i,i) = U(i,i) - L(i,k)*U(k,i)
    end do !Note if U(i,i) = 0 at the end of this, there is a problem.
    if (U(i,i) == 0) then
        print*,'Error: Matrix cannot be inverted.'
    end if
    do j = i+1,n
        U(i,j) = C(i,j)
        L(j,i) = C(j,i)
        do k = 1,i-1
            U(i,j) = U(i,j) - L(i,k)*U(k,j)
            L(j,i) = L(j,i) - L(j,k)*U(k,i)
        end do
        U(i,j) = U(i,j)/L(i,i)
        L(j,i) = L(j,i)/U(i,i)
    end do
end do
L(n,n) = 1
U(n,n) = C(n,n)
do k = 1,n-1
    U(n,n) = U(n,n) - L(n,k)*U(k,n)
end do
!Solve for Cinv from the LU decomposition, since L*U*Cinv = I
Y(1,1) = 1 !L(1,1)*Y(1,1) = 1, and L(1,1) = 1, so Y(1,1) = 1
Y(1,2:n) = 0 !(L(1,1)*Y(1,i) = 0 for i /= 1
do i = 2,n
    do j = 1,n
        if (i == j) then !Initialize to the value of I(i,j)
            Y(i,j) = 1
        else
            Y(i,j) = 0 !Note: I think the final result will be Y(i,j) = 0 for j>i as is the case for i=1. Look into this if possible.
        end if
        do k = 1,i-1
            Y(i,j) = Y(i,j) - L(i,k)*Y(k,j) !Subtract all the terms using Y(k,j) for k < i
        end do !I would then need to divide Y(i,j) by L(i,i), but L(i,i) = 1 is already known
    end do
end do
!Then solve U*Cinv = y for Cinv
Cinv(N,N) = Y(n,n)/U(n,n)
do i = n,1,-1
    do j = n,1,-1
        Cinv(i,j) = Y(i,j) !Initialize to the value of Y(i,j)
        do k = i+1,n
            Cinv(i,j) = Cinv(i,j) - U(i,k)*Cinv(k,j) !Subtract all the terms using Cinv(k,j) for i < k < bedlength
        end do
        Cinv(i,j) = Cinv(i,j)/U(i,i)
    end do
end do
matinv = Cinv
end function matinv

subroutine swap(a,b)
!Swaps two 1D arrays: a and b
implicit none
double precision,dimension(1:),intent(inout) :: a,b !The two arrays to swap.
double precision,dimension(:),allocatable :: temp !Temporary array
allocate(temp(size(a)))
temp = a
a = b
b = temp
end subroutine swap

function maxind(a)
!Finds the index of the maximum value of array a
implicit none
double precision,dimension(1:),intent(in) :: a
integer :: maxind
integer :: i
maxind = 1
do i = 1,size(a,1)
    if (a(i) > a(maxind)) then
        maxind = i
    end if
end do
end function maxind

subroutine init_random_seed()
!Initializes random number seeds.
implicit none
integer,dimension(:),allocatable :: seed
integer :: n,time,i
call random_seed(size = n)
allocate(seed(n))
call system_clock(count = time)
do i = 1,n
    seed(i) = time+37*(i-1)
end do
call random_seed(put=seed)
end subroutine init_random_seed

function line_point_dist(x,y,m,b)
    !Finds the distance from a point to a line
    implicit none
    double precision,intent(in) :: x,y,m,b !x,y = coordinates of point, m = slope of line, b = intercept of line
    double precision :: line_point_dist !Distance between the point and the line
    line_point_dist = abs(y-m*x-b)/sqrt(m**2+1.)
end function line_point_dist

function sec(x)
    !Finds the secant
    implicit none
    double precision,intent(in) :: x
    double precision :: sec
    sec = 1./cos(x)
end function sec

function csc(x)
    !Finds the cosecant
    implicit none
    double precision,intent(in) :: x
    double precision :: csc
    csc = 1./sin(x)
end function csc

function cot(x)
    !Finds the cotangent
    implicit none
    double precision,intent(in) :: x
    double precision :: cot
    cot = 1./tan(x)
end function cot

end module math