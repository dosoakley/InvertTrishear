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

module err_and_prob
!Module of subprograms for calculating the error or probability of a model.

contains

subroutine calc_bed_err(n,bedlength,thisbed,beds_output,slope,b,Cbed,Lc,restored_fit_params,Cbed2) !Calculate the error for an individual bed.
use parameters, only: RegSlope,intercepts,nparam_start_restored_fit,nrestored_fit_params
use options, only: FitType,ResultType,propagate,fit_dips,fit_beds,DatTypes,DatType_dips,n_bed_segs,n_bed_groups,bed_group_start,&
    state_for_Lc,fit_Lc
use data_module, only: nbeds,spher_var_R
use data_uncertainties, only: GetSigmaBed
use math, only: matinv,chol
use constants, only: pi
implicit none
integer, intent(in) :: n !The number of the bed being worked on.
integer, intent(in) :: bedlength !number of points in the bed being worked on
double precision, dimension(2,bedlength), intent(in) :: thisbed
double precision,dimension(nbeds),intent(out) :: beds_output !output for each bed
double precision,intent(out) :: slope,b !Slope and intercept of the bed
double precision,dimension(2,2,bedlength) :: Cbed!Use if we are propagating errors. Covariance matrices for x and y for each point.
double precision :: Lc !Correlation length; only used if ResultType==4 or 5.
double precision,dimension(1:),optional :: restored_fit_params !Use if restored dips or intercepts are being fit as a parameter.
double precision,dimension(2,2,bedlength),optional :: Cbed2 !For second (correlated) errors there are separate correlated and uncorrealted errors.
double precision,dimension(bedlength) :: R !residuals between model and data
double precision,dimension(bedlength) :: sigma,sigma2 !Array of sigmas for all values.
double precision :: Delta!Used in calculating best fit line to a restored bed
double precision,dimension(bedlength,bedlength) :: C !Covariance matrix for R. (for correlated probability along bed)
double precision,dimension(bedlength,bedlength) :: S !Matrix with diagonal elements equal to sigma and other elements equal to 0.
double precision,dimension(2,n_bed_segs-1) :: bendxy !x,y positions of bends in the restored bed.
double precision,dimension(n_bed_segs) :: bed_dips !Dips of the segments in a multi-segment restored bed.
integer :: group !Tells which group of beds this bed is in, if fitting dips for groups of beds with the same dip.
integer :: i,j !counters
!Calculate residuals for this bed
select case(FitType)
case (1) !Flat Line
    b = sum(thisbed(2,:))/bedlength !Just average y value
    slope = 0
case (2) !Known dip
    b = sum(thisbed(2,:)+RegSlope*thisbed(1,:))/bedlength !Calculate y intercept of best fit line with RegSlope as its slope
    slope = RegSlope
case (3) !Fit for dip
    !First find the lsq fit line of slope FitSlope and intercept b, using the equations from Taylor (1982)
    Delta = bedlength*sum(thisbed(1,:)**2)-sum(thisbed(1,:))**2
    slope = -(bedlength*sum(thisbed(1,:)*thisbed(2,:))-sum(thisbed(1,:))*sum(thisbed(2,:)))/Delta !Negative because of the weird reversed slope sign convention I'm using.
    b = (sum(thisbed(1,:)**2)*sum(thisbed(2,:))-sum(thisbed(1,:))*sum(thisbed(1,:)*thisbed(2,:)))/Delta
case (4) !Fit to a known line
    b = intercepts(n)
    slope = RegSlope
case (5) !Fit for dip and/or intercept of the restored bed
    select case(fit_dips)
    case(0)
        slope = RegSlope
    case(1)
        slope = tan(restored_fit_params(1)*pi/180) !This is the same opposite sign convention slope used with the regular RegSlope and RegDip.
    case(2)
        slope = tan(restored_fit_params(n)*pi/180)
    case(3)
        do i = 1,n_bed_groups
            if (n>=bed_group_start(i) .and. (i==n_bed_groups .or. i<bed_group_start(i+1))) then
                group = i
            end if
        end do
        slope = tan(restored_fit_params(group)*pi/180)
    end select
    if (fit_beds==0) then !Bed intercepts aren't parameters.
        b = sum(thisbed(2,:)+slope*thisbed(1,:))/bedlength !Calculate y intercept of best fit line with slope as its slope
    else !Bed intercepts are parameters.
        select case(fit_dips)
        case(0)
            b = restored_fit_params(n)
        case(1)
            b = restored_fit_params(n+1)
        case(2)
            b = restored_fit_params(n+nbeds+DatTypes(DatType_dips))
        case(3)
            b = restored_fit_params(n+n_bed_groups+DatTypes(DatType_dips))
        end select
    end if
case (6) !Multi-segment line
    bed_dips = -restored_fit_params((1+(n-1)*2*n_bed_segs):(n_bed_segs+(n-1)*2*n_bed_segs))*pi/180. !Negative because of the sign convention we use for bed dips.
    bendxy(1,:) = restored_fit_params((1+n_bed_segs+(n-1)*2*n_bed_segs):(2*n_bed_segs-1+(n-1)*2*n_bed_segs))
    bendxy(2,1) = restored_fit_params(n*2*n_bed_segs)
    do i = 2,(n_bed_segs-1)
        bendxy(2,i) = bendxy(2,i-1)+(bendxy(1,i)-bendxy(1,i-1))*tan(bed_dips(i))
    end do
case default
    print*,'Error: Invalid FitType'
end select
!Note: In cases where slope or b is calculated, we are not currently calculating uncertainties in those values.
if (FitType == 1) then !Flat line
    R = abs(thisbed(2,:)-b) !y of each point minus average y value
    sigma = sqrt(Cbed(2,2,:)) !Since we only care about the uncertainties in y.
    if (ResultType == 5) sigma2 = sqrt(Cbed2(2,2,:))
else if (FitType /= 6) then !Anything that fits to a single straight line segment.
    R = abs(thisbed(2,:)+slope*thisbed(1,:)-b)/sqrt(slope**2+1) !errors (dist from pt to line)
    sigma = sqrt(((slope**2)*Cbed(1,1,:)+Cbed(2,2,:)+2*slope*Cbed(1,2,:))/(slope**2+1.)) !It would be -2*slope*Cbed(1,2,:) if we used the regular sign convention for slope.
    if (ResultType == 5) sigma2 = sqrt(((slope**2)*Cbed2(1,1,:)+Cbed2(2,2,:)+2*slope*Cbed2(1,2,:))/(slope**2+1.))
else !Fit to a multi-segment line
    do j=1,bedlength
        do i = 1,n_bed_segs
            if(i == n_bed_segs) then !Closest to the rightmost segment.
                R(j) = abs(((thisbed(1,j)-bendxy(1,i-1))*tan(bed_dips(i))-thisbed(2,j)+bendxy(2,i-1))/sqrt(1+tan(bed_dips(i))**2))
                exit
            else if (thisbed(1,j)<(thisbed(2,j)-bendxy(2,i))/tan((pi+bed_dips(i)+bed_dips(i+1))/2)+bendxy(1,i)) then !Closest to ramp i. (Tests whether the point is to the left of the barrier between segments i and i+1.)
                R(j) = abs(((thisbed(1,j)-bendxy(1,i))*tan(bed_dips(i))-thisbed(2,j)+bendxy(2,i))/sqrt(1+tan(bed_dips(i))**2))
                exit
            end if
        end do
        if (propagate == 1) then
            sigma(j) = sqrt(((tan(bed_dips(i))**2)*Cbed(1,1,j)+Cbed(2,2,j)-2*tan(bed_dips(i))*Cbed(1,2,j))/(tan(bed_dips(i))**2+1.))
            if (ResultType == 5) sigma2(j) = sqrt(((tan(bed_dips(i))**2)*Cbed2(1,1,j)+Cbed2(2,2,j)-2*tan(bed_dips(i))*Cbed2(1,2,j))&
                /(tan(bed_dips(i))**2+1.))
        end if
    end do
end if
select case(ResultType)
    case(1) !RMS
        call bed_ms(bedlength,R,beds_output(n))
    case(2) !Uncorrelated probability
        call bed_prob(bedlength,R,sigma,beds_output(n))
    case(3) !Chi square
        call bed_chisq(bedlength,thisbed,b,-slope,beds_output(n)) !Using -slope because of the backwards slope convention I've used elsewhere. Really need to change this.
    case(4) !Correlated probability
        if (state_for_Lc==2 .or. fit_Lc==1) then !Restored state, need to calculate R matrix now
            call make_bed_correlation_matrix(bedlength,thisbed,sigma,Lc,C)  !Note I have not tested this with propagated errors. I think it should work the same way, though.
        else !Deformed state, R matrix already calculated, but still need to calculate C matrix.
            S(:,:) = 0
            do i = 1,bedlength
                S(i,i) = sigma(i)
            end do
            C = matmul(matmul(S,spher_var_R(n)%pts),S) !Calculate the covariance matrix: C = S*R*S.
        end if
        call bed_prob_cov(bedlength,R,C,beds_output(n))
    case(5) !Correlated and uncorrelated probabilities.
        if (state_for_Lc==2 .or. fit_Lc==1) then !Restored state, need to calculate R matrix now
            call make_bed_correlation_matrix(bedlength,thisbed,sigma2,Lc,C)  !Note I have not tested this with propagated errors. I think it should work the same way, though.
        else !Deformed state, R matrix already calculated, but still need to calculate C matrix.
            S(:,:) = 0
            do i = 1,bedlength
                S(i,i) = sigma2(i)
            end do
            C = matmul(matmul(S,spher_var_R(n)%pts),S) !Calculate the covariance matrix: C = S*R*S.
        end if
        do i = 1,bedlength
            C(i,i) = C(i,i)+sigma(i)**2
        end do
        call bed_prob_cov(bedlength,R,C,beds_output(n))
    case default
        print*,'Error: Invalid ResultType'
end select
end subroutine calc_bed_err

subroutine calc_beds_err(beds_output_indiv,output_beds) !Combine the errors from all beds to calculate a final error. Works for terraces too.
use options, only: ResultType
implicit none
double precision,dimension(1:),intent(in) :: beds_output_indiv !Output from each of the beds individually.
double precision,intent(out) :: output_beds !Final output for all the beds put together.
select case(ResultType)
case(1) !RMS
    call beds_RMS(beds_output_indiv,output_beds)
case(2) !Uncorrelated probability
    call beds_prob(beds_output_indiv,output_beds)
case(3) !Chi square
    call beds_chisq(beds_output_indiv,output_beds)
case(4,5) !Correlated probability
    call beds_prob(beds_output_indiv,output_beds)
end select
end subroutine calc_beds_err

subroutine calc_restored_bed_err(n,slope,b,output,sigma_bed_restored_x,sigma_bed_restored_y,restored_fit_params) !Calculates the error associated with points expected to be on the restored bed.
use data_module, only: restored_beds
use options, only:ResultType,FitType,n_bed_segs
use constants, only: pi
implicit none
integer,intent(in) :: n !The number of the restored bed being considered.
double precision,intent(in) :: slope,b !The slope (negative sign convention I've been using) and intercept of the restored bed.
double precision,intent(out) :: output !The error or probability
double precision,intent(in) :: sigma_bed_restored_x,sigma_bed_restored_y
double precision,dimension(1:),optional :: restored_fit_params !Use if restored dips or intercepts are being fit as a parameter.
integer :: npts !The number of points in the bed being considered
double precision,dimension(:,:),allocatable :: pts !The points in the bed
double precision,dimension(:),allocatable :: R !Errors for each point.
double precision,dimension(:),allocatable :: sigma !Uncertainty for each point.
double precision,dimension(2,n_bed_segs-1) :: bendxy !x,y positions of bends in the restored bed.
double precision,dimension(n_bed_segs) :: bed_dips !Dips of the segments in a multi-segment restored bed.
integer :: i,j !counters
npts = restored_beds(n)%npts
allocate(pts(2,npts),R(npts),sigma(npts))
pts = restored_beds(n)%pts
if (FitType /= 6) then !Anything that fits to a single straight line segment.
    R = abs(pts(2,:)+slope*pts(1,:)-b)/sqrt(slope**2+1) !errors (dist from pt to line)
else !Fit to a multi-segment restored bed.
    bed_dips = -restored_fit_params(1:n_bed_segs)*pi/180. !Negative because of the sign convention we use for bed dips.
    bendxy(2,:) = restored_fit_params((n_bed_segs+1):(2*n_bed_segs-1))
    bendxy(1,1) = restored_fit_params(2*n_bed_segs)
    do i = 2,(n_bed_segs-1)
        bendxy(1,i) = bendxy(1,i-1)+(bendxy(2,i)-bendxy(2,i-1))/tan(bed_dips(i))
    end do
    do j=1,npts
        do i = 1,n_bed_segs-1
            if (pts(1,j)<(pts(2,j)-bendxy(2,i))/tan((pi+bed_dips(i)+bed_dips(i+1))/2)+bendxy(1,i)) then !Closest to segment i.
                R(j) = abs(((pts(1,j)-bendxy(1,i))*tan(bed_dips(i))-pts(2,j)+bendxy(2,i))/sqrt(1+tan(bed_dips(i))**2))
                exit
            else if(i == n_bed_segs-1) then !Closest to lower ramp
                R(j) = abs(((pts(1,j)-bendxy(1,i))*tan(bed_dips(n_bed_segs))-pts(2,j)+bendxy(2,i))/&
                    sqrt(1+tan(bed_dips(n_bed_segs))**2))
            end if
        end do
    end do
end if
select case(ResultType)
    case(2,4,5) !Probability, we're not using correlation for this.
        sigma = sqrt(((slope**2)*sigma_bed_restored_x**2+sigma_bed_restored_y**2)/(slope**2+1.)) !There shouldn't be any covariance in sigma x, sigma y because we haven't moved them and don't have an option for a priori covariance.
        call bed_prob(npts,R,sigma,output)
    case default
        print*,'Error: Invalid ResultType for Restored Bed Points' !Since restored bed points requires you to also have beds, only probability is allowed, since you have more than one data type.
end select
end subroutine calc_restored_bed_err

subroutine calc_dip_err(RestoredDips,output_dips,sigma_dip_array,sigma_restored_dip,restored_fit_params)
use options, only: FitType,ResultType,propagate,fit_dips,n_bed_segs,n_bed_groups
use parameters, only: RegDipDeg
use data_module, only: ndips,nbeds
use constants
implicit none
double precision,dimension(ndips),intent(in) :: RestoredDips !values of dips after inverse fault movement is complete
double precision, intent(out) :: output_dips !The error to be returned by the subroutine
double precision,dimension(ndips),optional :: sigma_dip_array !Use if we are using probability.
double precision,dimension(ndips) :: R !residuals between model and data
double precision,dimension(ndips) :: sigma !Array of sigmas for all values.
double precision,intent(in) :: sigma_restored_dip !Use if restored dips or intercepts are being fit as a parameter.
double precision,dimension(1:),optional :: restored_fit_params !Use if restored dips or intercepts are being fit as a parameter.
double precision :: MeanDip !Average restored dip
!Calculate the residuals
select case(FitType)
case (1) !Flat Line
    R = abs(RestoredDips) !Residuals, in degrees
case (2) !Known dip
    R = abs(RestoredDips-RegDipDeg) !Residuals, in degrees
case (3) !Fit for dip
    MeanDip  = sum(RestoredDips)/ndips
    R = abs(RestoredDips-MeanDip) !Residuals, in degrees
case (4) !Specified line. For dips this is the same as known dip.
    R = abs(RestoredDips-RegDipDeg) !Residuals, in degrees
case (5) !Fit for dip and/or intercept of the restored bed
    select case(fit_dips)
    case(0)
        R = abs(RestoredDips-RegDipDeg)
    case(1)
        R = abs(RestoredDips-restored_fit_params(1))
    case(2)
        R = abs(RestoredDips-restored_fit_params(nbeds+1))
    case(3)
        R = abs(RestoredDips-restored_fit_params(n_bed_groups+1))
    end select
case(6) !Multi-segment restored bed
    R = abs(RestoredDips-restored_fit_params(nbeds*2*n_bed_segs+1))
case default
    print*,'Error: Invalid FitType'
end select
where(R>90.) R=180-R !Always use the acute angle
!Calculate the error
sigma = sqrt(sigma_dip_array**2+sigma_restored_dip**2)
select case(ResultType)
case(1) !RMS
    call dips_RMS(R,output_dips)
case(2) !Uncorrelated probability
    call dips_prob(R,sigma,output_dips)
case(3) !Chi-square statistic
    call dips_chisq(R,RestoredDips,output_dips)
case(4,5) !Correlated probability for beds. Same as 2 for dips.
    call dips_prob(R,sigma,output_dips)
case default
    print*,'Error: Invalid ResultType'
end select
end subroutine

subroutine calc_terr_err(n,length,pts,terr_output_indiv,Cterr,sigma_orig_elev) !Calculate the error for an individual terrace.
use parameters, only: terr_slope,terr_orig_elev
use options, only: ResultType,propagate
use data_module, only: nterr
implicit none
integer, intent(in) :: n !The number of the terrace being worked on.
integer, intent(in) :: length !number of points in the terrace being worked on
double precision, dimension(2,length), intent(in) :: pts
double precision,dimension(nterr),intent(out) :: terr_output_indiv !output for each bed
double precision,dimension(2,2,length),optional :: Cterr !Use for probability. Covariance matrix for (x,y) for each point on the terrace.
double precision,dimension(length) :: sigma_orig_elev !Use for probability. Uncertainty in the original elevation of the points.
double precision,dimension(length) :: R !residuals between model and data
double precision,dimension(length,length) :: CR !Covariance matrix for uncertainty in R (residuals).
!double precision,dimension(2) :: inedge_expect !A point at the expected inner edge formation elevation, along a line perpendicular to the terrace slope from the restored inner edge
double precision :: b !intercept of best fit line
integer :: i !A counter of points.
!Calculate residuals for this bed
!inedge_expect = (/pts(1,1)+(terr_orig_elev(n)-pts(2,1))*terr_slope,terr_orig_elev(n)/)
!sigma_inedge_expect = (/sigma_orig_elev(n)*terr_slope,sigma_orig_elev(n)/) !Uncertainty in the expected inner edge position.
!b = inedge_expect(2) + terr_slope*inedge_expect(1) !y intercept of  line through the expected inner edge position with the terrace slope. + b/c of the backwards slope sign convention.
!sigma_b = sqrt(sigma_inedge_expect(2)**2+(terr_slope*sigma_inedge_expect(1))**2)
!There are 2 ways to estimate the expected inner edge position, and thus b: 
!1) Assume the inner edge is at the restored x position but at the expected y position. 
!2) Calculate a point at the expected y, but along a line perpendicular to the terrace slope through the restored x.
!I'm just using the first option, but the commented out code above is for the second.
!Probably I should also include the uncertainty in pts(1,1) to calculate a sigma_b from sigma_orig_elev. This may affect the covariance matrix for sigmaR, since it will give other points a covariance with point 1.
b = terr_orig_elev(n) + terr_slope(n)*pts(1,1) !y intercept of  line through the expected inner edge position with the terrace slope. + b/c of the backwards slope sign convention. pts(1,1) is the restored x value of the point. It is at some y value but should be at terr_orig_elev.
R = abs(pts(2,:)+terr_slope(n)*pts(1,:)-b)/sqrt(terr_slope(n)**2+1) !errors (dist from pt to line)
select case(ResultType)
    case(1) !RMS
        call bed_ms(length,R,terr_output_indiv(n))
    case(2,4,5) !Probability. Correlation (case 4) is only considered for beds at the moment, not for terrace points.
        !CR = ((terr_slope(n)*Cterr(1,1,1))**2+sigma_orig_elev(n)**2)/(1+terr_slope(n)**2) !For all CR(i,j), i/=j, i/=1, and j/=1.
        !CR(1,1) = (Cterr(2,2,1)**2+sigma_orig_elev(n)**2)/(1+terr_slope(n)**2) !No sigma_x because xi-x1 = 0 when i = 1
        CR = ((terr_slope(n)**2)*Cterr(1,1,1)+sigma_orig_elev(n)**2)/(1+terr_slope(n)**2) !For all CR(i,j), i/=j, i/=1, and j/=1.
        CR(1,1) = (Cterr(2,2,1)+sigma_orig_elev(n)**2)/(1+terr_slope(n)**2) !No sigma_x because xi-x1 = 0 when i = 1
        do i = 2,length
            CR(1,i) = (sigma_orig_elev(n)**2)/(1+terr_slope(n)**2)
            CR(i,1) = (sigma_orig_elev(n)**2)/(1+terr_slope(n)**2)
            !CR(i,i) = (Cterr(2,2,i)**2+(terr_slope(n)**2)*(Cterr(1,1,i)**2+Cterr(1,1,1)**2) &
            CR(i,i) = (Cterr(2,2,i)+(terr_slope(n)**2)*(Cterr(1,1,i)+Cterr(1,1,1)) &
                +sigma_orig_elev(n)**2+2.*terr_slope(n)*Cterr(1,2,i))/(1+terr_slope(n)**2) !Cterr(1,2,i) should == Cterr(2,1,i). If not there's a problem.
        end do
        call bed_prob_cov(length,R,CR,terr_output_indiv(n))
    case(3) !Chi square
        call bed_chisq(length,pts,b,-terr_slope(n),terr_output_indiv(n)) !Using -slope because of the backwards slope convention I've used elsewhere. Really need to change this.
    case default
        print*,'Error: Invalid ResultType'
end select
end subroutine calc_terr_err

subroutine calc_fault_pts_err(R,output_fault,sigma_fault)
use data_module, only: nfaultpts
use options, only:ResultType
implicit none
double precision,dimension(nfaultpts),intent(in) :: R !The errors for each point on the fault (residuals between model and data).
double precision,intent(in) :: sigma_fault !Use for probability. Uncertainty in the fault points.
double precision,dimension(nfaultpts) :: sigma !Array of sigmas for all values.
double precision :: output_fault !The final output (RMS, probability, or chi squared).
select case(ResultType)
    case(1) !RMS
        output_fault = sqrt(sum(R**2)/nfaultpts)
    case(2) !Uncorrelated probability
        sigma = sigma_fault !No need to propagate errors in fault points.
        call bed_prob(nfaultpts,R,sigma,output_fault)
    case(3) !Chi square
        print*,'Error: Chi square not available for points on fault.'
    case(4,5) !Correlated probability. The correlation applies only to beds, not fault points.
        sigma = sigma_fault !No need to propagate errors in fault points.
        call bed_prob(nfaultpts,R,sigma,output_fault)
    case default
        print*,'Error: Invalid ResultType'
end select

end subroutine calc_fault_pts_err

subroutine bed_ms(bedlength,R,ms_err) !mean squared error for a single bed
implicit none
integer,intent(in) :: bedlength !number of points in the bed being worked on
double precision,dimension(bedlength),intent(in) :: R!the bed being worked on
double precision,intent(out) :: ms_err !mean square error for the bed
ms_err = sum(R**2)/bedlength
end subroutine Bed_ms

subroutine beds_RMS(ms_errs,RMS) !RMS error for more than one bed together.
use data_module, only: beds
implicit none
double precision,dimension(1:),intent(in) :: ms_errs !mean square errors for each bed
double precision,intent(out) :: RMS !RMS error for this model
integer :: n !The number of beds
integer :: i !A counter
double precision :: sum_sq_errs !Total sum of squares of errors for all points in all beds.
integer :: npts_total !Total number of points.
n = size(ms_errs)
!The original code just did RMS = sqrt(sum(ms_errs)/n), i.e. averaged the averages. This is only right if all beds are the same 
!length (or if you want to treat them each as equally important). What I've put in here instead is a quick, hacked-in way to fix
!that. Eventually, I should make this nicer, so I'm not just undoing an earlier calculation of the mean.
sum_sq_errs = 0.;
npts_total = 0;
do i = 1,n
    sum_sq_errs = sum_sq_errs+ms_errs(i)*beds(i)%npts
    npts_total = npts_total+beds(i)%npts
end do
RMS = sqrt(sum_sq_errs/npts_total)
!RMS = sqrt(sum(ms_errs)/n)
end subroutine beds_RMS

subroutine bed_prob(bedlength,R,sigma,lnp) !Uncorrelated bed probability
implicit none
integer,intent(in) :: bedlength !number of points in the bed being worked on
double precision,dimension(bedlength),intent(in) :: R!the bed being worked on
double precision,dimension(bedlength),intent(in) :: sigma !The uncertainty in the bed positions
double precision,intent(out) :: lnp !The natural logarithm of the probability (unnormalized)
lnp = -sum((R**2)/(2.*sigma**2)+log(sigma))
end subroutine bed_prob

subroutine make_bed_correlation_matrix(bedlength,bedpts,sigma,Lc,C)
implicit none
integer,intent(in) :: bedlength !number of points in the bed being worked on
double precision,dimension(2,bedlength),intent(in) :: bedpts !the x,y points of the bed being worked on
double precision,intent(in),dimension(bedlength) :: sigma !The uncertainty in the bed positions
double precision,intent(in) :: Lc !correlation length
double precision,dimension(bedlength,bedlength),intent(out) :: C !Covariance matrix
double precision,dimension(bedlength,bedlength) :: R,S !Matrices used in calculating C.
integer :: i,j !counters
double precision :: Lij,h !Distance between points i and j, and h = Lij/Lc
!Make S matrix
S(:,:) = 0
do i = 1,bedlength
    S(i,i) = sigma(i)
end do
!Make R matrix
do i = 1,bedlength
    Lij = 0
    do j = i,bedlength
        if (j>i) then !Add the next increment to Lij
            Lij = Lij + sqrt((bedpts(1,j)-bedpts(1,j-1))**2+(bedpts(2,j)-bedpts(2,j-1))**2)
        end if
        h = Lij/Lc
        if (h < 1.) then
            R(i,j) = 1.+0.5*(h**3-3.*h)
        else
            R(i,j) = 0
        end if
        R(j,i) = R(i,j) !Since the length between the two points is the same either way.
    end do
end do
!Calculate the covariance matrix: C = S*R*S
C = matmul(matmul(S,R),S) !Note: I'm pretty sure C will be symmetrical
end subroutine make_bed_correlation_matrix

subroutine make_spherical_variogram_R(beds_in,Lc,R_all)
!Takes in an array of beds and outputs an array of spherical variogram R matrices, formatted as an array of bed data types.
use data_module, only: bed
implicit none
type(bed),dimension(1:),intent(in) :: beds_in
double precision,intent(in) :: Lc !correlation length
type(bed),dimension(:),allocatable,intent(out) :: R_all !Array of all the R matrices. It's not really a bed, but the bed data type should work for this.
integer :: n,i,j !Counters
double precision :: Lij,h !Distance between points i and j, and h = Lij/Lc
integer :: nbeds !The number of beds in the beds_in array.
integer :: bedlength !number of points in the bed being worked on
double precision,dimension(:,:),allocatable :: bedpts !the x,y points of the bed being worked on
double precision,dimension(:,:),allocatable :: R !R matrix for a single bed.
type(bed) :: Rbed !R matrix as a bed.
nbeds = size(beds_in,1)
allocate(R_all(nbeds))
do n = 1,nbeds
    bedlength = beds_in(n)%npts
    allocate(bedpts(2,bedlength))
    allocate(R(bedlength,bedlength))
    bedpts = beds_in(n)%pts
    do i = 1,bedlength
        Lij = 0
        do j = i,bedlength
            if (j>i) then !Add the next increment to Lij
                Lij = Lij + sqrt((bedpts(1,j)-bedpts(1,j-1))**2+(bedpts(2,j)-bedpts(2,j-1))**2)
            end if
            h = Lij/Lc
            if (h < 1.) then
                R(i,j) = 1.+0.5*(h**3-3.*h)
            else
                R(i,j) = 0
            end if
            R(j,i) = R(i,j) !Since the length between the two points is the same either way.
        end do
    end do
    Rbed = bed(R,bedlength,beds_in(n)%ident,.false.)
    R_all(n) = Rbed
    deallocate(bedpts,R)
end do
end subroutine make_spherical_variogram_R

subroutine bed_prob_cov(bedlength,R,C,lnp) !Bed probability when there is a covariance matrix.
use math, only: matinv,chol
implicit none
integer,intent(in) :: bedlength !number of points in the bed being worked on
double precision,dimension(bedlength,1),intent(in) :: R!the errors in the bed being worked on
double precision,dimension(bedlength,bedlength) :: C,Cinv,L !The covariance matrix for R, its inverse, and its Cholesky decomposition.
double precision :: logsqrtdetC !The square root of the determinant of C
double precision,intent(out) :: lnp !The ln of the probability (unnormalized)
double precision,dimension(1,1) :: lnparray !ln(p) in the form of a 1x1 array
integer :: i !A counter of bed points.
Cinv = matinv(bedlength,C)
L = chol(bedlength,C)
logsqrtdetC = 0; !Initialize detC to 1 and then multiply by L(i,i)^2 is equivalent to initialize logsqrtdetC to 0 and then add log(L(i,i))
do i = 1,bedlength !I'm pretty sure this is right. Since C = LL^T, and diagonals of L and L^T are equal, det(C)=det(L)^2. And determinant of L = product of diagonal entries since L is triangular.
    logsqrtdetC = logsqrtdetC+log(L(i,i)) !This is log of sqrt of the determinant to keep it from growing huge and because that's what we use in the equation for probability anyway.
end do
lnparray = -0.5*matmul(matmul(transpose(R),Cinv),R)-logsqrtdetC
lnp = lnparray(1,1)
end subroutine bed_prob_cov

subroutine beds_prob(lnpbeds,lnpall) !Uncorrelated probability for multiple beds
implicit none
double precision,dimension(1:),intent(in) :: lnpbeds !vector of ln of probabilities for each bed
double precision,intent(out) :: lnpall !combined ln of probability of all the beds
lnpall = sum(lnpbeds) !ln(a*b) = ln(a)+ln(b)
end subroutine beds_prob

subroutine dips_RMS(R,RMS) !RMS error for dips
use data_module, only:ndips
implicit none
double precision,dimension(ndips),intent(in) :: R !Residuals
double precision,intent(out) ::RMS !RMS error for this model
RMS = sqrt(sum(R**2)/ndips) !Calculate the RMS
end subroutine dips_RMS

subroutine dips_prob(R,sigma,lnp) !Uncorrelated probability for dips
use data_module, only:ndips
implicit none
double precision,dimension(ndips),intent(in) :: R !Residuals
double precision,dimension(ndips),intent(in) :: sigma !The uncertainty in the dips
double precision,intent(out) :: lnp !The combined natural logarithm of probability for all the dip data
lnp = -sum((R**2)/(2*sigma**2)+log(sigma))
end subroutine dips_prob

subroutine bed_chisq(bedlength,bed,b,slope,chisq) !Chi square statistic for beds.
implicit none
integer,intent(in) :: bedlength !number of points in the bed being worked on
double precision,dimension(2,bedlength),intent(in) :: bed !this bed
double precision,intent(in) :: b !intercept of best fit line
double precision,intent(in) :: slope !slope of best fit line
double precision,intent(out) :: chisq ! The chi squared statistic for this bed.
double precision,dimension(bedlength) :: xm,ym !The closest x and y positions on the model (best fit) line to each point on the bed.
if (slope .eq. 0) then !Special case of when slope is zero, so we don't divide by zero.
    xm = bed(1,:)
    ym = b
else
    !xm = ((-slope)*bed(1,:)-bed(2,:))/(2*slope)
    !ym = slope*xm+b
    xm = ((bed(1,:)+slope*bed(2,:))-slope*b)/(slope**2+1)
    ym = (slope*(bed(1,:)+slope*bed(2,:))+b)/(slope**2+1)
end if
chisq = sum(((bed(1,:)-xm)**2)/abs(bed(1,:)+xm)) + sum(((bed(2,:)-ym)**2)/abs((bed(2,:)+ym))) !Calculate the chisq statistic
end subroutine bed_chisq

subroutine beds_chisq(chisq_beds,chisq_all) !Chi square for multiple beds
implicit none
double precision,dimension(1:),intent(in) :: chisq_beds !vector of chisq for each bed
double precision,intent(out) :: chisq_all !combined chisq of all the beds
chisq_all = sum(chisq_beds)
end subroutine beds_chisq

subroutine dips_chisq(R,RestoredDips,chisq_dips) 
use data_module, only:ndips
use parameters, only:RegDipDeg
use options, only:FitType
implicit none
double precision,dimension(ndips),intent(in) :: R !Residuals
double precision,dimension(ndips),intent(in) :: RestoredDips !The restored dips
double precision,intent(out) :: chisq_dips !The combined chi-square statistic  for all the dip data
double precision,dimension(ndips) :: ModelDips !Modeled dips
select case(FitType)
case (1) !Flat Line
    ModelDips = 0
case (2) !Known dip
    ModelDips = RegDipDeg
case (3) !Fit for dip
    ModelDips  = sum(RestoredDips)/ndips
case (4) !Specified line. For dips this is the same as known dip.
    ModelDips = RegDipDeg
case default
    print*,'Error: Invalid FitType'
end select
chisq_dips = sum((R**2)/(abs(RestoredDips+ModelDips+tiny(0.)))) !Since ModelDips = RestoredDips + R
end subroutine dips_chisq

subroutine prop_bend_error(C,v1,v2,vb,slope)
!This function propagates the uncertainty in the position of a point as it moves through a boundary between two constant velocity domains.
!Everything with dimension(2) is x coordinate or x component followed by y coordinate or y component.
double precision,dimension(2,2),intent(inout) :: C !The covariance matrix for x and y coordinates of the point.
double precision,dimension(2),intent(in) :: v1 !Velocity in the domain the point is coming from.
double precision,dimension(2),intent(in) :: v2 !Velocity in the domain the point is moving into.
double precision,dimension(2),intent(in) :: vb !Velocity at which the domain boundary is moving
double precision,intent(in) :: slope !The slope of the boundary. This is not allowed to change. Equals tan(boundary dip).
double precision,dimension(2,2) :: J !The Jacobian matrix.
J(1,1) = 1-slope*((v1(1)-v2(1))/(slope*(v1(1)-vb(1))-(v1(2)-vb(2)))) !dx/dx0
J(1,2) = (v1(1)-v2(1))/(slope*(v1(1)-vb(1))-(v1(2)-vb(2))) !dx/dy0
J(2,1) = -slope*((v1(2)-v2(2))/(slope*(v1(1)-vb(1))-(v1(2)-vb(2)))) !dy/dx0
J(2,2) = 1+((v1(2)-v2(2))/(slope*(v1(1)-vb(1))-(v1(2)-vb(2)))) !dy/dy0
C = matmul(matmul(J,C),transpose(J))
end subroutine prop_bend_error

end module err_and_prob
