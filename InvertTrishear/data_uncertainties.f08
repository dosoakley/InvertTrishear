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

module data_uncertainties
!This module contains variables and subroutines related to the uncertainty in data of different types.
!Note: I may want to rearrange things in the future so that different data types each have their own module, in which case I may want to split the uncertainty information up into those modules.

!Parameters relating to uncertainties:
double precision :: sigma_bed,sigma2_bed,Lc_fixed,sigma_dip_fixed,sigma_fault_fixed !For uncertainties. These too will be removed from this module when different result types are separated out.
double precision,dimension(:),allocatable :: sigma_bed_x_fixed,sigma_bed_y_fixed,sigma2_bed_x_fixed,sigma2_bed_y_fixed !x and y uncertainties in the positions of points in a bed when the uncertainty is fixed and not a model parameter.
double precision,dimension(:),allocatable :: sigma_terr_fixed !Uncertainty in terrace points. This can be different for different terraces.
double precision,dimension(:),allocatable :: sigma_orig_elev_fixed !Uncertainty in the restored inner edge elevation(s) of terraces (i.e. uncertainty in paleosealevel).
double precision :: sigma_restored_dip_fixed !Uncertainty in the restored dip values, that is not propagated through.
double precision :: sigma_bed_restored,sigma_bed_restored_x_fixed,sigma_bed_restored_y_fixed !Uncertainties in points expected to be on restored beds.
integer :: nparam_start_sigma_bed,nparam_start_sigma2_bed,nparam_start_sigma_dip,nparam_start_sigma_terr,nparam_start_sigma_fault !Parameter numbers of the first parameter with sigma values for each data type.
integer :: nparam_start_sigma_bed_restored  !Parameter numbers of the first parameter with sigma values for each data type.
integer :: nparam_Lc !Parameter number of the correlation length (Lc).

contains

subroutine GetSigmaBed(params,sigma_bed_x,sigma_bed_y)
    use options, only: fit_sigma,sigma_xy_diff,diff_sigma_bed,ResultType,n_sigma_bed_groups,sigma_bed_group_start
    use data_module, only: nbeds
    implicit none
    double precision,dimension(1:),intent(in) :: params
    double precision,dimension(1:),intent(out) :: sigma_bed_x,sigma_bed_y !x and y uncertainties in bed points.
    integer :: i
    integer :: end_ind !Index of the end of a bed group.
    if (ResultType==2 .or. ResultType>=4) then
        if (fit_sigma == 0) then
            sigma_bed_x = sigma_bed_x_fixed
            sigma_bed_y = sigma_bed_y_fixed
        else
            if (sigma_xy_diff ==0) then
                if (diff_sigma_bed == 1) then !Different uncertainties for each bed.
                    sigma_bed_x = params(nparam_start_sigma_bed:nparam_start_sigma_bed+nbeds-1)
                else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                    do i = 1,n_sigma_bed_groups
                        if (i<n_sigma_bed_groups) then
                            end_ind = sigma_bed_group_start(i+1)-1
                        else
                            end_ind = nbeds
                        end if
                        sigma_bed_x(sigma_bed_group_start(i):end_ind) = params(nparam_start_sigma_bed+i-1)
                    end do
                else
                    sigma_bed_x = params(nparam_start_sigma_bed)
                end if
                sigma_bed_y = sigma_bed_x
            else if(sigma_xy_diff == 1) then !Different x and y uncertainties.
                if (diff_sigma_bed == 1) then  !Different uncertainties for each bed.
                    sigma_bed_x = params(nparam_start_sigma_bed:nparam_start_sigma_bed+nbeds-1)
                    sigma_bed_y = params(nparam_start_sigma_bed+nbeds:nparam_start_sigma_bed+2*nbeds-1)
                else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                    do i = 1,n_sigma_bed_groups
                        if (i<n_sigma_bed_groups) then
                            end_ind = sigma_bed_group_start(i+1)-1
                        else
                            end_ind = nbeds
                        end if
                        sigma_bed_x(sigma_bed_group_start(i):end_ind) = params(nparam_start_sigma_bed+i-1)
                        sigma_bed_y(sigma_bed_group_start(i):end_ind) = params(nparam_start_sigma_bed+n_sigma_bed_groups+i-1)
                    end do
                else
                    sigma_bed_x = params(nparam_start_sigma_bed)
                    sigma_bed_y = params(nparam_start_sigma_bed+1)
                end if
            end if
        end if
    else
        sigma_bed_x = 0
        sigma_bed_y = 0
    end if
end subroutine GetSigmaBed

subroutine GetSigma2Bed(params,sigma2_bed_x,sigma2_bed_y)
    !2nd sigma, for use if there are both correlated and uncorrelated errors.
    use options, only: fit_sigma,sigma_xy_diff,diff_sigma_bed,ResultType,n_sigma_bed_groups,sigma_bed_group_start
    use data_module, only: nbeds
    implicit none
    double precision,dimension(1:),intent(in) :: params
    double precision,dimension(1:),intent(out) :: sigma2_bed_x,sigma2_bed_y !x and y uncertainties in bed points.
    integer :: i
    integer :: end_ind !Index of the end of a bed group.
    if (ResultType==5) then
        if (fit_sigma == 0) then
            sigma2_bed_x = sigma2_bed_x_fixed
            sigma2_bed_y = sigma2_bed_y_fixed
        else
            if (sigma_xy_diff ==0) then
                if (diff_sigma_bed == 1) then !Different uncertainties for each bed.
                    sigma2_bed_x = params(nparam_start_sigma2_bed:nparam_start_sigma2_bed+nbeds-1)
                else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                    do i = 1,n_sigma_bed_groups
                        if (i<n_sigma_bed_groups) then
                            end_ind = sigma_bed_group_start(i+1)-1
                        else
                            end_ind = nbeds
                        end if
                        sigma2_bed_x(sigma_bed_group_start(i):end_ind) = params(nparam_start_sigma2_bed+i-1)
                    end do
                else
                    sigma2_bed_x = params(nparam_start_sigma2_bed)
                end if
                sigma2_bed_y = sigma2_bed_x
            else if(sigma_xy_diff == 1) then !Different x and y uncertainties.
                if (diff_sigma_bed == 1) then  !Different uncertainties for each bed.
                    sigma2_bed_x = params(nparam_start_sigma2_bed:nparam_start_sigma2_bed+nbeds-1)
                    sigma2_bed_y = params(nparam_start_sigma2_bed+nbeds:nparam_start_sigma2_bed+2*nbeds-1)
                else if (diff_sigma_bed==2) then !Different uncertainties for bed groups.
                    do i = 1,n_sigma_bed_groups
                        if (i<n_sigma_bed_groups) then
                            end_ind = sigma_bed_group_start(i+1)-1
                        else
                            end_ind = nbeds
                        end if
                        sigma2_bed_x(sigma_bed_group_start(i):end_ind) = params(nparam_start_sigma2_bed+i-1)
                        sigma2_bed_y(sigma_bed_group_start(i):end_ind) = params(nparam_start_sigma2_bed+n_sigma_bed_groups+i-1)
                    end do
                else
                    sigma2_bed_x = params(nparam_start_sigma2_bed)
                    sigma2_bed_y = params(nparam_start_sigma2_bed+1)
                end if
            end if
        end if
    else
        print*,'Error: GetSigma2Bed should only be used if ResultType==5.'
        sigma2_bed_x = 0
        sigma2_bed_y = 0
    end if
end subroutine GetSigma2Bed

subroutine GetLc(params,Lc)
!Get the correlation length from the parameters list if using correlated errors and fitting for the correlation length.
use options, only: ResultType,fit_Lc
implicit none
double precision,dimension(1:),intent(in) :: params
double precision,intent(out) :: Lc
if (ResultType >= 4) then !Only in this case are we using Lc at all.
    if (fit_Lc == 0) then
        Lc = Lc_fixed
    else
        Lc = params(nparam_Lc)
    end if
else
    Lc = 0
end if
    
end subroutine

subroutine GetSigmaDip(params,sigma_dip,sigma_restored_dip)
    use options,only: fit_sigma,use_sigma_restored_dip,ResultType
    implicit none
    double precision,dimension(1:),intent(in) :: params
    double precision,intent(out) :: sigma_dip,sigma_restored_dip
    if (ResultType==2 .or. ResultType==4) then
        if (fit_sigma == 0) then
            sigma_dip = sigma_dip_fixed
            sigma_restored_dip = sigma_restored_dip_fixed
        else
            sigma_dip = params(nparam_start_sigma_dip)
            if (use_sigma_restored_dip == 1) then
                sigma_restored_dip = params(nparam_start_sigma_dip+1)
            else
                sigma_restored_dip = sigma_restored_dip_fixed
            end if
        end if
    else
        sigma_dip = 0
        sigma_restored_dip = 0
    end if
end subroutine GetSigmaDip

subroutine GetSigmaTerr(params,sigma_terr,sigma_orig_elev)
    use options,only: fit_sigma,ResultType
    use data_module, only: nterr
    implicit none
    double precision,dimension(1:),intent(in) :: params
    double precision,dimension(nterr),intent(out) :: sigma_terr,sigma_orig_elev
    if (ResultType==2 .or. ResultType==4) then
        if (fit_sigma == 0) then
            sigma_terr = sigma_terr_fixed
            sigma_orig_elev = sigma_orig_elev_fixed
        else
            sigma_terr = params(nparam_start_sigma_terr:nparam_start_sigma_terr+nterr-1)
            sigma_orig_elev = params(nparam_start_sigma_terr+nterr:nparam_start_sigma_terr+2*nterr-1)
        end if
    else
        sigma_terr = 0
        sigma_orig_elev = 0
    end if
end subroutine GetSigmaTerr

subroutine GetSigmaFault(params,sigma_fault)
    use options,only: fit_sigma,ResultType
    implicit none
    double precision,dimension(1:),intent(in) :: params
    double precision,intent(out) :: sigma_fault
    if (ResultType==2 .or. ResultType==4) then
        if (fit_sigma == 0) then
            sigma_fault = sigma_fault_fixed
        else
            sigma_fault = params(nparam_start_sigma_fault)
        end if
    else
        sigma_fault = 0
    end if
end subroutine GetSigmafault

subroutine GetSigmaBedRestored(params,sigma_bed_restored_x,sigma_bed_restored_y)
    use options, only: fit_sigma,sigma_xy_diff,ResultType
    implicit none
    double precision,dimension(1:),intent(in) :: params
    double precision,intent(out) :: sigma_bed_restored_x,sigma_bed_restored_y !x and y uncertainties in bed points.
    if (ResultType==2 .or. ResultType==4) then
        if (fit_sigma == 0) then
            sigma_bed_restored_x = sigma_bed_restored_x_fixed
            sigma_bed_restored_y = sigma_bed_restored_y_fixed
        else
            if (sigma_xy_diff==0) then
                sigma_bed_restored_x = params(nparam_start_sigma_bed_restored)
                sigma_bed_restored_y = sigma_bed_restored_x
            else
                sigma_bed_restored_x = params(nparam_start_sigma_bed_restored)
                sigma_bed_restored_y = params(nparam_start_sigma_bed_restored+1)
            end if
        end if
    else
        sigma_bed_restored_x = 0
        sigma_bed_restored_y = 0
    end if
end subroutine GetSigmaBedRestored

end module data_uncertainties