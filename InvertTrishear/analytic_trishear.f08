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

module analytic_trishear_module
!Functions and subroutines for the analytic / semi-analytic trishear solution in polar coordinates.
contains

!Note: These functions currently fail (produce NaNs) in the special cases P/S = 1 and P/S = 0. Specific code is needed to address this.
!Also, there is some kind of a problem occurring with the program freezing up when u is at or near m or -m.
!The problem seems to involve one or both of the functions repeating itself indefinitely. Adding return statements might help some, but not enough.
!In particular, I think the trishear_slip function is causing problems.

subroutine analytic_trishear(x,y,C,m,PoverS,rem_slip,slip_sign,slip)
    !use options, only: errs_id
    implicit none
    double precision,intent(inout) :: x,y !Position in trishear coordinates
    double precision,dimension(2,2),intent(inout) :: C !Covariance matrix for x and y
    double precision,intent(in) :: m,PoverS !m = tan(phi), and PoverS = fault propagation to slip ratio
    double precision,intent(in) :: rem_slip !Remaining maximum possible slip
    integer,intent(in) :: slip_sign !Tells if slip is forwards (positive) or backwards (negative). Should be +1 or -1.
    double precision,intent(out) :: slip !Amount of slip the point undergoes in the trishear zone.
    double precision :: r0,u0,r,u !New coordinates, initial and final positions.
    double precision :: hw_slip,fw_slip !Slip necessary to reach the hanging wall or footwall
    double precision :: c1 !A constant
    double precision :: W,W0 !A term that gets used a lot, to simplify equations.
    double precision,dimension(2,2) :: J !The Jacobian matrix
    double precision :: dc1du0,dSdu0,dudS,dudu0,dudr0,drdu_changingu0,drdu_changingr0,drdu0,drdr0 !Various derivates needed for error propagation
    double precision :: a,b,c3
    !It might help to also output something to tell if the point ends up in the fw, hw, or still in the trishear zone. I may able implement that later if it seems necessary.
    !Transform coordinates
    r0 = sqrt(x**2+y**2)
    u0 = y/x
    if (u0>m .and. u0-m<1e-6) then !This corrects for some rounding errors that can occur. I'm not sure if I need the <1e-6 part in this or the -m version.
        u0 = m
        !    u0 = m-1e-1 !1e-1 is too large an amount to add/subtract
    else if(u0<-m .and. -u0-m<1e-6) then
        u0 = -m
        !    u0 = -m+1e-1
    end if
    J = reshape((/x/r0,-y/x**2,y/r0,1/x/),(/2,2/)) !dr/dx, du/dx, dr/dy,du/dy, corresponding to J11,J21,J12,J22 respectively in column major order.
    C = matmul(matmul(J,C),transpose(J))
    !Calculate new r, u and slip.
    u = u0
    r = r0
    W0 = m**2+2*m*(1-2*PoverS)*u0+u0**2 !This term gets used a lot
    c1 = r0*W0/sqrt(1+u0**2)
    if ((slip_sign == -1 .and. PoverS<1) .or. (slip_sign == 1 .and. PoverS>1)) then !If the point can leave the trishear zone.
        !print*,'a'
        hw_slip = trishear_slip(u,m,m,PoverS,c1)
        if (slip_sign*hw_slip>0 .and. abs(hw_slip)<=abs(rem_slip)) then
            !print*,'b'
            slip = hw_slip
            u = m
        else
            !print*,'c'
            fw_slip = trishear_slip(u,-m,m,PoverS,c1)
            if (slip_sign==1 .and. slip_sign*fw_slip>0 .and. abs(fw_slip)<=abs(rem_slip)) then !I'm pretty sure you can only enter the footwall if slip is forward. Velocity at footwall boundary is vx=-v0*PoverS, vy = 0, so it's entirely about the sign of v0.
                !print*,'d'
                slip = fw_slip
                u = -m
            else
                !print*,'e'
                u = trishear_u(u,rem_slip,m,PoverS,c1)
                slip = rem_slip
            end if
        end if
        !print*,'d'
    else
!        write(errs_id,*),'a'
!        write(errs_id,*),u,u0,rem_slip,m,PoverS,c1,r
!        a = 1. !For an equation of the form 1/(a*x^2+b*x*c)^2.
!        write(errs_id,*),a
!        b = 2*m*(1-2*PoverS)
!        write(errs_id,*),b
!        c3 = m**2
!        write(errs_id,*),c3
!        write(errs_id,*),b**2-4*a*c3
!        write(errs_id,*),(1/sqrt(b**2-4*a*c3))*log(abs((2*a*u0+b-sqrt(b**2-4*a*c3))/(2*a*u0+b+sqrt(b**2-4*a*c3))))
!        write(errs_id,*),(2/sqrt(4*a*c3-b**2))*atan((2*a*u0+b)/sqrt(4*a*c3-b**2))
        u = trishear_u(u,rem_slip,m,PoverS,c1)
!        write(errs_id,*),'b'
        slip = rem_slip
    end if
    W = m**2+2*m*(1-2*PoverS)*u+u**2 !This term gets used a lot
    r = c1*sqrt(1+u**2)/W
    !Propagate error in the location of the point
    dc1du0 = 2*r0*(m*(1-2*PoverS)+u0)/sqrt(1+u0**2)-u0*c1/(1+u0**2)
    dSdu0 = -(slip/c1)*dc1du0-4*m*c1*(u0**2+2*m*(1-2*PoverS)*u0+m**2)**(-2) !Not sure if this should have a negative sign (on both terms), but giving it a try. It seems to be necessary.
    !dSdu0 = dc1du0*(slip/c1+4*m*(u0**2+2*m*(1-2*PoverS)*u0+m**2)**(-2))
    dudS = (-1/(4*m*c1))*(u**2+2*m*(1-2*PoverS)*u+m**2)**2
    dudu0 = dudS*dSdu0
    dudr0 = -(slip/r0)*dudS !This, like dSdu0, seems to need a mysterious negative sign. Since both involve dS/dc1, I think they're related.
    !drdu = c1*u/(sqrt(1+u**2)*W)-c1*sqrt(1+u**2)*(2*u+2*m*(1-2*PoverS))/W**2+(dc1du0/dudu0)*sqrt(1+u**2)/W+(1/dudr0)*W0/W
    drdu_changingu0 = c1*u/(sqrt(1+u**2)*W)-c1*sqrt(1+u**2)*(2*u+2*m*(1-2*PoverS))/W**2+(dc1du0/dudu0)*sqrt(1+u**2)/W; !We need the last term on this because we're calculating dr/du for a given slip, rather than for a u0
    drdu0 = drdu_changingu0*dudu0
    !drdu0 = drdu*dudu0
    drdu_changingr0 = c1*u/(sqrt(1+u**2)*W)-c1*sqrt(1+u**2)*(2*u+2*m*(1-2*PoverS))/W**2+(1/dudr0)*W0*sqrt(1+u**2)/(W*sqrt(1+u0**2)); !The last term is dc1/du with u0 fixed and r0 a function of u
    drdr0 = drdu_changingr0*dudr0 !For this and drdu0, I could consolidate the drduc_changingr0 (or u0) and this into a single line and leave out the drdu terms as separate variables.
    !drdr0 = drdu*dudr0
    J = reshape((/drdr0,dudr0,drdu0,dudu0/),(/2,2/)) !The Jacobian matrix.
    C = matmul(matmul(J,C),transpose(J))
    !Transform back to regular trishear coordinates.
    x = r/sqrt(u**2+1);
    y = r*u/sqrt(u**2+1);
    if (u==m) then !Rounding errors can result in y not quite equalling +/-m*x when a point should be on the trishear boundary, which can cause problems with later tests for whether it's in the hw or fw.
        y = m*x
    else if (u==-m) then
        y = -m*x
    end if
    J = (1./sqrt(u**2+1.))*reshape((/1.0d0,u,-r*u/(u**2+1),r-(r*u**2)/(u**2+1.)/),(/2,2/)) !dx/dr, dy/dr, dx/du, dy/du, corresponding to J11,J21,J12,J22 respectively in column major order.
    C = matmul(matmul(J,C),transpose(J))
end subroutine analytic_trishear

subroutine analytic_trishear_dip(x,y,dip,sigmadip,m,PoverS,rem_slip,slip_sign,slip)
    implicit none
    double precision,intent(inout) :: x,y !Position in trishear coordinates
    double precision,intent(inout) :: dip !The dip at the point.
    double precision,intent(inout) :: sigmadip !Uncertainty in the dip. So far we are ignoring uncertainty in the points.
    double precision,intent(in) :: m,PoverS !m = tan(phi), and PoverS = fault propagation to slip ratio
    double precision,intent(in) :: rem_slip !Remaining maximum possible slip
    integer,intent(in) :: slip_sign !Tells if slip is forwards (positive) or backwards (negative). Should be +1 or -1.
    double precision,intent(out) :: slip !Amount of slip the point undergoes in the trishear zone.
    double precision :: r0,u0,r,u !New coordinates, initial and final positions.
    double precision :: dip0 !Initial dip
    double precision :: hw_slip,fw_slip !Slip necessary to reach the hanging wall or footwall
    double precision :: c1 !A constant
    double precision :: W,W0 !A term that gets used a lot, to simplify equations.
    double precision :: dr0du0,dc1du0,dSdu,du0dS,dc1du,drdu !Derivatives needed to calculate dip
    double precision :: ddr0du0ddip0,ddc1du0ddip0,ddu0dSddip0,ddc1duddip0,ddrduddip0,ddipddip0 !Various derivates needed for error propagation
    !It might help to also output something to tell if the point ends up in the fw, hw, or still in the trishear zone. I may able implement that later if it seems necessary.
    !Transform coordinates
    r0 = sqrt(x**2+y**2)
    u0 = y/x
    if (u0>m .and. u0-m<1e-6) then !This corrects for some rounding errors that can occur.
        u0 = m
        !    u0 = m-1e-1 !1e-1 is too large an amount to add/subtract
    else if(u0<-m .and. -u0-m<1e-6) then
        u0 = -m
        !    u0 = -m+1e-1
    end if
    !Calculate new r, u and slip.
    u = u0
    r = r0
    c1 = r*(m**2+2*m*(1-2*PoverS)*u+u**2)/sqrt(1+u**2)
    if ((slip_sign == -1 .and. PoverS<1) .or. (slip_sign == 1 .and. PoverS>1)) then !If the point can leave the trishear zone.
        hw_slip = trishear_slip(u,m,m,PoverS,c1)
        if (slip_sign*hw_slip>0 .and. abs(hw_slip)<=abs(rem_slip)) then
            slip = hw_slip
            u = m
        else
            fw_slip = trishear_slip(u,-m,m,PoverS,c1)
            if (slip_sign==1 .and. slip_sign*fw_slip>0 .and. abs(fw_slip)<=abs(rem_slip)) then !I'm pretty sure you can only enter the footwall if slip is forward. Velocity at footwall boundary is vx=-v0*PoverS, vy = 0, so it's entirely about the sign of v0.
                slip = fw_slip
                u = -m
            else
                u = trishear_u(u,rem_slip,m,PoverS,c1)
                slip = rem_slip
            end if
        end if
    else
        u = trishear_u(u,rem_slip,m,PoverS,c1)
        slip = rem_slip
    end if
    W = m**2+2*m*(1-2*PoverS)*u+u**2 !This term gets used a lot
    W0 = m**2+2*m*(1-2*PoverS)*u0+u0**2
    r = c1*sqrt(1+u**2)/W
    !Calculate dip
    dip0 = -dip !Because of the sign convention I've been using for the dips.
    dr0du0 = r0/(tan(dip0)-u0)+r0*u0/(u0**2+1);
    dc1du0 = dr0du0*W0/sqrt(1+u0**2)+r0*(2*m*(1-2*PoverS)+2*u0)/sqrt(1+u0**2)-r0*W0*u0/(1+u0**2)**(3./2.);
    dSdu = -4*m*c1*W**(-2);
    du0dS = -((slip/c1)*dc1du0+4*m*c1*W0**(-2))**(-1) !This seems to need a negative sign, but I'm not sure why. 
    dc1du = dc1du0*du0dS*dSdu;
    drdu = dc1du*sqrt(1+u**2)/W+c1*u/(sqrt(1+u**2)*W)-c1*sqrt(1+u**2)*(2*u+2*m*(1-2*PoverS))/W**2; !This is dr/du for known u0 and slip and unknown dip, with r0 as a function of u0 and dip.
    dip = atan((u*drdu+r-(r*u**2)/(u**2+1))/(drdu-r*u/(u**2+1)));
    !Propagate error in the dip.
    ddr0du0ddip0 = -r0/(((tan(dip0)-u0)**2)*(cos(dip0))**2)
    ddc1du0ddip0 = ddr0du0ddip0*W0/sqrt(1+u0**2)
    ddu0dSddip0 = (slip/c1)*ddc1du0ddip0/((slip/c1)*dc1du0+4*m*c1*W0**(-2))**2 !Changed this to + instead of - in front, because du0dS had - instead of + in front.
    ddc1duddip0 = ddc1du0ddip0*du0dS*dSdu+dc1du0*ddu0dSddip0*dSdu
    ddrduddip0 = ddc1duddip0*sqrt(1+u**2)/W
    ddipddip0 = (1/(tan(dip)**2+1))*(-r/(drdu-r*u/(u**2+1))**2)*ddrduddip0
    sigmadip = ddipddip0*sigmadip
    !Convert back to regular trishear coordinates
    x = r/sqrt(u**2+1);
    y = r*u/sqrt(u**2+1);
    if (u==m) then !Rounding errors can result in y not quite equalling +/-m*x when a point should be on the trishear boundary, which can cause problems with later tests for whether it's in the hw or fw.
        y = m*x
    else if (u==-m) then
        y = -m*x
    end if
    dip = -dip !Back to the sign convention that the rest of the program uses.
end subroutine analytic_trishear_dip

function trishear_slip(u0,ufinal,m,PoverS,c1)
!use options, only: errs_id
!Determines the slip necessary to get to some ufinal.
implicit none
double precision :: u0,ufinal,m,PoverS,c1
double precision :: a,b,c,c2,part1,part1_0
double precision :: trishear_slip
!write(errs_id,*),'a'
!print*,'abc'
a = 1. !For an equation of the form 1/(a*x^2+b*x*c)^2.
b = 2*m*(1-2*PoverS)
c = m**2
if (4*a*c-b**2>0) then
    part1 = (2/sqrt(4*a*c-b**2))*atan((2*a*ufinal+b)/sqrt(4*a*c-b**2))
    part1_0 = (2/sqrt(4*a*c-b**2))*atan((2*a*u0+b)/sqrt(4*a*c-b**2)) !Need this for calculating c2.
else if (4*a*c-b**2<0) then
    part1 = (1/sqrt(b**2-4*a*c))*log(abs((2*a*ufinal+b-sqrt(b**2-4*a*c))/(2*a*ufinal+b+sqrt(b**2-4*a*c))))
    part1_0 = (1/sqrt(b**2-4*a*c))*log(abs((2*a*u0+b-sqrt(b**2-4*a*c))/(2*a*u0+b+sqrt(b**2-4*a*c)))) !Need this for calculating c2.
else
    print*,'Error: P/S = 0 or 1'
end if
c2 = (4*m*c1)*((2*a*u0+b)/((4*a*c-b**2)*(a*u0**2+b*u0+c))+((2*a)/(4*a*c-b**2))*part1_0)
trishear_slip = -(4*m*c1)*((2*a*ufinal+b)/((4*a*c-b**2)*(a*ufinal**2+b*ufinal+c))+((2*a)/(4*a*c-b**2))*part1)+c2
return
end function trishear_slip

function trishear_u(u0,slip,m,PoverS,c1)
!use options, only: errs_id
!Calculates the u value that results from some given slip in the trishear zone.
use parameters, only: increment
implicit none
double precision :: u0,slip,m,PoverS,c1
double precision ::f,fprime,u,S
double precision :: trishear_u
integer :: i
real :: rn
f = 100.
CALL RANDOM_NUMBER(rn)
u = -m+2*m*rn !A random starting guess in the range [-m,m]
i = 1
!print*,'a'
do while (abs(f)>abs(increment))
    !write(errs_id,*),i
    S = trishear_slip(u0,u,m,PoverS,c1)
    f = S-slip
    fprime = -(4*m*c1)*((u**2+2*m*(1-2*PoverS)*u+m**2)**(-2))
    !fdoubleprime = (8*m*c1)*((u**2+2*m*(1-2*PoverS)*u+m**2)**(-3))*(2*u+2*m*(1-2*PoverS))
    u = u-(f/fprime)/(1-(f/fprime)*(-(2*u+2*m*(1-2*PoverS))/(u**2+2*m*(1-2*PoverS)*u+m**2))) !Halley's method. The last part is f''/(2*f') reduced.
    !u = u-f/fprime !Newton's method.
    if (u>m .or. u<-m) then
        !Try a new guess.
        CALL RANDOM_NUMBER(rn)
        u = -m+2*m*rn
    end if
    i = i+1
    if (i>1e6) then
        print*,'Too many iterations: exiting loop'
        exit
    end if
end do
trishear_u = u
return
end function trishear_u

end module analytic_trishear_module