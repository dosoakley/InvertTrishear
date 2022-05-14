%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ dipposze,dipsze ] = trishear_func_dip_decol( dipposze,dipsze,v0,m,total_slip,increment,PoverS,s,decolze,ramp_angle )
%Trishear velocity field
%   Calculates trishear velocity field and change of dip using the Zehnder and Allmendinger
%   (2000) solution and Oakley and Fisher (2015) change of dip equations 
%   and moves beds accordingly. It uses the eta-zeta
%   coordinate system with the origin at the initial fault tip, the zeta
%   axis parallel to the fault, and the eta axis perpendicular to it. It
%   also uses an x, y coordinate system that starts out the same as
%   eta,zeta but moves with the fault tip. This variant of the funciton
%   has a ramp coming up from a flat decollement. Trishear is only used around
%   the fault tip, while for the bend from flat to ramp, a simple kink is
%   used.

decol_ang = max(ramp_angle,pi-ramp_angle); %decollement angle from ramp
bisect_ang = (1/2)*decol_ang; %bisector angle
tan_bisect = tan(bisect_ang);
tan_ramp = tan(ramp_angle);
sin_ramp = sin(ramp_angle);
cos_ramp = cos(ramp_angle);
half_v0 = v0/2; %v0/2
v0term = half_v0/s; %v0/(2s) constant term that goes in front in dip equation
slip_sign = sign(v0); %+ if moving forward, - if moving back (restoring)
exp1 = (1+s)/(2*s); %exponent in 1st term for change in dip
exp2 = (1-s)/(2*s); %exponent in 2nd term for change in dip
xexp = 1/s; %exponent in vx term
yexp = (1+s)/s; %exponent in vy term
total_slip_abs = abs(total_slip); %magnitude of the total slip
sqrtm = sqrt(m); %square root of m
ndips = length(dipsze);
for j = 1:ndips
    x = dipposze(1,j);
    y = dipposze(2,j); 
    dip = dipsze(j);
    slip = 0; %slip so far        
    loc = 0; %location: ramp = 1, flat = 2, tri = 3, fw = 4            
    while abs(slip) < total_slip_abs
        rem_slip = total_slip - slip; %remaining slip
        xd = decolze(1) - PoverS*slip; %x for decollement start
        if loc == 0 %need to find out location
            if (y>0 && y>=m*x && y<=tan_bisect*(x-xd))    %ramp hanging wall
                loc = 1;
            elseif (y>=tan_ramp*(xd-x) && y>=m*x && y>0)  %flat hanging wall
                loc = 2;
            elseif (x<xd || y<=-m*x) %foot wall
                loc = 4;
            else %trishear zone
                loc = 3;
            end
        end
        if loc == 1    %ramp hanging wall
            rem_slip_x = rem_slip-rem_slip*PoverS; %rem_slip as change in x
            if slip_sign == -1 && y >= tan_bisect*(x+rem_slip-xd) %goes into flat
                hw_slip = y/tan_bisect-x+xd;
                loc = 2;
                %Solve for change in dip assuming slip conservation across syncline
                phi = ramp_angle; %don't really need 2 names for the same thing   
                theta = -dip; %Since going backwards, theta is dip with ramp as horizontal, but negative since phi is other way now
                gamma1 = pi/2-phi/2-theta;
                gamma2 = atan((tan(phi/2)-sin(theta)/(cos(phi/2)*cos(phi/2+theta)))^-1); % Precalculating at least tan(phi/2) and cos(phi/2) would help.
                if gamma2 < 0 %gamma2 should be in the range [0,pi]
                    gamma2 = pi+gamma2;
                end
                delta = pi-gamma1-gamma2; %change in dip
                dip = dip+delta;
            elseif y <= m*(x+rem_slip_x) %goes into trishear zone
                hw_slip = (y/m-x)/(1-PoverS);
                loc = 3;
            else %stays in ramp hw
                hw_slip = rem_slip;
            end
            x = x+hw_slip-hw_slip*PoverS;
            slip = slip + hw_slip;
        elseif loc == 2  %flat hanging wall
            if y>tan_bisect*(x+rem_slip*cos_ramp-xd)+rem_slip*sin_ramp %stays on flat
                fl_slip = rem_slip;
                x = x+fl_slip*cos_ramp-fl_slip*PoverS;
                y = y-fl_slip*sin_ramp;
            else %goes into hanging wall ramp
                fl_slip = slip_sign*(y-tan_bisect*(x-xd))/(sin_ramp+cos_ramp*tan_bisect);
                x = x+fl_slip*cos_ramp-fl_slip*PoverS;
                y = tan_bisect*(x-(xd-fl_slip*PoverS));
                loc = 1;
                %Solve for change in dip assuming slip conservation across syncline
                phi = ramp_angle; %don't really need 2 names for the same thing   
                theta = dip-ramp_angle; %Just rotate back to xy coord. system.
                gamma1 = pi/2-phi/2-theta;
                gamma2 = atan((tan(phi/2)-sin(theta)/(cos(phi/2)*cos(phi/2+theta)))^-1); % Precalculating at least tan(phi/2) and cos(phi/2) would help.
                if gamma2 < 0 %gamma2 should be in the range [0,pi]
                    gamma2 = pi+gamma2;
                end
                delta = pi-gamma1-gamma2; %change in dip
                dip = dip-delta;
            end
            slip = slip+fl_slip;
        elseif loc == 4 %foot wall               
            if slip_sign == 1 %forward motion, can't leave footwall
                fw_slip = rem_slip; %slip while in footwall
            else %inversion, can leave fw when trishear zone comes back to it
                if x-rem_slip*PoverS>0 && y > -m*(x-rem_slip*PoverS) %enters trishear zone                      
                    fw_slip = (x+y/m)/PoverS;
                    loc = 3;
                else
                    fw_slip = rem_slip;
                end
            end
            x = x-fw_slip*PoverS;
            slip = slip+fw_slip;
        else %trishear zone, loc = 3
            slip = slip+increment;
            dip = dip+(v0term/x)*(sqrtm*((abs(y)/(m*x))^(exp1))*cos(dip)+sign(y)*(1/sqrtm)*((abs(y)/(m*x))^(exp2))*sin(dip))^2;
            vx = half_v0*(sign(y)*(abs(y)/(m*x))^xexp+1)-PoverS*increment;
            vy = half_v0*(m/(1+s))*((abs(y)/(m*x))^yexp-1);
            x = x+vx;
            y = y+vy;
            loc = 0;
        end
    end
    dipposze(1,j) = x+PoverS*total_slip;
    dipposze(2,j) = y;
    dipsze(j) = dip;
end


end