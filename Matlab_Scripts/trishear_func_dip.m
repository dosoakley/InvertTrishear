%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ dipposze,dipsze ] = trishear_func_dip( dipposze,dipsze,v0,m,total_slip,increment,PoverS,s )
%Trishear_func_dip
%Calculates the change of dip of a point as it moves through the Trishear
%zone. Based on the velocity equations of Zehnder and Allmendinger (2000)
%and the equations for change in dip from Oakley and Fisher (2015).
%dipposze should be a 2D matrix, not 3D.

half_v0 = v0/2; %v0/2
v0term = half_v0/s; %v0/(2s) constant term that goes in front in dip equation
exp1 = (1+s)/(2*s); %exponent in 1st term for change in dip
exp2 = (1-s)/(2*s); %exponent in 2nd term for change in dip
xexp = 1/s; %exponent in vx term
yexp = (1+s)/s; %exponent in vy term
total_slip_abs = abs(total_slip); %magnitude of the total slip
slip_sign = sign(v0); %+ if moving forward, - if moving back (restoring)
prop_sign = sign(PoverS);
sqrtm = sqrt(m); %square root of m
ndips = length(dipsze);
for j = 1:ndips
    x = dipposze(1,j);
    y = dipposze(2,j);  
    dip = dipsze(j);
    slip = 0;
    hw = 0; %tells whether the point has been in the hanging wall yet
    fw = 0; %tells whether the point has been in the footwall yet        
    while abs(slip) < total_slip_abs
        rem_slip = total_slip - slip; %remaining slip
        if hw==0 && y>0 && y>=m*x %above trishear zone (hanging wall)
            rem_slip_x = rem_slip-rem_slip*PoverS; %rem_slip as change in x
            if y >= m*(x+rem_slip_x)
                hw_slip = rem_slip; %slip while in the hw
            else
                hw_slip = (y/m-x)/(1-PoverS);
            end
            x = x+hw_slip-hw_slip*PoverS;
            slip = slip + hw_slip;
            hw = 1;
        elseif fw==0 && y<=-m*x %below trishear zone (foot wall)
            if slip_sign == prop_sign %forward motion, can't leave footwall
                fw_slip = rem_slip; %slip while in footwall
            else %inversion, can leave fw when trishear zone comes back to it
                if y > -m*(x-rem_slip*PoverS) %enters trishear zone                      
                    fw_slip = (x+y/m)/PoverS;
                else
                    fw_slip = rem_slip;
                end
            end
            x = x-fw_slip*PoverS;
            slip = slip+fw_slip;
            fw = 1;
        else %trishear zone
            slip = slip+increment;
            dip = dip+(v0term/x)*(sqrtm*((abs(y)/(m*x))^(exp1))*cos(dip)+sign(y)*(1/sqrtm)*((abs(y)/(m*x))^(exp2))*sin(dip))^2;
            vx = half_v0*(sign(y)*(abs(y)/(m*x))^xexp+1)-PoverS*increment;
            vy = half_v0*(m/(1+s))*((abs(y)/(m*x))^yexp-1);
            x = x+vx;
            y = y+vy;
        end
    end
    dipposze(1,j) = x+PoverS*total_slip;
    dipposze(2,j) = y;
    dipsze(j) = dip;
end

end