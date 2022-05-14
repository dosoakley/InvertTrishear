%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ bedsze ] = trishear_func( bedsze,v0,m,total_slip,increment,PoverS,s )
%Trishear velocity field
%   Calculates trishear velocity field using the Zehnder and Allmendinger
%   (2000) solution and moves beds accordingly. It uses the eta-zeta
%   coordinate system with the origin at the initial fault tip, the zeta
%   axis parallel to the fault, and the eta axis perpendicular to it. It
%   also uses an x, t coordinate system that starts out the same as
%   eta,zeta but moves with the fault tip

half_v0 = v0/2; %v0/2
xexp = 1/s; %exponent in vx term
yexp = (1+s)/s; %exponent in vy term
total_slip_abs = abs(total_slip); %magnitude of the total slip
slip_sign = sign(v0); %+ if moving forward, - if moving back (restoring)
prop_sign = sign(PoverS);
for j = 1:size(bedsze,3)
    for k = 1:size(bedsze,2)
        x = bedsze(1,k,j);
        y = bedsze(2,k,j);
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
                        if (fw_slip == 0) %This is to avoid the case where x = 0 and y = 0 which first goes here then to trishear, where 0/0 causes a NaN
                            fw_slip = increment;
                        end
                    else
                        fw_slip = rem_slip;
                    end
                end
                x = x-fw_slip*PoverS;
                slip = slip+fw_slip;
                fw = 1;
            else %trishear zone
                slip = slip+increment;
                vx = half_v0*(sign(y)*(abs(y)/(m*x))^xexp+1)-PoverS*increment;
                vy = half_v0*(m/(1+s))*((abs(y)/(m*x))^yexp-1);
                x = x+vx;
                y = y+vy;
            end
        end
        bedsze(1,k,j) = x+PoverS*total_slip;
        bedsze(2,k,j) = y;
    end
end

end

