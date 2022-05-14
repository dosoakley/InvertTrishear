%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ptsze] = trishear_func_bend(ptsze,v0,m,total_slip,increment,PoverS,s,bendze,ramp_angle)
ramp_diff = ramp_angle(1)-ramp_angle(2); %Difference between the two ramp angles.
bend_ang = max(abs(ramp_diff),pi-abs(ramp_diff)); %full angle of the bend.
bisect_ang = 0.5*bend_ang; %bisector angle
tan_bisect = tan(bisect_ang); %Just calculate once = faster
tan_diff = tan(ramp_diff);
sin_diff = sin(ramp_diff);
cos_diff = cos(ramp_diff);
half_v0 = v0/2; %v0/2
xexp = 1/s; %exponent in vx term
yexp = (1+s)/s; %exponent in vy term
total_slip_abs = abs(total_slip);
bendz = bendze(1,1); %zeta position of the bend
if (v0 >= 0) %determine the sign of the slip
    slip_sign = 1;
else
    slip_sign = -1;
end
if (bendz <= 0) %Determine which segment the tip  is in to start.
    tip_seg = 1;
else
    tip_seg = 2;
end
if ((tip_seg == 1 && slip_sign == 1) || (tip_seg == 2 && slip_sign == -1)) %Find out how many segments the tip passes through
    nseg = 1;
else
    dist = sqrt(bendze(1,1)^2 + bendze(2,1)^2); %Note bendze(2,1) should be 0
    if (dist > total_slip_abs*PoverS)
        nseg = 1;
    else
        nseg = 2;
    end
end
if (nseg == 1)
    slipseg = total_slip_abs;
else
    slipseg = [dist/PoverS,total_slip_abs-dist/PoverS];
end
for n = 1:nseg
    if (n==2) %On the second segment. Need to transform points.
        %Note: I probably need to rotate the uncertainties too.
        tip_seg = tip_seg-slip_sign; %Should change 2 to 1 or vice versa
        ptsze = rot2D(ptsze,ramp_angle(tip_seg+slip_sign)-ramp_angle(tip_seg)); %The angle should be negative if going from 2 to 1, positive if 1 to 2; 1>2
        bendz = 0; %At this point of change, the fault tip is at the bend
    end
    for j = 1:size(ptsze,3)
        for k = 1:size(ptsze,2)
            x = ptsze(1,k,j);
            y = ptsze(2,k,j);
            slip = 0; %slip so far
            loc = 0; %location: ramp = 1, flat = 2, tri = 3, fw = 4
            while (abs(slip)<slipseg(n))
                rem_slip = slip_sign*slipseg(n) - slip; %remaining slip
                xd = bendz-PoverS*slip; %x for bend start
                if (loc == 0) %need to find out location
                    if (y>0 && y>=m*x && (tip_seg==2 || y<=tan_bisect*(x-xd))) %upper ramp hanging wall (or lower ramp if only one segment)
                        loc = 1;
                    elseif (tip_seg==1 && y>=tan_diff*(xd-x) && y>=m*x && y>0)  %lower hanging wall
                        loc = 2;
                    elseif ((tip_seg==1 && x<xd) || y<=-m*x) %foot wall
                        loc = 4;
                    else %trishear zone
                        loc = 3;
                    end
                end
                if (loc == 1) %ramp hanging wall
                    if (slip_sign == 1) %Forward motion
                        if (y <= m*(x+rem_slip*(1-PoverS))) %goes into trishear zone
                            hw_slip = (y/m-x)/(1-PoverS);
                            loc = 3;
                        else %stays in ramp hw
                            hw_slip = rem_slip;
                        end
                    else %Backward motion / thrust inversion
                        slip2fl = y/tan_bisect-x+xd; %Slip required to reach the hanging wall flat
                        slip2tri = (y/m-x)/(1-PoverS); %Slip required to reach the trishear zone
                        if (PoverS > 1. && (slip2tri > slip2fl || tip_seg==2) && slip2tri > rem_slip) %Goes into trishear zone, > b/c all slips are negative (but not if slip_sign=1).
                            hw_slip = slip2tri;
                            loc = 3; % Go to trishear zone
                        elseif (tip_seg==1 && slip2fl > rem_slip) %Goes into flat
                            hw_slip = slip2fl;
                            loc = 2; %Go to flat
                        else %Stays in ramp hanging wall
                            hw_slip = rem_slip;
                        end
                    end
                    x = x+hw_slip-hw_slip*PoverS;
                    slip = slip + hw_slip;
                elseif (loc == 2) %flat hanging wall
                    if (y>tan_bisect*(x+rem_slip*cos_diff-xd)+rem_slip*sin_diff) %stays on flat
                        fl_slip = rem_slip;
                        x = x+fl_slip*cos_diff-fl_slip*PoverS;
                        y = y-fl_slip*sin_diff;
                    else %goes into hanging wall ramp
                        fl_slip = slip_sign*(y-tan_bisect*(x-xd))/(sin_diff+cos_diff*tan_bisect);
                        x = x+fl_slip*cos_diff-fl_slip*PoverS;
                        y = tan_bisect*(x-(xd-fl_slip*PoverS));
                        loc = 1;
                    end
                    slip = slip+fl_slip;
                elseif (loc == 4) %foot wall
                    if (slip_sign == 1) %forward motion, can't leave footwall
                        fw_slip = rem_slip; %slip while in footwall
                    else %inversion, can leave fw when trishear zone comes back to it
                        if (x-rem_slip*PoverS>0 && y > -m*(x-rem_slip*PoverS)) %enters trishear zone
                            fw_slip = (x+y/m)/PoverS;
                            loc = 3;
                        else
                            fw_slip = rem_slip;
                        end
                    end
                    x = x-fw_slip*PoverS;
                    slip = slip+fw_slip;
                else %trishear zone, loc = 3
                    if (y >= 0)
                        signy = 1;
                    else
                        signy = -1;
                    end
                    slip = slip+increment;
                    vx = half_v0*((signy*(abs(y)/(m*x))^xexp)+1)-PoverS*increment;
                    vy = half_v0*(m/(1+s))*(((abs(y)/(m*x))^yexp)-1);
                    x = x+vx;
                    y = y+vy;
                    loc = 0;
                end
            end
            ptsze(1,k,j) = x; %This is not really transforming back to zeta, eta but keeping with the new tip
            ptsze(2,k,j) = y;
        end
    end
end
end

function pts_rot = rot2D(pts,angle) %Rotates n points in a 2xnxm array through a given angle
angle_cos = cos(angle);
angle_sin = sin(angle);
pts_rot=zeros(size(pts,1),size(pts,2),size(pts,3));
for j = 1:size(pts,3)
    for i = 1:size(pts,2)
        pts_rot(1,i,j) = angle_cos*pts(1,i,j)-angle_sin*pts(2,i,j);
        pts_rot(2,i,j) = angle_sin*pts(1,i,j)+angle_cos*pts(2,i,j);
    end
end
end