%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ pts,dips ] = trishear_func_dip_multi_bend( pts,dips,tip_start,tip_seg_start,...
    tip_seg_end,v0,phi,m,increment,PoverS,s,bendxy,ramp_angle,ramp_dir,...
    axis_dips,R,slipseg )
%trishear_func_multi_bend Fault with a series of bends.
%   This function moves points according to the kinematics of trishear
%   fault propagation folding for a fault consisting of multiple segments,
%   with changes in ramp angle, phi, or P/S between them. Segments are
%   numbered from the top (segment number 1) down.
%   pts = the points to be moved
%   dips = the dips at the locations in pts.
%   tip_start = [x,y] position to start the fault tip at
%   tip_seg_start = segment number that the fault tip starts in
%   tip_seg_end = segment number that the fault tip ends in
%   v0 = slip velocity = increment
%   phi = trishear apical angle, vector of length nsegments
%   m = tan(phi)
%   increment = increment of slip = v0
%   PoverS = propagation to slip ratio, vector of length nsegments
%   s = center concentration factor for trishear
%   bendxy = [x,y] positions of all fault bends. 2 x nsegments vector
%   ramp_angle = fault dip, vector of length nsegments
%   ramp_dir = direction of fault dip
%   axis_dips = dips of fold axes
%   R = ratio of slip between segments
%   slipseg = amount of slip in each fault segment.

nbends = size(bendxy,2);
nsegs = nbends+1;
nangles = 1;
for i = 2:nsegs %Count how many different ramp angles there are.
    if ramp_angle(i) ~= ramp_angle(i-1)
        nangles = nangles+1;
    end
end

tip_loc = zeros(1,nsegs);
ramp_angle_loc = zeros(1,nangles);
R_loc = zeros(1,nangles);
bendxy_loc = zeros(2,nangles-1);
axis_slopes = zeros(1,nangles-1);
i = 1;
ramp_angle_loc(1) = ramp_angle(1);
R_loc(1) = R(1);
tip_loc(1) = 1;
for n = 2:nsegs
    if (ramp_angle(n) ~= ramp_angle(n-1))
        ramp_angle_loc(i+1) = ramp_angle(n);
        R_loc(i+1) = R(n);
        bendxy_loc(:,i) = bendxy(:,n-1);
        axis_slopes(i) = tan(axis_dips(n-1));
        i = i+1;
    end
    tip_loc(n) = i; %This lets us translate a tip segment to the corresponding value of the "loc" variable.
end
sin_ramp_loc = sin(ramp_angle_loc);
cos_ramp_loc = cos(ramp_angle_loc);
tan_ramp_loc = tan(ramp_angle_loc);

sin_ramp = sin(ramp_angle);
cos_ramp = cos(ramp_angle);
half_v0 = v0/2; %v0/2
xexp = 1/s; %exponent in vx term
yexp = (1+s)/s; %exponent in vy term
v0term = half_v0/s; %v0/(2s) constant term that goes in front in dip equation
exp1 = (1+s)/(2*s); %exponent in 1st term for change in dip
exp2 = (1-s)/(2*s); %exponent in 2nd term for change in dip
sqrtm = sqrt(m); %square root of m
if (v0 >= 0) %determine the sign of the slip
    slip_sign = 1;
else
    slip_sign = -1;
end
nseg_tip = abs(tip_seg_end-tip_seg_start)+1; %Number of segments tip will pass through.
nbends = size(bendxy,2);
nsegs = nbends+1;
for n = tip_seg_start:-slip_sign:tip_seg_end %Loop through each segment that the tip goes through.
    for k = 1:length(dips)
        pt = pts(:,k);
        if (n == tip_seg_start)
            tip = tip_start;
        else
            tip = bendxy(:,n+(slip_sign-1)/2);
        end
        slip = 0; %slip so far
        loc = 0; %location: 0 = unknown, 1 to nsegs = fault segments, nsegs+1 = trishear, nsegs+2 = footwall
        dip = dips(k);
        rem_slip = slipseg(n);
        while (rem_slip*slip_sign>0)
            if (loc == 0) %need to find out location
                %I also need to deal with what happens if two fold axes intersect.
                if ((ramp_dir*(pt(1)-tip(1))>ramp_dir*(pt(2)-tip(2))/tan(ramp_angle(n)+ramp_dir*phi(n))) &&...
                        ((pt(2)-tip(2))>(pt(1)-tip(1))*tan(ramp_angle(n)-ramp_dir*phi(n)))) %Trishear zone
                    loc = nangles+1;
                else
                    for i = tip_loc(n):nangles-1 %Segment i is always the segment above bend i.
                        if (pt(1)*ramp_dir>bendxy_loc(1,i)*ramp_dir && pt(2)-bendxy_loc(2,i)<(pt(1)-bendxy_loc(1,i))*tan_ramp_loc(i))
                            loc = nangles+2; %Footwall, below segment i
                            break
                        elseif (ramp_dir*(pt(1)-bendxy_loc(1,i))>ramp_dir*(pt(2)-bendxy_loc(2,i))/axis_slopes(i)) %In the hangingwall
                            loc = i; %Above segment i
                            break
                        end
                    end
                    if (loc ==0) %If we still haven't assigned it to anywhere, the only remaining place is above or below the lowest ramp.
                        if (pt(2)-bendxy_loc(2,nangles-1)>(pt(1)-bendxy_loc(1,nangles-1))*tan_ramp_loc(nangles))
                            loc = nangles; %Hanging wall, lowest segment.
                        else
                            loc = nangles+2; %Footwall, lowest segment.
                        end
                    end
                end
            elseif (loc == nangles+1) %Trishear
                ptze(:,1) = pt;
                ptze = xy_to_ze2(ptze,tip(1),tip(2),ramp_angle(n)); %Rotate into trishear coordinates.
                dip = (ramp_dir*dip+min(ramp_angle(n),pi-ramp_angle(n))); %Convert dips to ze coord. system
                x = ptze(1,1);
                y = ptze(2,1);
                slip = 0;
                while (slip == 0 || (y<=m(n)*x && y>=-m(n)*x))
                    if (y >= 0)
                        signy = 1;
                    else
                        signy = -1;
                    end
                    %First slip
                    slip = slip+increment;
                    dip = dip+(R(n)*v0term/x)*(signy*sqrtm(n)*((abs(y)/(m(n)*x))^(exp1))*cos(dip)+(1/sqrtm(n))...
                        *((abs(y)/(m(n)*x))^(exp2))*sin(dip))^2;
                    vx = R(n)*half_v0*((signy*(abs(y)/(m(n)*x))^xexp)+1);
                    vy = R(n)*half_v0*(m(n)/(1+s))*(((abs(y)/(m(n)*x))^yexp)-1);
                    x = x+vx-R(n)*PoverS(n)*increment;
                    y = y+vy;
                    if (abs(slip)>abs(rem_slip))
                        break
                    end
                end
                ptze(:,1) = [x+R(n)*PoverS(n)*slip,y]; %The +R(n)*PoverS(n)*slip is because we are transforming back using the same tip position we used to transform into trishear coordinates.
                ptze = ze_to_xy2(ptze,tip(1),tip(2),ramp_angle(n));%Rotate out of trishear coordinates.
                pt = ptze(:,1);
                dip = ramp_dir*(dip - min(ramp_angle(n),pi-ramp_angle(n))); %Convert dip back into section coordinates.
                loc = 0;
            elseif(loc==nangles+2) %In the footwall
                vb = v0*R(n)*PoverS(n)*[cos_ramp(n);sin_ramp(n)]; %Velocity of trishear boundary.
                slip_tri = calc_slip(pt,tip,v0,[0.0d0;0.0d0],vb,tan(ramp_angle(n)-ramp_dir*phi(n)));
                if((n==nsegs || pt(1)*ramp_dir>bendxy(1,n)*ramp_dir) && slip_sign*slip_tri<slip_sign*rem_slip  && slip_sign*slip_tri>0) %Goes into trishear segment.
                    slip = slip_tri;
                    loc = nangles+1;
                else %Stays in FW.
                    slip = rem_slip;
                end
            elseif (loc==tip_loc(n)) %In the tip segment.
                v = v0*R(n)*[cos_ramp(n);sin_ramp(n)];
                vb = v0*R(n)*PoverS(n)*[cos_ramp(n);sin_ramp(n)];
                slip_tri = calc_slip(pt,tip,v0,v,vb,tan(ramp_angle(n)+ramp_dir*phi(n)));
                if (slip_sign*slip_tri<slip_sign*rem_slip && slip_sign*slip_tri>0)
                    slip = slip_tri;
                    loc = nangles+1;
                elseif(slip_sign==-1 && loc<nangles)
                    slip_next = calc_slip(pt,bendxy_loc(:,loc),v0,v,[0.0;0.0],axis_slopes(loc));
                    if (slip_next>rem_slip) %Greater than because slip_sign = -1
                        slip = slip_next;
                        vnext = v0*R_loc(loc+1)*[cos_ramp_loc(loc+1);sin_ramp_loc(loc+1)];
                        dip = calc_dip_change(dip,v,vnext,[0,0],axis_slopes(loc));
                        loc = loc+1;
                    else
                        slip = rem_slip;
                    end
                else
                    slip = rem_slip;
                end
                pt = pt+v*(slip/v0);
            else %In one of the fault segments
                v = v0*R_loc(loc)*[cos_ramp_loc(loc);sin_ramp_loc(loc)];
                if (loc == nangles && slip_sign == -1) %In the last segment and moving backwards. Can't go anywhere else.
                    slip = rem_slip;
                else
                    if(slip_sign==1) %Forwards slip, can slip into segment above.
                        slip_next = calc_slip(pt,bendxy_loc(:,loc-1),v0,v,[0.0,0.0],axis_slopes(loc-1));
                    else    %Backwards slip, can slip into segment below.
                        slip_next = calc_slip(pt,bendxy_loc(:,loc),v0,v,[0.0,0.0],axis_slopes(loc));
                    end
                    if(slip_next*slip_sign<rem_slip*slip_sign && slip_sign*slip_next>0)
                        slip = slip_next;
                        vnext = v0*R_loc(loc-slip_sign)*[cos_ramp_loc(loc-slip_sign),sin_ramp_loc(loc-slip_sign)];
                        dip = calc_dip_change(dip,v,vnext,[0,0],axis_slopes(loc-(1+slip_sign)/2));
                        loc = loc-slip_sign;
                    else
                        slip = rem_slip;
                    end
                end
                pt = pt+v*(slip/v0);
            end
            rem_slip = rem_slip - slip; %remaining slip
            tip(1) = tip(1)+slip*R(n)*PoverS(n)*cos_ramp(n);
            tip(2) = tip(2)+slip*R(n)*PoverS(n)*sin_ramp(n);
            slip = 0; %Otherwise, when we go to loc = 0, we reuse the same slip twice.
        end
        pts(:,k) = pt;
        dips(k) = dip;
    end
end

function slip = calc_slip(pt,ptb,v0,v1,vb,slope)
    %This function calculates the slip necessary for a point moving in a constant velocity domain to reach a boundary.
    %Everything with dimension(2) is x coordinate or x component followed by y coordinate or y component.
    slip = ((pt(1)-ptb(1))*slope-pt(2)+ptb(2))/((1./v0)*(v1(2)-vb(2)-(v1(1)-vb(1))*slope));
end

function newdip = calc_dip_change(dip,v1,v2,vb,slope)
    %This function calculates the changed dip for a point crossing a boundary between two constant velocity domains.
    %Everything with dimension(2) is x coordinate or x component followed by y coordinate or y component.
    dip = -dip; %An unfortunate consequence of the sign convention I have been using, with dips positive down to the right.
    tan_dip = tan(dip);
    % secsq_dip = (sec(dip))^2;
    A = v1(1)*tan_dip-v1(2);
    B = vb(1)*slope-vb(2);
    C = slope-tan_dip;
    W = A*slope-B*tan_dip+C*v2(2);
    X = A-B+C*v2(1);
    % Y = v1(1)*secsq_dip*slope-B*secsq_dip-v2(2)*secsq_dip;
    % Z = v1(1)*secsq_dip-v2(1)*secsq_dip;
    dip = atan(W/X); %New dip
    % sigma_dip = abs((Y*X-W*Z)/(W**2+X**2))*sigma_dip !New uncertainty in dip.
    newdip = -dip; %Switch back to the sign convention used in the rest of the program.
end

end

