%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ dippos_out,dips_out ] = parallel_FPF_func_dips( dippos,dips,v0,...
    total_slip,PoverS,bendxy,ramp_angle,gamma,gamma1,gamma_star,kappa,...
    beta1,beta2,tip_start)
%parallel_FPF_func_dip Parallel Fault Propagation Fold dips
%   Moves and rotates dips according to the rules of parallel fault
%   propagation folding, assuming initially flat beds.
bend = bendxy(:,1); %Only need this to be a 1D array.
ef_constant = sin(ramp_angle(1))*(sin(gamma_star-beta1)/(sin(gamma_star)*sin(2.*gamma_star-beta1))); %To be used in calculating the distance ef or the velocity of point e.
if (v0 >= 0) %determine the sign of the slip
    slip_sign = 1;
else
    slip_sign = -1;
end
%Calculate the slip ratios. All expresed relative to slip on the lowest segment.
R0 = sin(gamma1+ramp_angle(2))/sin(gamma1+ramp_angle(1));
R1 = R0*sin(gamma1+ramp_angle(1))/sin(gamma1+gamma);
R2 = R0*sin(beta2)/sin(beta2-ramp_angle(1)+gamma);
if (slip_sign == 1) %Forwards slip
    pte_start = tip_start; %Note: This would not work for normal faults, if they were allowed.
else %Backwards slip
    pte_start = tip_start-PoverS*R0*total_slip*ef_constant*[cos(kappa);sin(kappa)]; %minus because total_slip is negative.
end
%Calculate velocities for each domain
v1 = v0*[cos(ramp_angle(2));sin(ramp_angle(2))];
v2 = R0*v0*[cos(ramp_angle(1));sin(ramp_angle(1))];
v3 = R1*v0*[cos(gamma);sin(gamma)];
v4 = R2*v0*[cos(gamma);sin(gamma)];
v5 = [0;0];
%Calculate velocities for each boundary
vb = [0;0]; %Velocity of the bend: doesn't move.
ve = R0*v0*PoverS*[cos(ramp_angle(1))+ef_constant*cos(kappa);sin(ramp_angle(1))+ef_constant*sin(kappa)]; %Velocity of point e.
vt = R0*v0*PoverS*[cos(ramp_angle(1));sin(ramp_angle(1))]; %Velocity of the fault tip.
%Calculate slopes for each boundary
gam1_slope = -tan(gamma1); %Slope of axes with dip gamma, dipping to the right.
gam_slope = tan(gamma); %Slope of axes with dip gamma1, dipping to the left.
kap_slope = tan(kappa); %Slope of boundary with dip kappa.
for k = 1:size(dippos,2)
    tip = tip_start; %Where the tip starts moving from. (Generally the final tip position for restoring a fault-propagation fold.)
    pt = dippos(:,k); %The point to be moved.
    dip = dips(k); %The dip at the point to be moved.
    pte = pte_start; %Location e in Figure 11.6  of Allmendinger et al. (2012) and Fig. 16 of Suppe and Medwedeff (1990)
    slip = 0; %slip so far
    rem_slip = total_slip; %Remaining slip
    loc = 0; %location, zone numbers from Figure 11.6, Allmendinger et al. (2012), 5 = footwall, 0 = unknown
    loc_passed = zeros(1,5); %Tells if a point has already passed through one of the 5 zones.
    while (rem_slip*slip_sign>0. && isnan(pt(1))==0)
        %disp(loc)
        tip = tip + (slip/v0)*vt;
        pte = pte + (slip/v0)*ve;
        if (loc == 1) %Domain 1
            loc_passed(1) = 1;
            if (slip_sign == 1) %Forward motion
                slip2 = calc_slip(pt,bend,v0,v1,vb,gam1_slope); %Slip needed to reach domain 2
                if (loc_passed(2)==0 && rem_slip > slip2) %Goes into domain 2
                    slip = slip2;
                    dip = calc_dip_change(dip,v1,v2,vb,gam1_slope);
                    loc = 2;
                else
                    slip = rem_slip; %Stays in domain 1
                end
            else %Backwards motion
                slip = rem_slip; %Stays in domain 1
            end
            pt = pt+(slip/v0)*v1;
        elseif (loc == 2) %Domain 2
            loc_passed(2) = 1;
            if (slip_sign == 1) %Forward motion
                slip3 = calc_slip(pt,pte,v0,v2,ve,gam1_slope); %Slip needed to reach domain 3
                if (loc_passed(3)==0 && rem_slip > slip3 && pt(2)+(slip3/v0)*v2(2) > pte(2)+(slip3/v0)*ve(2) && slip3>0) %Goes into domain 3.  I changed 0 to 1e-13 to address an infinite loop bug probably caused by rounding.
                    slip = slip3;
                    dip = calc_dip_change(dip,v2,v3,ve,gam1_slope);
                    loc = 3;
                else
                    slip4 = calc_slip(pt,tip,v0,v2,vt,kap_slope); %Slip needed to reach domain 4
                    if (loc_passed(4)==0 && rem_slip>slip4 && pt(2)+(slip4/v0)*v2(2) < pte(2)+(slip4/v0)*ve(2) && slip4>0) %Goes into domain 4
                        slip = slip4;
                        dip = calc_dip_change(dip,v2,v4,vt,kap_slope);
                        loc = 4;
                    else
                        slip = rem_slip; %Stays in domain 2
                    end
                end
            else %Backward motion
                slip1 = calc_slip(pt,bend,v0,v2,vb,gam1_slope);
                slip3 = calc_slip(pt,pte,v0,v2,ve,gam1_slope);
                slip4 = calc_slip(pt,tip,v0,v2,vt,kap_slope);
                if (loc_passed(4)==0 && rem_slip < slip4 && slip4<0 && pt(2)+(slip4/v0)*v2(2) < pte(2)+(slip4/v0)*ve(2) && slip1<slip4) %Goes into domain 4.
                    slip = slip4;
                    dip = calc_dip_change(dip,v2,v4,vt,kap_slope);
                    loc = 4;
                elseif (loc_passed(3)==0 && rem_slip < slip3 && slip3<0 && pt(2)+(slip3/v0)*v2(2) > pte(2)+(slip3/v0)*ve(2) ...
                        && slip1<slip3) %Goes into domain 3.
                    slip = slip3;
                    dip = calc_dip_change(dip,v2,v3,ve,gam1_slope);
                    loc = 3;
                elseif (loc_passed(1)==0 && rem_slip < slip1) %Goes into domain 1. < because slips are negative
                    slip = slip1;
                    dip = calc_dip_change(dip,v2,v1,vb,gam1_slope);
                    loc = 1;
                else
                    slip = rem_slip; %Stays in domain 2
                end
            end
            pt = pt+(slip/v0)*v2;
        elseif (loc == 3) %Domain 3
            loc_passed(3) = 1;
            if (slip_sign == 1) %Forward motion
                slip2 = calc_slip(pt,pte,v0,v3,ve,gam1_slope); %Slip needed to reach domain 2
                slip4 = calc_slip(pt,pte,v0,v3,ve,gam_slope); %Slip needed to reach domain 4.
                if (loc_passed(2)==0 && rem_slip > slip2 && (slip4>slip2 || slip4<0 || loc_passed(4)==1) && slip2 > 0) %Goes into domain 2. I changed 0 to 1e-13 to address an infinite loop bug probably caused by rounding.
                    slip = slip2;
                    dip = calc_dip_change(dip,v3,v2,ve,gam1_slope);
                    loc = 2;
                elseif (loc_passed(4)==0 && rem_slip > slip4 && slip4>0) %Goes into domain 4
                    slip = slip4;
                    loc = 4;
                    dip = calc_dip_change(dip,v3,v4,ve,gam_slope);
                else
                    slip = rem_slip; %Stays in domain 3
                end
                %Slip in domain 3 is parallel to the domain 3-4 boundary, so the point cannot enter domain 4 if slip (and thus propagation) is forward.
            else %Backward motion
                slip2 = calc_slip(pt,pte,v0,v3,ve,gam1_slope); %Slip needed to reach domain 2.
                slip4 = calc_slip(pt,pte,v0,v3,ve,gam_slope); %Slip needed to reach domain 4.
                if (loc_passed(2)==0 && rem_slip < slip2 && (slip4<slip2 || slip4>0 || loc_passed(4)==1) && slip2<0) %Goes into domain 2. I'm pretty sure this can only happen if P/S <1, otherwise material flows the other way, but I don't think it's quite that simple.
                    slip = slip2;
                    dip = calc_dip_change(dip,v3,v2,ve,gam1_slope);
                    loc = 2;
                elseif (loc_passed(4)==0 && rem_slip < slip4 && slip4<0) %Goes into domain 4. Pte must move along a line less steep than gamma.
                    slip = slip4;
                    dip = calc_dip_change(dip,v3,v4,ve,gam_slope);
                    loc = 4;
                else
                    slip = rem_slip;
                end
            end
            pt = pt+(slip/v0)*v3;
        elseif (loc == 4) %Domain 4
            loc_passed(4) = 1;
            if (slip_sign == 1); %Forward motion
                slip3 = calc_slip(pt,pte,v0,v4,ve,gam_slope); %Slip needed to reach domain 3.
                if (loc_passed(3)==0 && rem_slip>slip3 && pt(2)+(slip3/v0)*v4(2) > pte(2)+(slip3/v0)*ve(2) && slip3 > 0)
                    slip = slip3;
                    dip = calc_dip_change(dip,v4,v3,ve,gam_slope);
                    loc = 3;
                else
                    slip2 = calc_slip(pt,tip,v0,v4,vt,kap_slope); %Slip needed to reach domain 2.
                    if (loc_passed(2)==0 && rem_slip > slip2 && pt(2)+(slip2/v0)*v4(2) < pte(2)+(slip2/v0)*ve(2) && slip2 > 0)
                        slip = slip2;
                        dip = calc_dip_change(dip,v4,v2,vt,kap_slope);
                        loc = 2;
                    else
                        slip = rem_slip;
                    end
                end
            else %Backward motion
                slip2 = calc_slip(pt,tip,v0,v4,vt,kap_slope); %Slip needed to reach domain 2.
                slip3 = calc_slip(pt,pte,v0,v4,ve,gam_slope); %Slip needed to reach domain 3.
                slip5 = calc_slip(pt,tip,v0,v4,vt,gam_slope); %Slip needed to reach domain 5 (footwall).
                if (loc_passed(2)==0 && rem_slip < slip2 && (slip5<slip2 || slip5>0 || loc_passed(5)==1) ...
                        && pt(2)+(slip2/v0)*v4(2) < pte(2)+(slip2/v0)*ve(2) && slip2<0)
                    slip = slip2;
                    dip = calc_dip_change(dip,v4,v2,vt,kap_slope);
                    loc = 2;
                elseif (loc_passed(3)==0 && rem_slip < slip3 &&(slip5<slip3 || slip5>0 || loc_passed(5)==1) ...
                        && pt(2)+(slip3/v0)*v4(2) > pte(2)+(slip3/v0)*ve(2) && slip3<0)
                    slip = slip3;
                    dip = calc_dip_change(dip,v4,v3,ve,gam_slope);
                    loc = 3;
                elseif (loc_passed(5)==0 && rem_slip < slip5 && slip5<0)
                    slip = slip5;
                    dip = calc_dip_change(dip,v4,v5,vt,gam_slope);
                    loc = 5;
                else
                    slip = rem_slip;
                end
            end
            pt = pt+(slip/v0)*v4;
        elseif (loc == 5) %Footwall (domain 5)
            loc_passed(5) = 1;
            if (slip_sign == 1) %Forward motion
                slip4 = calc_slip(pt,tip,v0,v5,vt,gam_slope); %Slip needed to reach domain 4.
                if (loc_passed(4)==0 && rem_slip > slip4 && pt(2)+(slip4/v0)*v5(2) > tip(2)+(slip4/v0)*vt(2))
                    slip = slip4;
                    dip = calc_dip_change(dip,v5,v4,vt,gam_slope);
                    loc = 4;
                else
                    slip = rem_slip;
                end
            else
                slip = rem_slip;
            end
        else %Need to find out where the point is.
            slip = 0; %So no more slip gets subtracted from rem_slip while we're just deciding.
            if (loc_passed(5)==0 && (pt(1)<bend(1) && pt(2)-bend(2)<(pt(1)-bend(1))*tan(ramp_angle(2))) || ...
                    (pt(1)>bend(1) && pt(1)<tip(1) && pt(2)-bend(2)<(pt(1)-bend(1))*tan(ramp_angle(1))) || ...
                    (pt(1)>tip(1) && pt(2)-tip(2)<(pt(1)-tip(1))*gam_slope)) %Footwall
                loc = 5;
            elseif (loc_passed(1)==0 && pt(2)-bend(2)<(pt(1)-bend(1))*gam1_slope) %Domain 1
                loc = 1;
            elseif (loc_passed(2)==0 && (pt(2)<pte(2) && pt(1)-tip(1)<(pt(2)-tip(2))/kap_slope) || ...
                    (pt(2)>pte(2) && pt(1)-pte(1)<(pt(2)-pte(2))/gam1_slope)) %Domain 2
                loc = 2;
            elseif (loc_passed(3)==0 && pt(2)>pte(2) && pt(1)-pte(1)<(pt(2)-pte(2))/gam_slope) %Domain 3
                loc = 3;
            elseif (loc_passed(4)==0 && pt(1)-tip(1)<(pt(2)-tip(2))/gam_slope) %Domain 4
                loc = 4;
            else
                disp('Error: point is not in any of the domains') %This shouldn't happen.
            end
        end
        rem_slip = rem_slip - slip; %remaining slip
    end
    dippos(:,k) = pt; %Store the moved point back in the pts array.
    dips(k) = dip; %Store the changed dip.
end
dippos_out = dippos;
dips_out = dips;

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