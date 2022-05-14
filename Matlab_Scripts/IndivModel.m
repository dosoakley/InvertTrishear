%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

%Analyze an individual model.
%You will need to create a parameter file containing just the parameters
%for this one model. Only a single column is necessary.

%This only works for the first three fault types. It also was developed
%with an older version of the program and may require modification to work
%with the latest version.

%Read in the data
NDatTypes = InputGoodOnly('Enter number of data types: ',[1,2,3]);
ordinals = {'1st','2nd','3rd'};
DatTypes = zeros(1,3); %Tells whether each data type is present or not.
for n = 1:NDatTypes
    DatType = InputGoodOnly(['Enter ',ordinals{n},' data type:\n(1) beds\n(2) dips\n(3)fault points\n'],[1,2,3]);
    switch DatType
        case 1 %Beds
            DatTypes(1) = 1;
            %open a file with bed geometry
            file = -1;
            while file == -1
                file_name = input('input bedding file name (with extension): ','s');
                file = fopen(file_name);
                if file == -1
                    disp('No such file.')
                end
            end
            beds = ReadBeds(file,1);
            nbeds = size(beds,3);
            fclose(file);
            disp(['Read ',num2str(nbeds),' beds'])
        case 2 %Dips
            DatTypes(2) = 1;
            %open a file with dips and their positions
            file_name = input('Input dip file name (with extension): ','s');
            dip_profile = load(file_name);
            if size(dip_profile,2)>size(dip_profile,1) && size(dip_profile,2)>4 %In case it was formatted other way
                dip_profile = dip_profile';
            end
            dip_pos = [dip_profile(:,1)';dip_profile(:,2)'];
            dips = dip_profile(:,3)';
            rad_dips = dips*pi/180; %dips in radians
        case 3 %Fault Points
            DatTypes(3) = 1;
            %open a file containing fault points
            file_name = input('Input fault points file name (with extension): ','s');
            fault_pts = load(file_name);
            fault_pts = fault_pts';
    end
end

%Find out the type of fault
FaultType = InputGoodOnly(['Enter Fault Type: \n(1) Straight Fault\n(2)'...
    ' Ramp from Horizontal Detachment\n(3) Ramp with a bend in it\n'],[1,2,3]);

%Import fault parameters
file = -1;
while file == -1
    file_name = input('Input individual model parameter file name (with extension): ','s');
    file = fopen(file_name);
    if file == -1
        disp('No such file.')
    end
end
params = textscan(file,'%f %f %f');
fclose(file);
nextra = 0; %Keeps track of extra parameters that are read only in some cases.
tipx = params{1}(1); %x coordinate of fault tip
tipy = params{1}(2); %y coordinate of fault tip
total_slip = params{1}(3); %total slip on the fault
ramp_angle = params{1}(4); %angle of fault ramp to horizontal
ramp_angle = ramp_angle*pi/180; %convert to radians
ramp_dir = sign(tan(ramp_angle)); %1=dips left, -1 = dips right
ramp_acute = min(ramp_angle,pi-ramp_angle);
phi = params{1}(5); %trishear angle; assumed symmetric
phi = phi*pi/180; %convert to radians
m = tan(phi); %for use in trishear equations
if phi == pi/2;
    disp('Error: phi = 90 degrees, not allowed, using 89 degrees')
    m = tan(89*pi/180);
end
PoverS = params{1}(6); %propagation / slip
s = params{1}(7); % s term in velocity field equation, 1 for linear
if FaultType == 2 %Detachment
    decol_y = params{1}(8);
    nextra = 1;
end
if FaultType == 3 %Fault with a bend
    ramp_angle = [ramp_angle,(pi/180)*params{1}(8)];
    ramp_acute = [ramp_acute,min(ramp_angle(2),pi-ramp_angle(2))];
    bend_y = params{1}(9);
    nextra = 2;
end
slip_sense = params{1}(8+nextra); %1 for reverse or -1 for normal
increment = params{1}(9+nextra); %slip increment in which trishear is applied
total_slip = total_slip*slip_sense;
increment = increment*slip_sense;
v0 = increment; %same thing, just phrased as a velocity

%Find out what tip tipx and tipy refer to
TipType = InputGoodOnly('Tip Position is: (1) Initial, or (2) Final ',[1,2]);
if (FaultType == 1 || FaultType == 2) %Straight ramp or detachment
    if TipType == 1 %Initial
        tipx_init = tip;
        tipy_init = tipy;
        tipx_final = tipx+total_slip*PoverS*cos(ramp_angle);
        tipy_final = tipy+total_slip*PoverS*sin(ramp_angle);
    else %Final
        tipx_final = tipx;    
        tipy_final = tipy;
        tipx_init = tipx-total_slip*PoverS*cos(ramp_angle);
        tipy_init = tipy-total_slip*PoverS*sin(ramp_angle);
    end
    %If there's a detachment, calculate the x coordinate of where it intersects
    if FaultType == 2
        decol_x = tipx_final-(tipy_final-decol_y)/tan(ramp_angle); %x coordinate of where decollement meets ramp
    end
else %bend
    %Determine which segment the tip (tipx,tipy) is in. Note: this can be initial or final tip depending which one tipx,tipy represent.
    if (tipy >= bend_y)
        tip_seg = 1;
    else
        tip_seg = 2;
    end
    %Calculate the bend x coordinate
    bend_x = tipx-(tipy-bend_y)/tan(ramp_angle(tip_seg)); %x coordinate of where decollement meets ramp
    %Calculate initial and final tipx and tipy (final = present-day, initial = when trishear motion began). Final must be higher than initial.
    if (TipType == 1) %If tipx, tipy are initial
        tipx_init = tipx;
        tipy_init = tipy;
        tip_seg_init = tip_seg;
        %Note: If we allowed tipxfinal to be lower than tipxinit, this wouldn't be so simple. Same for if TipToSolve = 2 below
        if (tip_seg==1) %Starts and stays in upper segment
            tipx_final = tipx+abs(total_slip)*PoverS*cos(ramp_angle(1)); %Remember that if doing an inversion, total_slip and slip_sense will be negative
            tipy_final = tipy+abs(total_slip)*PoverS*sin(ramp_angle(1));
            tip_seg_final = 1;
        else
            dist = sqrt((tipy-bend_y)^2+(tipx-bend_x)^2);
            if (abs(total_slip)*PoverS <= dist) %Starts and stays in lower segment
                tipx_final = tipx+abs(total_slip)*PoverS*cos(ramp_angle(2));
                tipy_final = tipy+abs(total_slip)*PoverS*sin(ramp_angle(2));
                tip_seg_final = 2;
            else %Goes from lower into upper segment
                tipx_final = tipx+dist*cos(ramp_angle(2)) + (abs(total_slip)*PoverS-dist)*cos(ramp_angle(1));
                tipy_final = tipy+dist*sin(ramp_angle(2)) + (abs(total_slip)*PoverS-dist)*sin(ramp_angle(1));
                tip_seg_final = 1;
            end
        end
    else %If tipx, tipy are final
        tipx_final = tipx;
        tipy_final = tipy;
        tip_seg_final = tip_seg;
        if (tip_seg == 2) %Starts and stays in lower segment
            tipx_init = tipx-abs(total_slip)*PoverS*cos(ramp_angle(2));
            tipy_init = tipy-abs(total_slip)*PoverS*sin(ramp_angle(2));
            tip_seg_init = 2;
        else
            dist = sqrt((tipy-bend_y)^2+(tipx-bend_x)^2); %Distance tip will have to propagate to change segments.
            if (abs(total_slip)*PoverS <= dist) %Starts and stays in upper segment
                tipx_init = tipx-abs(total_slip)*PoverS*cos(ramp_angle(1));
                tipy_init = tipy-abs(total_slip)*PoverS*sin(ramp_angle(1));
                tip_seg_init = 1;
            else %tip_init in lower segment, tip_final in upper segment
                tipx_init = tipx-dist*cos(ramp_angle(1)) - (abs(total_slip)*PoverS-dist)*cos(ramp_angle(2));
                tipy_init = tipy-dist*sin(ramp_angle(1)) - (abs(total_slip)*PoverS-dist)*sin(ramp_angle(2));
                tip_seg_init = 2;
            end
        end
    end
end


%Find out what type of fit to do
FitType = InputGoodOnly('Type of fit:\n(1) Fit to flat line.\n(2) Fit to known dip.\n(3) Fit to known line. ',[1,2,3]);
switch FitType
    case 1 %Flat Line
        regional_dip_deg = 0;
        regional_dip = 0;
        regional_slope = 0;
    case 2 %Known Dip
        regional_dip_deg = input('Enter regional dip (degrees): (Negative is down to the left, positive right) ');
        regional_dip = regional_dip_deg*pi/180;
        regional_slope = tan(regional_dip);
    case 3 %Known Line
        regional_dip_deg = input('Enter regional dip (degrees): (Negative is down to the left, positive right) ');
        regional_dip = regional_dip_deg*pi/180;
        regional_slope = tan(regional_dip);
        if DatTypes(1) == 1 %Has beds
            intercepts = [];
            while length(intercepts) ~= nbeds
                intercepts = input(['Enter all ',num2str(nbeds),' bed intercepts as a vector: ']);
                if length(intercepts) ~= nbeds
                    disp('Error: Wrong number of beds')
                end
            end
        end
end

%Restore to flat
if FaultType == 2 %If there's a detachment, calculate decolze.
    decolze = xy_to_ze2([decol_x;decol_y],tipx_final,tipy_final,ramp_angle);
elseif FaultType == 3 %Or for a bend, calculate bendze
    bendze = xy_to_ze2([bend_x;bend_y],tipx_final,tipy_final,ramp_angle(tip_seg_final));
end
if DatTypes(1) == 1 %for bedding data
    %Calculate trishear velocity field and move beds accordingly
    if FaultType == 1 %Straight ramp
        bedsze = xy_to_ze2(beds,tipx_final,tipy_final,ramp_angle);
        bedsze = trishear_func4(bedsze,-v0,m,-total_slip,-increment,PoverS,s);
        beds_restored = ze_to_xy2(bedsze,tipx_final,tipy_final,ramp_angle);
    elseif FaultType == 2 %Detachment
        bedsze = xy_to_ze2(beds,tipx_final,tipy_final,ramp_angle);
        bedsze = trishear_func_decol4(bedsze,-v0,m,...
            -total_slip,-increment,PoverS,s,decolze,ramp_acute);
        beds_restored = ze_to_xy2(bedsze,tipx_final,tipy_final,ramp_angle);
    else %Bend
        bedsze = xy_to_ze2(beds,tipx_final,tipy_final,ramp_angle(tip_seg_final));
        bedsze = trishear_func_bend(bedsze,-v0,m,-total_slip,-increment,...
            PoverS,s,bendze,ramp_acute);
        beds_restored = ze_to_xy2(bedsze,tipx_init,tipy_init,ramp_angle(tip_seg_init));
    end
    % Compare beds to a straight line
    sq_errors = [];
    bed_depths = zeros(size(bedsze,3),1); %Holds the depth to the beds (for flat beds) or their intercepts (for dipping beds)
    for i = 1:size(bedsze,3)
    	bed_to_fit = beds_restored(:,:,i); %take just the bed being fit
    	bed_to_fit_x = bed_to_fit(1,isnan(bed_to_fit(1,:))==0);%remove NaN
    	bed_to_fit_y = bed_to_fit(2,isnan(bed_to_fit(2,:))==0);%remove NaN                       
        if FitType == 1 %Flat Line
            ymean = mean(bed_to_fit_y); %flat line has the form y = mean of y vals
            error = abs(bed_to_fit_y-ymean);
            bed_depths(i) = ymean;
        elseif FitType == 2 %Known regional dip
            bedlength = size(bed_to_fit,2); %Length of the bed
            b = sum(bed_to_fit_y+regional_slope*bed_to_fit_x)/bedlength; %Calculate y intercept of best fit line with RegSlope as its slope
            error = abs(bed_to_fit_y+regional_slope*bed_to_fit_x-b)/sqrt(regional_slope^2+1); %Errors (dist from pt to line)
            bed_depths(i) = b;
        elseif FitType == 3 %Known line
            b = intercepts(i); %Entered by user above.
            error = abs(bed_to_fit_y+regional_slope*bed_to_fit_x-b)/sqrt(regional_slope^2+1); %Errors (dist from pt to line)
            bed_depths(i) = b;
        end
        sq_errors = [sq_errors,error.^2]; %mean squared errors for each bed
    end
    RMS_beds_restored = sqrt(mean(sq_errors)); %root mean squared errors for all beds
end
if DatTypes(2) == 1 %For dip data
    if FaultType == 1 %Straight Fault
        dipposze = xy_to_ze2(dip_pos,tipx_final,tipy_final,ramp_angle); %dip positions in zeta, eta
        dipsze = (ramp_dir*rad_dips+ramp_acute); %dips in ze coord. system
        [dipposze,modeldipsze] = trishear_func_dip(dipposze,dipsze,-v0,m,-total_slip,-increment,PoverS,s);
        dippos_restored = ze_to_xy2(dipposze,tipx_final,tipy_final,ramp_angle);
        modeldips = ramp_dir*(modeldipsze - ramp_acute); %convert back to world coordinate system
    elseif FaultType == 2 %Decollement
        dipposze = xy_to_ze2(dip_pos,tipx_final,tipy_final,ramp_angle); %dip positions in zeta, eta
        dipsze = (ramp_dir*rad_dips+ramp_acute); %dips in ze coord. system
        [dipposze,modeldipsze] = trishear_func_dip_decol(dipposze,dipsze,-v0,m,-total_slip,-increment,PoverS,s,decolze,ramp_acute);
        dippos_restored = ze_to_xy2(dipposze,tipx_final,tipy_final,ramp_angle);
        modeldips = ramp_dir*(modeldipsze - ramp_acute); %convert back to world coordinate system
    else %Bend
        dipposze = xy_to_ze2(dip_pos,tipx_final,tipy_final,ramp_angle(tip_seg_final)); %dip positions in zeta, eta
        dipsze = (ramp_dir*rad_dips+ramp_acute(tip_seg_final)); %dips in ze coord. system
        [dipposze,modeldipsze] = trishear_func_dip_bend(dipposze,dipsze,-v0,m,-total_slip,-increment,PoverS,s,bendze,ramp_acute);
        dippos_restored = ze_to_xy2(dipposze,tipx_init,tipy_init,ramp_angle(tip_seg_init));
        modeldips = ramp_dir*(modeldipsze - ramp_acute(tip_seg_init)); %convert back to world coordinate system
    end
    modeldips = modeldips*180/pi; %convert to degrees
    modeldips(modeldips>90) = modeldips(modeldips>90)-180; %fix dips > 90
    modeldips(modeldips<-90) = modeldips(modeldips<-90)+180; %fix dips < -90
    %Flat or known regional dip
    R = abs(modeldips-regional_dip_deg); %residuals    
    R(abs(R)>90) = 180-abs(R(abs(R)>90));
    RMS_dips_restored = sqrt(sum(R.^2)/length(R));
end

%Calculate Fault Point Errors (No separate forward or reverse.)
if DatTypes(3) == 1 %Fault point data
    switch FaultType
        case(1) %Straight ramp
            errors = abs(((fault_pts(1,:)-tipx)*tan(ramp_angle)-fault_pts(2,:)+tipy)/sqrt(1+tan(ramp_angle)^2));
        case(2) %Detachment
            errors = zeros(1,size(fault_pts,2));
            for i = 1:size(fault_pts,2)
                if fault_pts(2,i)>-ramp_dir*(fault_pts(1,i)-decol_x)*tan((pi-ramp_acute)/2)+decol_y %Closer to ramp
                    errors(i) = abs(((fault_pts(1,i)-tipx)*tan(ramp_angle)-fault_pts(2,i)+tipy)/sqrt(1+tan(ramp_angle)^2));
                else %Closer to flat
                    errors(i) = abs(fault_pts(2,i)-decol_y);
                end
            end
        case(3) %Bend
            errors = zeros(1,size(fault_pts,2));
            for i = 1:size(fault_pts,2)
                if fault_pts(2,i)>-ramp_dir*(fault_pts(1,i)-bend_x)*tan((pi-ramp_acute(1)-ramp_acute(2))/2)+bend_y %Closer to upper ramp
                    errors(i) = abs(((fault_pts(1,i)-bend_x)*tan(ramp_angle(1))-fault_pts(2,i)+bend_y)/sqrt(1+tan(ramp_angle(1))^2));
                else %Closer to lower ramp
                    errors(i) = abs(((fault_pts(1,i)-bend_x)*tan(ramp_angle(2))-fault_pts(2,i)+bend_y)/sqrt(1+tan(ramp_angle(2))^2));
                end
            end
        case(default)
            disp('Error: Incorrect FaultType')
    end
    RMS_fault = sqrt(mean(errors.^2));
end

%Forward Model based on inverse results
if FaultType == 2 %If there's a detachment, calculate decolze.
    decolze = xy_to_ze2([decol_x;decol_y],tipx_init,tipy_init,ramp_angle);
elseif FaultType == 3 %Or for a bend, calculate bendze
    bendze = xy_to_ze2([bend_x;bend_y],tipx_init,tipy_init,ramp_angle(tip_seg_init));
end
if DatTypes(1) == 1 %for bedding data
    %Find the closest point on the restored bed to each restored bed point
    closest_pts = zeros(size(beds_restored));
    if FitType == 1 %Flat Beds
        for i = 1:size(beds,3)
            closest_pts(1,:,i) = beds_restored(1,:,i);
            closest_pts(2,:,i) = bed_depths(i);
        end
        closest_pts(isnan(beds_restored)==1) = NaN; %The ones that were NaN before should still be NaN.
    else %Regional Dip
        for i = 1:size(beds,3)
            closest_pts(1,:,i) = ((beds_restored(1,:,i)-regional_slope*beds_restored(2,:,i))...
                -regional_slope*bed_depths(i))/(regional_slope^2+1); %This seems to be wrong.
            closest_pts(2,:,i) = (regional_slope*(-beds_restored(1,:,i)+regional_slope*beds_restored(2,:,i))...
                +bed_depths(i))/(regional_slope^2+1);
        end
    end
    %Calculate trishear velocity field and move beds accordingly
    if FaultType == 1 %Straight ramp
        bedsze = xy_to_ze2(closest_pts,tipx_init,tipy_init,ramp_angle);
        bedsze = trishear_func4(bedsze,v0,m,total_slip,increment,PoverS,s);
        beds_forward = ze_to_xy2(bedsze,tipx_init,tipy_init,ramp_angle);
    elseif FaultType == 2 %Detachment
        bedsze = xy_to_ze2(closest_pts,tipx_init,tipy_init,ramp_angle);
        bedsze = trishear_func_decol4(bedsze,v0,m,...
            total_slip,increment,PoverS,s,decolze,ramp_acute);
        beds_forward = ze_to_xy2(bedsze,tipx_init,tipy_init,ramp_angle);
    else %Fault with a bend in it
        bedsze = xy_to_ze2(closest_pts,tipx_init,tipy_init,ramp_angle(tip_seg_init));
        bedsze = trishear_func_bend(bedsze,v0,m,total_slip,increment,...
            PoverS,s,bendze,ramp_acute);
        beds_forward = ze_to_xy2(bedsze,tipx_final,tipy_final,ramp_angle(tip_seg_final));
    end
    % Compare beds to data
    sq_errors = [];
    for i = 1:size(bedsze,3)
        bed_to_fit = beds_forward(:,:,i); %take just the bed being fit
    	bed_to_fit_x = bed_to_fit(1,isnan(bed_to_fit(1,:))==0);%remove NaN
    	bed_to_fit_y = bed_to_fit(2,isnan(bed_to_fit(2,:))==0);%remove NaN  
        bed_orig = beds(:,:,i); %take just the bed being fit
    	bed_orig_x = bed_orig(1,isnan(bed_orig(1,:))==0);%remove NaN
    	bed_orig_y = bed_orig(2,isnan(bed_orig(2,:))==0);%remove NaN
        errors = sqrt((bed_to_fit_x-bed_orig_x).^2+(bed_to_fit_y-bed_orig_y).^2);
        sq_errors = [sq_errors,errors.^2]; %mean squared errors for each bed
    end
    RMS_beds_forward = sqrt(mean(sq_errors)); %root mean squared errors for all beds
end
if DatTypes(2) == 1 %For dip data
    if FaultType == 1 %Straight Ramp
        dipposze = xy_to_ze2(dippos_restored,tipx_init,tipy_init,ramp_angle); %dip positions in zeta, eta
        dipsze = ones(size(dips))*(ramp_dir*regional_dip*pi/180+ramp_acute); %dips in ze coord. system
        [dipposze,modeldipsze] = trishear_func_dip(dipposze,dipsze,v0,m,total_slip,increment,PoverS,s);
        modeldips = ramp_dir*(modeldipsze - ramp_acute); %convert back to world coordinate system
    elseif FaultType == 2 %Fault with Detachment
        dipposze = xy_to_ze2(dippos_restored,tipx_init,tipy_init,ramp_angle); %dip positions in zeta, eta
        dipsze = ones(size(dips))*(ramp_dir*regional_dip*pi/180+ramp_acute); %dips in ze coord. system
        [dipposze,modeldipsze] = trishear_func_dip_decol(dipposze,dipsze,v0,m,total_slip,increment,PoverS,s,decolze,ramp_acute);
        modeldips = ramp_dir*(modeldipsze - ramp_acute); %convert back to world coordinate system
    else %Bend
        dipposze = xy_to_ze2(dippos_restored,tipx_init,tipy_init,ramp_angle(tip_seg_init)); %dip positions in zeta, eta
        dipsze = ones(size(dips))*(ramp_dir*regional_dip*pi/180+ramp_acute(tip_seg_init)); %dips in ze coord. system
        [dipposze,modeldipsze] = trishear_func_dip_bend(dipposze,dipsze,v0,m,total_slip,increment,PoverS,s,bendze,ramp_acute);
        modeldips = ramp_dir*(modeldipsze - ramp_acute(tip_seg_final)); %convert back to world coordinate system
    end
    modeldips = modeldips*180/pi; %convert to degrees
    modeldips(modeldips>90) = modeldips(modeldips>90)-180; %fix dips > 90
    modeldips(modeldips<-90) = modeldips(modeldips<-90)+180; %fix dips < -90
    %Flat or known regional dip
    R = abs(modeldips-dips); %residuals    
    R(abs(R)>90) = 180-abs(R(abs(R)>90));
    RMS_dips_forward = sqrt(sum(R.^2)/length(R));
end

%Create forward modelled beds to plot or save
if DatTypes(1) == 1
    %make the flat beds
    disp('Create forward modeled beds:')
    xmin = 1;
    xmax = 0;
    while xmin>xmax
        xmin = input(   'Left boundary: ');
        xmax = input(   'Right boundary: ');
        if xmin>xmax
            disp('Error: Right boundary must be greater than left.')
        end
    end
    step = input(   'x interval between points: ');
    model_beds = zeros(2,(xmax-xmin)/step+1,size(beds,3));
    for n = 1:size(beds,3)
        model_beds(1,:,n) = xmin:step:xmax;
        model_beds(2,:,n) = bed_depths(n)-regional_slope*(xmin:step:xmax); %- because of sign convention on regional slope
    end
    if FaultType == 2 %If there's a detachment, calculate decolze.
        decolze = xy_to_ze2([decol_x;decol_y],tipx_init,tipy_init,ramp_angle);
    end
    if FaultType == 2 %If there's a detachment, calculate decolze.
        decolze = xy_to_ze2([decol_x;decol_y],tipx_init,tipy_init,ramp_angle);
    elseif FaultType == 3 %Or for a bend, calculate bendze
        bendze = xy_to_ze2([bend_x;bend_y],tipx_init,tipy_init,ramp_angle(tip_seg_init));
    end
    %fold the beds
    %Calculate trishear velocity field and move beds accordingly
    if FaultType == 1 %Straight ramp
        model_bedsze = xy_to_ze2(model_beds,tipx_init,tipy_init,ramp_angle);
        model_bedsze = trishear_func4(model_bedsze,v0,m,total_slip,increment,PoverS,s);
        model_beds = ze_to_xy2(model_bedsze,tipx_init,tipy_init,ramp_angle);
    elseif FaultType == 2 %Detachment
        model_bedsze = xy_to_ze2(model_beds,tipx_init,tipy_init,ramp_angle);
        model_bedsze = trishear_func_decol4(model_bedsze,v0,m,...
            total_slip,increment,PoverS,s,decolze,ramp_acute);
        model_beds = ze_to_xy2(model_bedsze,tipx_init,tipy_init,ramp_angle);
    else %Bend
        model_bedsze = xy_to_ze2(model_beds,tipx_init,tipy_init,ramp_angle(tip_seg_init));
        model_bedsze = trishear_func_bend(model_bedsze,v0,m,total_slip,increment,...
            PoverS,s,bendze,ramp_acute);
        model_beds = ze_to_xy2(model_bedsze,tipx_final,tipy_final,ramp_angle(tip_seg_final));
    end
end

%Let the user choose things to do.
while 1
donext = InputGoodOnly(['Options: \n(1) Display RMS errors.\n(2) Display restored bed elevations.\n'...
    '(3) Plot forward modeled beds.\n(4) Plot restored beds\n'...
    '(5) Save forward modeled beds as text file for Move.\n(6) Quit\n'],[1:6]);
switch donext
    case 1 %RMS errors
        disp('Restored RMS Errors')
        if DatTypes(1) == 1
            disp(['   Beds: ',num2str(RMS_beds_restored)])
        end
        if DatTypes(2) == 1
            disp(['   Dips: ',num2str(RMS_dips_restored)])
        end
        disp('Forward RMS Errors')
        if DatTypes(1) == 1
            disp(['   Beds: ',num2str(RMS_beds_forward)])
        end
        if DatTypes(2) == 1
            disp(['   Dips: ',num2str(RMS_dips_forward)])
        end
        if DatTypes(3) == 1
            disp(['Fault Points Errors: ',num2str(RMS_fault)])
        end
    case 2 %Restored bed elevations
        if FitType == 1 %Flat
            disp('Restored Bed Elevations:')
        else %Regional dip
            disp('Restored Bed y Intercepts:')
        end
        for i = 1:size(beds,3)
            disp(['Bed ',num2str(i),': ',num2str(bed_depths(i))])
        end
    case 3 %Plot forwarded modeled beds
        %plot the trishear zone
        x = 1:10:1e2;
        if FaultType ~= 3 %Straight or Decol
            plot(x*cos(ramp_angle)-m*x*sin(ramp_angle)+tipx_final,x*sin(ramp_angle)+...
                m*x*cos(ramp_angle)+tipy_final,'g',x*cos(ramp_angle)+m*x*sin(ramp_angle)...
                +tipx_final,x*sin(ramp_angle)-m*x*cos(ramp_angle)+tipy_final,'g')
        else %Bend
            plot(x*cos(ramp_angle(tip_seg_final))-m*x*sin(ramp_angle(tip_seg_final))+...
                tipx_final,x*sin(ramp_angle(tip_seg_final))+m*x*cos(ramp_angle(tip_seg_final))+...
                tipy_final,'g',x*cos(ramp_angle(tip_seg_final))+m*x*sin(ramp_angle(tip_seg_final))...
                +tipx_final,x*sin(ramp_angle(tip_seg_final))-m*x*cos(ramp_angle(tip_seg_final))+tipy_final,'g')
        end
        xlabel('horizontal distance')
        ylabel('vertical distance')
        title('Forward Model from Restored Bed Elevations')
        hold on
        %Plot the beds
        faultslope = (tipy_final-tipy_init)/(tipx_final-tipx_init);
        faultintercept = tipy_final-faultslope*tipx_final;
        for i = 1:size(beds,3)
            ybed = model_beds(2,:,i);
            xbed = model_beds(1,:,i);
            if tipy_final < min(beds(2,:,i))
                plot(xbed,ybed)
            else
                hanging = model_beds(:,ybed>xbed*faultslope+faultintercept,i); %This doesn't work quite right for faults with non-horizontal bends or for footwall points below a detachment.
                foot = model_beds(:,ybed<xbed*faultslope+faultintercept,i);
                plot(foot(1,:),foot(2,:),'b',hanging(1,:),hanging(2,:),'b')
            end
        end
        %plot the fault
        switch FaultType
            case(1) %Straight Fault
                ylimits = get(gca,'YLim');
                base = min([ylimits,tipy_init,tipy_final]);
                ramp_base = tipx_init-((tipy_init-base)/tan(ramp_angle));
                plot([ramp_base,tipx_init,tipx_final],[base,tipy_init,tipy_final],'r')
            case(2) %Fault With Detachment
                xlimits = get(gca,'XLim'); %This will return two numbers.
                if ramp_dir == 1
                    edge = min([xlimits,decol_x]);
                else
                    edge = max([xlimits,decol_x]);
                end
                plot([edge,decol_x,tipx_final],[decol_y,decol_y,tipy_final],'r')
            case(3) %Fault with a bend
                xlimits = get(gca,'XLim'); %This will return two numbers.
                if ramp_dir == 1
                    edge_x = min([xlimits,bend_x]);
                else
                    edge_x = max([xlimits,bend_x]);
                end
                edge_y = (edge_x-bend_x)*tan(ramp_angle(2))+bend_y;
                if tipy_final > bend_y
                    plot([edge_x,bend_x,tipx_final],[edge_y,bend_y,tipy_final],'r')
                else
                    plot([edge_x,tipx_final],[edge_y,tipy_final],'r')
                end
        end
        hold off
    case 4 %Plot restored beds
        %plot the trishear zone
        x = 1:10:1e2;
        if FaultType ~= 3 %Straight or Decol
            plot(x*cos(ramp_angle)-m*x*sin(ramp_angle)+tipx_init,x*sin(ramp_angle)+...
                m*x*cos(ramp_angle)+tipy_init,'g',x*cos(ramp_angle)+m*x*sin(ramp_angle)...
                +tipx_init,x*sin(ramp_angle)-m*x*cos(ramp_angle)+tipy_init,'g')
        else %Bend
            plot(x*cos(ramp_angle(tip_seg_init))-m*x*sin(ramp_angle(tip_seg_init))+...
                tipx_init,x*sin(ramp_angle(tip_seg_init))+m*x*cos(ramp_angle(tip_seg_init))+...
                tipy_init,'g',x*cos(ramp_angle(tip_seg_init))+m*x*sin(ramp_angle(tip_seg_init))...
                +tipx_init,x*sin(ramp_angle(tip_seg_init))-m*x*cos(ramp_angle(tip_seg_init))+tipy_init,'g')
        end
        xlabel('horizontal distance')
        ylabel('vertical distance')
        title('Restored Bed Points')
        hold on
        %Plot the beds
        faultslope = (tipy_final-tipy_init)/(tipx_final-tipx_init);
        faultintercept = tipy_init-faultslope*tipx_init;
        for i = 1:size(beds,3)
            ybed = beds_restored(2,:,i);
            xbed = beds_restored(1,:,i);
            if tipy_init < min(beds(2,:,i))
                plot(xbed,ybed)
            else
                hanging = beds_restored(:,ybed>xbed*faultslope+faultintercept,i);
                foot = beds_restored(:,ybed<xbed*faultslope+faultintercept,i);
                plot(foot(1,:),foot(2,:),'b',hanging(1,:),hanging(2,:),'b')
            end
        end
        %plot the fault
        switch FaultType
            case(1) %Straight Fault
                ylimits = get(gca,'YLim');
                base = min([ylimits,tipy_init]);
                ramp_base = tipx_init-((tipy_init-base)/tan(ramp_angle));
                plot([ramp_base,tipx_init],[base,tipy_init],'r')
            case(2) %Fault With Detachment
                xlimits = get(gca,'XLim'); %This will return two numbers.
                if ramp_dir == 1
                    edge = min([xlimits,decol_x]);
                else
                    edge = max([xlimits,decol_x]);
                end
                plot([edge,decol_x,tipx_init],[decol_y,decol_y,tipy_init],'r')
            case(3) %Fault with a bend
                xlimits = get(gca,'XLim'); %This will return two numbers.
                if ramp_dir == 1
                    edge_x = min([xlimits,bend_x]);
                else
                    edge_x = max([xlimits,bend_x]);
                end
                edge_y = (edge_x-bend_x)*tan(ramp_angle(2))+bend_y;
                if tipy_init > bend_y
                    plot([edge_x,bend_x,tipx_init],[edge_y,bend_y,tipy_init],'r')
                else
                    plot([edge_x,tipx_init],[edge_y,tipy_init],'r')
                end
        end
        hold off
    case 5 %Save beds
        save_name = [input('Enter filename: ','s'),'.txt'];
        disp('Enter coordinates of left post for cross section: ')
        xpost = input('x (Easting): ');
        ypost = input('y (Northing): ');
        trend = input('Enter trend of line of section (in +x direction): ');
        trend_geom = 90-trend;%Convert from geographic to geometric coordinates.
        if trend_geom < -180
            trend_geom = 360+trend_geom;
        end
        trend_geom = trend_geom*pi/180;
        savefile = fopen(save_name,'w');
        for i = 1:size(beds,3)
            fprintf(savefile,'%s\r\n',[]);
            z = model_beds(2,:,i); %Elevation
            xsec = model_beds(1,:,i); %x in section coordinates
            x = xpost+xsec*cos(trend_geom); %Easting
            y = ypost+xsec*sin(trend_geom); %Northing
            for j = 1:length(x)
                fprintf(savefile,'%g	%g	%g\r\n',x(j),y(j),z(j));
            end
        end
        fprintf(savefile,'%s\r\n',[]);
        fclose(savefile);
    case 6
        return
end
end