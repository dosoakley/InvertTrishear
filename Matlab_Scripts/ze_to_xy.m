%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ beds ] = ze_to_xy( bedsze,tipx,tipy,ramp_angle )
%For Trishear, converts zeta, eta to x,y
%   Converts zeta, eta coordinates with zeta parallel to the fault and
%   eta perpendicular to it and the origin at the initial fault tip into 
%   standard x,y coordinates with the x axis horizontal and y vertical.
%   This version takes ramp_angle <= 90 and reflects if the ramp dips to
%   the left.

beds = zeros(size(bedsze)); %to hold the beds in x y coordinates
if ramp_angle > pi/2
    dip_dir = -1;
    ramp_angle = pi - ramp_angle;
else
    dip_dir = 1;
end
ramp_sin = sin(ramp_angle);
ramp_cos = cos(ramp_angle);
for i = 1:size(bedsze,3)
    %rotate
    for j=1:size(bedsze,2)
        beds(:,j,i) = [ramp_cos,-ramp_sin;...
            ramp_sin,ramp_cos]*bedsze(:,j,i);
    end
    %reflect
    if dip_dir == -1
        beds(1,:,i) = -beds(1,:,i);
    end
    %translate
    beds(:,:,i) = beds(:,:,i)+[tipx*ones(1,size(beds,2));...
        tipy*ones(1,size(beds,2))];

end
end

