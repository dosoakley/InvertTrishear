%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ bedsze ] = xy_to_ze( beds,tipx,tipy,ramp_angle )
%For Trishear, converts x,y to zeta, eta
%   Converts standard x,y coordinates with the x axis horizontal and y
%   vertical to zeta, eta coordinates with zeta parallel to the fault and
%   eta perpendicular to it and the origin at the initial fault tip.
%   This version takes ramp_angle <= 90 and reflects if the ramp dips to
%   the left.

bedsze = zeros(size(beds)); %to hold the beds in eta zeta coordinates
if ramp_angle > pi/2
    dip_dir = -1;
    ramp_angle = pi - ramp_angle;
else
    dip_dir = 1;
end
ramp_sin = sin(-ramp_angle);
ramp_cos = cos(-ramp_angle);
for i = 1:size(beds,3)
    %translate
    bedsze(:,:,i) = beds(:,:,i)-[tipx*ones(1,size(beds,2));...
        tipy*ones(1,size(beds,2))];
    %reflect
    if dip_dir == -1
        bedsze(1,:,i) = -bedsze(1,:,i);
    end
    %rotate
    for j=1:size(beds,2)
        bedsze(:,j,i) = [ramp_cos,-ramp_sin;...
            ramp_sin,ramp_cos]*bedsze(:,j,i);
    end
end

end

