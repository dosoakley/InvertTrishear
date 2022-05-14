%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ beds ] = ReadBeds( file )
%ReadBeds Read in bedding or contact data
%   This function reads a bedding or contact data file into Matlab,
%   producing a 3D matrix containing the points that define the beds.
%
    bed=[];
    beds = []; %for all the new beds that get loaded in
    i = 0; %counter, for points in a bed
    b = 1; %another counter, for beds
    while 1
        i = i+1;
        line = fgets(file); %get line
        if line ~= -1
            line = str2num(line); %convert to a number
        end
        if (isempty(line) == 0 && line(1)~=-1)
            bed(1,i) = line(1); %x coordinate
            bed(2,i) = line(end); %y coordinate
        else
            i = 0; %reset for a new bed
            if isempty(bed) == 0 %for some reason it sometimes reads empty beds
                if size(bed,2) > size(beds,2)
                    difference = size(bed,2)-size(beds,2);
                    beds = [beds,zeros(2,difference,size(beds,3))/0];
                elseif size(bed,2) < size(beds,2)
                    difference = size(beds,2)-size(bed,2);
                    bed = [bed,zeros(2,difference)/0];
                end
                beds(:,:,b) = bed;
                b = b+1; % on to a new bed
            end
            bed = []; %clear it for a new bed
        end
        if line == -1 %when you reach the end
            break
        end
    end


end

