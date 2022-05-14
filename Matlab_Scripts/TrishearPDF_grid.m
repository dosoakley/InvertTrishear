%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ p ] = TrishearPDF_grid( p,steps,param,normalize )
%TRISHEARPDF_grid - Find probability distribution of trishear parameters
%Attempts to find the probability density function, or more
%accurately the probability over discrete ranges, for a trishear parameter
%using results of a trishear model run with InvertTrishear. It differs
%from the simple histogram for a cutoff RMS value by instead using all the
%data and weighting by the RMS error in the data.
%
%p = multidimensional array of probability densities
%steps = values of steps used for trishear parameters
%param = parameter to use, 1 = tipx, 2 = tipy, etc.
%normalize = whether or not to normalize; 1 = yes, 0 = no

for j = 1:ndims(p)
    if j ~= param
        siz = size(p,j); %size of the j dimension of p
        if siz ~= 1
            %Integrate using Simpson quadrature
            p = QuadIntegrate(p,j,steps(j));
        end
    end
end
p = squeeze(p);
if param == 2
    p = p'; %So the data will be in the first dimension
end
%Normalize
if normalize == 1
    A = QuadIntegrate(p,1,steps(param));
    C = 1/A;
    p = p*C;
elseif normalize ~= 0
    disp('Error: normalize option should be 1 or 0')
end


end

