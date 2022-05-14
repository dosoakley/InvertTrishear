%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ p ] = TrishearPDF2D_grid( p,steps,param1,param2,normalize )
%TRISHEAR PDF - Find 2D probability distribution of trishear parameters
%Attempts to fund the probability density function, or more
%accurately the probability over discrete ranges, for 2 trishear parameters
%using results of a trishear model run with Trishear. It differs
%from the simple histogram for a cutoff RMS value by instead using all the
%data and weighting by the RMS error in the data.
%
%p = multidimensional array of probability densities
%steps = values of steps used for trishear parameters
%param1 and param2 = parameters to use, 1 = tipx, 2 = tipy, etc.
%normalize = whether or not to normalize; 1 = yes, 0 = no

for j = 1:ndims(p)
    if (j ~= param1 && j ~= param2)
        siz = size(p,j); %size of the j dimension of p
        if siz ~= 1
            p = QuadIntegrate(p,j,steps(j));
        end
    end
end
%Normalize
if normalize == 1
    A = QuadIntegrate(p,param1,steps(param1));
    A = QuadIntegrate(A,param2,steps(param2));
    k = 1/A;
    p = p*k;
    p = squeeze(p);
elseif normalize ~= 0
    disp('Error: normalize option should be 1 or 0')
end


end

