%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ z ] = QuadIntegrate( y,dim,step )
%QuadIntegrate Function for integration by quadrature.
%   %Uses Simpson quadrature to numerically integrate the probability density
%functions produced from trishear models (or really evenly spaced results
%of any multidimensional function). If the number of points is odd, the 
%last step is integrated using the trapezoid rule. To integrate across
%multiple dimensions, run the function repeatedly on the results of the
%previous integration.
% y = the matrix of values to be integrated
% dim = the dimension to integrate across.
% step = the step size between evaluation points. The distance between
%        endpts is 2*step
% z = the output matrix

perm = [dim:max(ndims(y),dim) 1:dim-1];%Matrix to speciy order of dimensions for permute
y = permute(y,perm);%Put the dimension to be integrated first
n = size(y,1); %number of points on the integration variable
endpts = 1:2:min(2*floor(n/2)+1,n); %end points of integration intervals
midpts = 2:2:endpts(end)-1; %mid points of integration intervals
if endpts(end) ~= n %An extra point, not in an interval, use trapezoidal rule
    trappts = [n-1,n];
else
    trappts = [];
end
steps = step*ones(1,length(midpts));
z = (2*steps/6)*(y(endpts(1:end-1),:)+4*y(midpts,:)+y(endpts(2:end),:))+step*0.5*sum(y(trappts,:)); %calculate the numerical integral using matrix multiplication
siz = size(y);
siz(1) = 1; %since that dimension has been integrated across
z = reshape(z,siz); %Above z = step made it 2 dimensional. Now return to multidimensional
z = ipermute(z,perm);
end

