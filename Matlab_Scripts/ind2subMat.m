%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ subs ] = ind2subMat( siz,ind )
%ind2subMat Like ind2sub, but with output as a single matrix
%   Converts index "ind" of a matrix of size "siz" into subscript notation,
%   and stores the results in a single matrix.

n = length(siz);
subs = zeros(1,n);
ndx = ind;
k = [1 cumprod(siz(1:end-1))];
for i = n:-1:1
  vi = rem(ndx-1, k(i)) + 1;         
  vj = (ndx - vi)/k(i) + 1; 
  subs(i) = vj; 
  ndx = vi;     
end

end