%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ params ] = Params_from_Errs_Grid( errors,mins,steps,varargin )
%PARAMS_FROM_ERRS Recalculates parameters from errors
%  %This program recalculates model parameters for each run given the error
%matrix produced by the InvertTrishear program and the model parameters file
%used. It can also optionally convert only those specified by a list of
%indexes in 'errstosave'.
%errors = the multidimensional error matrix produced by a grid search.
%mins = the minimum values for each parameter
%steps = the step sizes for each parameter
%varargin = optional list of indixes for the errors grid, so as to save 
%   parameters only for those models.

%Calculate and save parameters with their errors
nparams = length(mins);
if nargin == 4
    errstosave = varargin{1};
    nsave = numel(errstosave);
    errvals = errors(errstosave); %error values, for output
else
    nsave = numel(errors);
    errvals = errors;
end
params = zeros(nsave,nparams+1); %matrix to hold the minima parameters
for i = 1:nsave
    if nargin == 3
        ind = i;
    else
        ind = errstosave(i);
    end
    subsc = ind2subMat(size(errors),ind);
    if length(subsc) < length(mins)
        subsc = [subsc,ones(1,length(mins)-length(subsc))];
    end
    params(i,:) = [mins+(subsc-1).*steps,errvals(i)];
end
end