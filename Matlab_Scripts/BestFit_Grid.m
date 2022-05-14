%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ best,subs ] = BestFit_Grid( errors )
%BestFit Finds the best fit trishear model from an errors matrix
%   Finds the minimum error trishear model from grid search results
%   and returns the indices of the best fitting model and its error value.
%   errors = the multidimensional array of results from a grid search.
%   best = error value (RMS, probability, or whatever was returend) for the
%       best fit model
%   subs = subscripts of the best fit model in the array.

best = min(min(min(min(min(min(min(min(errors))))))));
ind = find(errors == best);
if length(ind)> 1 % in case there is more than one equally good fit
    disp('More than one equally good fit.')
    if length(ind) <=50
        disp('Indices are:')
        disp(ind)
    else
        disp([num2str(length(ind)),' fits are equally good.'])
    end
    disp('Using first value only')
    ind = ind(1);
end
subs = ind2subMat(size(errors),ind); %subscripts
end