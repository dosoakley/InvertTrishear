%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ errors ] = ReadGridtxtErrs( filename,mins,maxs,steps )
%ReadGridtxtErrs Reads the output of the InvertTrishear program
%   Opens a textfile containing the RMS errors or probabilities produced by
%   the InvertTrishear program and converts them into an array to be used
%   by Matlab.
%   filename = name of the file being opened
%   mins = list of the minimum values of all model parameters
%   maxs = list of the maximum values of all model parameters
%   steps = list of the step sizes for all model parameters

file = fopen(filename);
if file == -1 %If the file doesn't exist, leave the function
    disp('Error: No such file.')
    return
end
errssize = (maxs-mins)./steps+1; %dimensions of the errors matrix
errssize = round(errssize); %In case of small discrepancies that occur at some point. A large discrepancy is a problem hopefully noted before now.
errors = zeros(errssize); %create the errors matrix
nmodels = numel(errors); %total number of models
nparams = length(mins); %total number of parameters
if nparams ~= length(maxs) || nparams ~= length(mins)
    disp('Error: Mins, Maxs, and Steps not all the same length')
    return
end
for i = 1:nmodels
    strerror = fgets(file); %Get this error as a string
    error = str2double(strerror); %Convert error to a number
    subsc = ones(1,nparams);%initialize subscript matrix
    kback = [1, cumprod(errssize(end:-1:2))];
    kback = fliplr(kback);
    ndx = i; %keeps track of remaining number
    for j = 1:nparams,
        vi = rem(ndx-1, kback(j)) + 1;         
        vj = (ndx - vi)/kback(j) + 1; 
        subsc(j) = vj; 
        ndx = vi;  
    end
    ind = 1; %Matlab index notation
    k = [1 cumprod(errssize(1:end-1))];
    for j = 1:nparams %convert subscript to index
        ind = ind + (subsc(j)-1)*k(j);
    end
    errors(ind) = error; %Put this error into the error matrix
end

fclose(file);
end

