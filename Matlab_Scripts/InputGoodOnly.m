%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ output ] = InputGoodOnly( words,goodvals,vargin )
%InputGoodOnly Only accepts good values
%   Runs input and then checks that the results are from a list of
%   acceptable results. If not, it asks for the input again. "words" is the
%   text to go into input, "goodvals" is the list of acceptable input
%   values, and "vargin", if present, must be s. Note that, if the input is
%   a string, then goodvals must be a cell array of acceptable strings.

good = 0;
while good == 0
    if nargin == 2
        output = input(words);
        if isempty(output) == 1
            good = 0;
        elseif output ~= goodvals
            good = 0;
            disp('Invalid Entry')
        else
            good = 1;
        end
    else
        output = input(words,vargin);
        for i = 1:length(goodvals)
            if length(output)==length(goodvals{i}) && max(output==goodvals{i})==1
                good = 1;
                break
            end
        end
        if good == 0
            disp('Invalid Entry')
        end
    end
end
    
end

