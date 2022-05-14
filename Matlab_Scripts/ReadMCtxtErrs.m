%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ results ] = ReadMCtxtErrs(filename,n)
%Reads a text file of results from a Monte Carlo simulation as saved by the
%InvertTrishear program.
%filename = name of the file being opened
%n = number of entries per row in the file (number of parameters + 1)
%results = an array of the results.
file = fopen(filename);
if file == -1 %If the file doesn't exist, leave the function
    disp('Error: No such file.')
    return
end
format_str = '%f';
for i = 1:n-1
    format_str = [format_str,' %f'];
end
results_cell = textscan(file,format_str);
fclose(file);
results = zeros(size(results_cell{1},1),n);
for i = 1:n
    results(:,i) = results_cell{i};
end


end