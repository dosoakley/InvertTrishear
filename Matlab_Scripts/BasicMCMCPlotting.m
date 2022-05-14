%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function  [output] = BasicMCMCPlotting(results,nbins)
%BasicMCMC Plotting
%Basic Analysis for MCMC run.
%Plots histograms and Markov chain paths.
%results = array containing the results of the MCMC run, as read in by
    %ReadMCtxtErrs.
%nbins = number of bins to use in each histogram.
close all
n = size(results,2);
for i = 2:n
    figure
    hist(results(:,i),nbins)
end
for i = 2:n
    figure
    plot(results(:,i))
end

disp('means: ')
for i = 2:n
    disp(mean(results(:,i)))
end

disp('standard deviations: ')
for i = 2:n
    disp(std(results(:,i)))
end

output = 'Finished';

end