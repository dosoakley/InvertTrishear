%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

%MCMC_plot_general.
%Plot figures of MCMC result histograms as subplots within a single figure.
%To use, first read in the results text file produced by InvertTrishear,
%assigning it to a matlab variable called "results", and save that variable
%as a .mat file. Then change the variables under the "Change these as
%needed" heading.

%This is for fault types 2 or 3 (horizontal detachment or fault with a
%single bend). For other options, it may require further modification.

%Change these as needed
filename = 'errs_apt.mat'; %The file name. You need to load your "results" variable as a matlab file.
nbins = 100; %number of bins to use in histograms.
burn = 2e4; %The number of models to remove as a burn-in period.
nparams = 7; %Number of parameters to plot
xdim = 4; %Number of plots in x direction
ydim = 2; %Number of plots in y direction
names = {'tip x', 'tip y','total slip','ramp angle (deg)','phi (deg)','P/S','detachment depth'}; %Parameter names
skip_s = 1; %Tells it to skip making a plot for the parameter s, if you kept that at 1. If not, change to 0, make nparams 8, and add s in names.

%You shouldn't have to change the rest
load(filename)
results = results(burn+1:end,:);
if skip_s == 1
    results2 = zeros(size(results,1),nparams+1);
    results2(:,1:7) = results(:,1:7);
    results2(:,8:end) = results(:,9:end);
    results = results2;
    clear('results2')
end
for i = 1:nparams
    subplot(ydim,xdim,i)
    hist(results(:,i+1),nbins)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','k')
    xlabel(names(i),'FontSize',12)
    if mod(i,xdim) == 1
        ylabel('n models','FontSize',12)
    end
end