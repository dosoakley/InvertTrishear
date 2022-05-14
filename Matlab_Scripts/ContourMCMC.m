%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ x,y,h ] = ContourMCMC( results,var1,var2,nbins,pl )
%ContourMCMC Contour plot of MCMC results
%   results = the results matrix from the MCMC run
%   var1 = the index of the first variable, tipx = 1, tipy = 2, etc
%   var2 = the index of the second variable
%   nbins = number of bins on each side
%   pl = whether or not to produce a plot (0 = no, 1 (or anything else) =
%       yes

h = zeros(nbins,nbins); %The histogram bins
data = [results(:,var1+1),results(:,var2+1)]; %The data for the two variables to plot
xmin = min(data(:,1));
ymin = min(data(:,2));
xmax = max(data(:,1));
ymax = max(data(:,2));
xrange = xmax-xmin; %range of x values
yrange = ymax-ymin; %range of y values
xstep = xrange/nbins; %step size in x direction
ystep = yrange/nbins; %step size in y direction
for i = 1:nbins
    for j = 1:nbins
        xl = xmin+(i-1)*xstep;
        xr = xmin+i*xstep;
        yl = ymin+(j-1)*ystep;
        yr = ymin+j*ystep;
        bin_data = data(find(data(:,1)>=xl & (data(:,1)<xr) |...
            (i==nbins & data(:,1)==xr)),:);
        bin_data = bin_data(find(bin_data(:,2)>=yl & (bin_data(:,2)<yr | ...
            (j == nbins & bin_data(:,2)==yr))),:);
        h(i,j) = size(bin_data,1);
    end
end
h = h'; %Since x should be the second dimension and y the first.
x = xmin+xstep/2:xstep:xmax-xstep/2;
y = ymin+ystep/2:ystep:ymax-ystep/2;
if pl ~= 0
    contour(x,y,h)
end


end

