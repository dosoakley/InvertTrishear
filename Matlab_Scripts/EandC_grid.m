%Matlab script written by David Oakley for use with the program
%InvertTrishear. If using in a publication, please acknowledge.

function [ E,C ] = EandC_grid( density,mins,maxs,steps )
%EandC Expectation and Covariance matrix for trishear parameters for a grid
%search
%   This function returns the expected value for each trishear
%   parameter and the covariance matrix, as calculated from equations 6.11
%   and 6.12 respectively of Tarantola and Valette (1982), "Inverse
%   Problems = Quest for Information."
%   density = Multidimensional array of probability densities for the grid search.
%   mins = minimum values of each parameter
%   maxs = maximum values of each parameter
%   steps = step sizes for each parameter.
%   E = expected values for each parameter.
%   C = covariance matrix.

nparams = ndims(density); %Total number of parameters
E = zeros(1,nparams); %Vector for expected values
C = zeros(nparams,nparams); %Vector for correlation coefficients

%1st calculate E and C(i,i) for each parameter
for i = 1:nparams
    if mins(i) == maxs(i) %If there was only one value tried
        E(i) = mins(i); %The only value for this not fit parameter
        C(i,i) = 0; %No variance
    else
        p = density;
        %Calculate the PDF for parameter i
        for k = 1:nparams
            if k ~= i
                siz = size(p,k); %size of the k dimension of p
                if siz ~= 1
                    %Integrate using Simpson quadrature
                    p = QuadIntegrate(p,k,steps(k));
                end
            end
        end
        p = squeeze(p);
        if i == 2
            p = p'; %So the data will be in the first dimension
        end
        %Normalize
        A = QuadIntegrate(p,1,steps(i));
        k = 1/A;
        p = p*k;
        %Calculate E(i) and C(i,i)
        vals = mins(i):steps(i):maxs(i); %the values of the parameter
        E(i) = QuadIntegrate(vals.*p',2,steps(i));
        C(i,i) = QuadIntegrate((vals.^2).*p',2,steps(i))-(E(i))^2; %Variance
    end
end

%Now calculate C for each pair of different parameters
for i = 1:nparams
    for j = i+1:nparams
        if (mins(i) == maxs(i)) || (mins(j) == maxs(j)) %If only one value        
            C(i,j) = 0; %No covariance 
            C(j,i) = 0;
        else            
            %Calculate the 2D pdf            
            p = density;            
            for k = 1:nparams            
                if k ~= i && k~=j                
                    siz = size(p,k); %size of the k dimension of p
                    if siz ~= 1
                        %Integrate using Simpson quadrature
                        p = QuadIntegrate(p,k,steps(k));
                    end
                end
            end
            %Normalize
            A = QuadIntegrate(p,i,steps(i));
            A = QuadIntegrate(A,j,steps(j));
            k = 1/A;
            p = p*k;
            p = squeeze(p);
            %Calculate C(i,j) and C(j,i)
            valsi = ones(size(p));%Values of parameter i
            for col = 1:size(p,2) %Fill in the columns of valsi
                valsi(:,col) = (mins(i):steps(i):maxs(i))';
            end
            valsj = mins(j):steps(j):maxs(j); %Values of parameter j
            Cpart = QuadIntegrate(valsi.*p,1,steps(i));
            Cpart = QuadIntegrate(valsj.*Cpart,2,steps(j));
            C(i,j) = Cpart - E(i)*E(j);
            C(j,i) = C(i,j);
        end
    end
end

end