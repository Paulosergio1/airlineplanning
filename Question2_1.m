
%%  Initialization    
clc
clearvars
close all
warning('off','MATLAB:lang:badlyScopedReturnValue')
warning('off','MATLAB:xlswrite:NoCOMServer')

%%  Determine input
filn        =   [pwd '/AE4423_Datasheets.xlsx'];

Demand      =   xlsread(filn,'Group 8', 'C15:V34');
Airport_data=   xlsread(filn,'Group 8', 'C6:V9');

Pop2017    =   xlsread(filn,'General', 'C4:C23');
GDP2017    =   xlsread(filn,'General', 'G4:G23');

%%  Determine distance between airports
Airport_distance = zeros(20,20); 
for i = 1:20
    for j = 1:20
        airportDistance_ij = arclen(i,j,Airport_data);
        if airportDistance_ij ~= Inf
            Airport_distance(i,j) = airportDistance_ij;
        end
    end
end

%% Determine error using Non-Linear Least Squares https://nl.mathworks.com/help/optim/ug/nonlinear-curve-fitting-with-lsqcurvefit.html
xdata = [Pop2017 GDP2017 Airport_distance];
ydata = Demand;
a0 = [2 2 2 2];

hybridopts = optimset('MaxFunEvals', 1200);

predicted = @(a,xdata) a(1) * ((xdata(:,1)*transpose(xdata(:,1)).^(a(2))) * (xdata(:,2)*transpose(xdata(:,2)).^(a(3)))) / ((1.42 * xdata(:,3:22)).^(a(4)));
[ahat,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(predicted,a0,xdata,ydata,[-Inf;-Inf;-Inf;-Inf],[Inf;Inf;Inf;Inf],hybridopts);



%% Function for determining the demand between two airports    
function out = demand(airport_i,airport_j,Airport_data,Population,GDP,k,b1,b2,b3)
    out = k * (((Population(airport_i)*Population(airport_j))^b1 * (GDP(airport_i)*GDP(airport_j) )^b2) / (1.42 * arclen(airport_i,airport_j,Airport_data))^b3);
end


%% Function for determining the great circle distance between two airports
function out = arclen(airport_i,airport_j,Airport_data)
    delta_sigma = 2*asin(sqrt((sin(deg2rad((Airport_data(1,airport_i)-Airport_data(1,airport_j))/2)))^2+cos(deg2rad(Airport_data(1,airport_i)))*cos(deg2rad(Airport_data(1,airport_j)))*(sin(deg2rad((Airport_data(2,airport_i)-Airport_data(2,airport_j))/2)))^2));
    %[arclen, azimuth] = distance(Airport_data(1:2,(1:end-1)), Airport_data(1:2,(2:end)), spheroid);
    out = 6371*delta_sigma;
end


%% Estimate the demands
%parameters = [0.3 0.7 0.7 0.7];
%demandEstimates = zeros(20,20);
%demandDifferences = zeros(20,20);
%for i = 1:20
%    for j = 1:20
%        airportDistance_ij = demand(i,j,Airport_data,Pop2017,GDP2017,parameters(1),parameters(2),parameters(3),parameters(4));
%        if airportDistance_ij ~= Inf
%            demandEstimates(i,j) = airportDistance_ij;
%            demandDifferences(i,j) = Demand(i,j) - airportDistance_ij;
%        end
%    end
%end

