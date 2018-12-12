
%%  Initialization    
clc
clearvars
close all
warning('off','MATLAB:lang:badlyScopedReturnValue')
warning('off','MATLAB:xlswrite:NoCOMServer')

%%  Determine input
filn        =   [pwd '/AE4423_Datasheets.xlsx'];
filn2       =   [pwd '/Group8_results.xlsx'];

Demand               =   xlsread(filn,'Group 8', 'C15:V34');
Airport_data         =   xlsread(filn,'Group 8', 'C6:V9');
Airport_data_incl_US =   xlsread(filn,'Group 8', 'C6:Z9');


Pop2010_incl_US    =   xlsread(filn,'General', 'B4:B27')/1000;
GDP2010_incl_US    =   xlsread(filn,'General', 'F4:F27');

Pop2017    =   xlsread(filn,'General', 'C4:C23')/1000;
Pop2017_US =   xlsread(filn,'General', 'C24:C27')/1000;
GDP2017    =   xlsread(filn,'General', 'G4:G23');
GDP2017_US =   xlsread(filn,'General', 'G24:G27');

fuelfactor2017 = 1.42;
fuelfactor2022 = 1.6;

%%  Determine distance between airports
 Airport_distance_incl_US = zeros(24,24); 
for i = 1:24
    for j = 1:24
        airportDistance_ij = arclen(i,j,Airport_data_incl_US);
        if airportDistance_ij ~= Inf
            Airport_distance_incl_US(i,j) = airportDistance_ij;
        end
    end
end

%Remove US values for parameter estimation with just European destinations
Airport_distance = Airport_distance_incl_US(1:20,1:20);


%%  Convert data to right format, removing rows with 0 demand (when flying from i to i)
ydata = log(reshape(Demand,[400,1]));
ydata(ydata==-Inf) = 1e-15;

popdata = log(reshape(Pop2017*transpose(Pop2017),[400,1]));
popdata(popdata==-Inf) = 1e-15;

GDPdata = log(reshape(GDP2017*transpose(GDP2017),[400,1]));
GDPdata(GDPdata==-Inf) = 1e-15;

distdata = log(fuelfactor2017*reshape(Airport_distance,[400,1]));
distdata(distdata==-Inf) = 1e-15;

xdata = [popdata GDPdata distdata];


%%  Linear regression to estimate the coefficients
mdl = fitlm(xdata,ydata,'linear','RobustOpts','on');
log_k = mdl.Coefficients.Estimate(1); %Natural logarithm of k
b1 = mdl.Coefficients.Estimate(2);
b2 = mdl.Coefficients.Estimate(3);
b3 = mdl.Coefficients.Estimate(4);

%%  Estimate demands for 2017 (Only European destinations)
demandEstimates2017 = zeros(20,20);
for i = 1:20
    for j = 1:20
        if i ~= j
            lndemandEstimates_ij = log_k + b1 * log(Pop2017(i)*Pop2017(j)) + b2 * log(GDP2017(i)*GDP2017(j)) + b3 * log(fuelfactor2017 * Airport_distance(i,j));
            demandEstimates2017(i,j) = exp(lndemandEstimates_ij);
        end
    end
end

%%  Estimate population and GPD for 2022, assuming linear variation (y = ax + b)
a_pop = ([Pop2017; Pop2017_US] - Pop2010_incl_US)/7;
Pop2022 = a_pop*12 + Pop2010_incl_US;

a_GDP = ([GDP2017; GDP2017_US] - GDP2010_incl_US)/7;
GDP2022 = a_GDP*12 + GDP2010_incl_US;

%%  Estimate European demands for 2020
demandEstimates2022 = zeros(24,24);
for i = 1:24
    for j = 1:24
        if i ~= j
            lndemandEstimates_ij = log_k + b1 * log(Pop2022(i)*Pop2022(j)) + b2 * log(GDP2022(i)*GDP2022(j)) + b3 * log(fuelfactor2022 * Airport_distance_incl_US(i,j));
            demandEstimates2022(i,j) = exp(lndemandEstimates_ij);
            if i > 20 || j > 20
               demandEstimates2022(i,j) = 10 * demandEstimates2022(i,j); 
            end
        end
    end
end


%%  Write ouput to excel file
xlswrite(filn2,round(demandEstimates2022),'Demands2022')

%% Function for determining the great circle distance between two airports
function out = arclen(airport_i,airport_j,Airport_data)
    delta_sigma = 2*asin(sqrt((sin(deg2rad((Airport_data(1,airport_i)-Airport_data(1,airport_j))/2)))^2+cos(deg2rad(Airport_data(1,airport_i)))*cos(deg2rad(Airport_data(1,airport_j)))*(sin(deg2rad((Airport_data(2,airport_i)-Airport_data(2,airport_j))/2)))^2));
    %[arclen, azimuth] = distance(Airport_data(1:2,(1:end-1)), Airport_data(1:2,(2:end)), spheroid);
    out = 6371*delta_sigma;
end

