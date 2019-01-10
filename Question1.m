function Multicommodity ()
%%  Initialization
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
    savepath
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\examples\src\matlab');
    savepath
    clc
	clearvars
    close all
    warning('off','MATLAB:lang:badlyScopedReturnValue')
    warning('off','MATLAB:xlswrite:NoCOMServer')

    %%  Determine input
%   Select input file and sheet
    filn        =   [pwd '/AE4423_Datasheets.xlsx'];
    
    Demand      =   xlsread(filn,'Group 8', 'C15:V34');
    Airport_data=   xlsread(filn,'Group 8', 'C6:V9');
    [~,Airport_name] =   xlsread(filn,'Group 8', 'C5:Z5'); 
    
    fleet       =   xlsread(filn,'Group 8', 'B12:F12');
    
    Population  =   xlsread(filn,'General', 'B4:C27');
    GDOP        =   xlsread(filn,'General', 'F4:G27');
 
    ACData      =   xlsread(filn,'Group 8', 'B116:F124');
    
    Nodes      =    length(Demand(1,:));
    
    %%  Initiate CPLEX model
%   Create model
        model                   =   'Airline_planning';
        cplex                   =   Cplex(model);
        cplex.Model.sense       =   'maximize';

%       Number of aircraft types
        k_ac                       =   3;
        
%       Load factor
        LF = 0.75;
        
%       Variable g which determines if the airport is a hub or not
        g=ones(Nodes,1);
        g(3,1)=0;
        
%       Utilisation time in hours
        time_used=10*7;
        
        %Fuel price
        Fuel_price              = 1.42; %USD/gallon
        
%   Decision variables
    
    %%  Objective Function
        DV_Xnodes                  =   Nodes*Nodes; % number of nodes for direct passenger
        DV_Wnodes                  =   Nodes*Nodes; % number of nodes for transfer passenger
        DV_Znodes                  =   Nodes*Nodes*k_ac;% number of nodes for flights with ac type 

        %   Decision variables
        DV=DV_Xnodes+DV_Wnodes+DV_Znodes
        
        
        obj                     =   ones(DV,1);
        lb                      =   zeros(DV, 1);                                 % Lower bounds
        ub                      =   inf(DV, 1);                                   % Upper bounds
        ctype                   =   char(ones(1, (DV)) * ('I'));                  % Variable types 'C'=continuous; 'I'=integer; 'B'=binary
        
        l = 1;                                      % Array with DV names  (OPTIONAL, BUT HELPS READING THE .lp FILE)
        for i =1:Nodes % objective function values for direct passengers
            for j = 1:Nodes                    % of the x_{ij}^k variables
                if i == j
                    obj(l,1)      = 1e-15;
                else
                    obj(l,1)      = (5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043)*(arclen(i,j,Airport_data));
                end
                NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '   '];
                l = l + 1;
            end
        end
        
        for i =1:Nodes % objective function values for trasfer passenger
            for j = 1:Nodes                    % of the x_{ij}^k variables
                if i == j 
                    obj(l,1)      = 1e-15;
                else
                    obj(l,1)      = (5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043)*(arclen(i,j,Airport_data));
                end
                NameDV (l,:)  = ['W_' num2str(i,'%02d') ',' num2str(j,'%02d') '   '];
                l = l + 1;
            end
        end
        
        for k=1:k_ac % Objective function values for flights per ac type
            for i =1:Nodes
                for j = 1:Nodes                    % of the x_{ij}^k variables
                    if g(i)==0 || g(j)==0 %Economies of scale, for flights arriving or departing at hub
                       obj(l,1)    =   -0.7*(ACData(8,k)*arclen(i,j,Airport_data)/ACData(1,k)+ACData(9,k)*Fuel_price/(1.5)*arclen(i,j,Airport_data)+ACData(7,k));
                    else
                        obj(l,1)    =   -(ACData(8,k)*arclen(i,j,Airport_data)/ACData(1,k)+ACData(9,k)*Fuel_price/(1.5)*arclen(i,j,Airport_data)+ACData(7,k));
                    end
                    NameDV (l,:)  = ['z_' num2str(i,'%02d') ',' num2str(j,'%02d') ',' num2str(k, '%02d') ];
                    l = l + 1;
                end
            end
        end
       
        % cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
        cplex.addCols(obj, [], lb, ub, ctype, NameDV);
        cplex.writeModel([model '.lp']);
        
    %%  Fixed cost calculation
        fixed_cost=0;
        for k=1:k_ac
            fixed_cost=fixed_cost+fleet(k)*ACData(6,k);
        end
        
    %%  Constraints
        % Aurcraft cannot be used more than 10 hours a day, so 70 hours a
        % weekr
        for k=1:k_ac
            C_time_ac=zeros(1,DV);
            for i=1:Nodes
                for j=1:Nodes
                    distance=arclen(i,j,Airport_data);
                    TAT=ACData(3,k)/60; % Changing TAT from minutes to hours
                    if g(j)==0
                        TAT=max([1 TAT*2]); % TAT at the hub
                    end
                    C_time_ac(varindex(i,j,k,'z',Nodes))=distance/ACData(1,k)+TAT;
                end
            end
            rightvariable=time_used*fleet(k);
            cplex.addRows(0, C_time_ac, rightvariable, sprintf('Timeusedac%d',k));
        end
        
    %   Passengers not more than demand
        for i = 1:Nodes
            for j = 1:Nodes
                C_dem = zeros(1,DV);
                C_dem(varindex(i,j,0,'x',Nodes)) = 1;
                C_dem(varindex(i,j,0,'w',Nodes)) = 1;
                cplex.addRows(0, C_dem, Demand(i,j), sprintf('DemandConstraint%d_%d',i,j));
            end
        end
    %   No transfer if one the airports is hub 
        for i = 1:Nodes
            for j = 1:Nodes
                C_hub     =   zeros(1, DV);
                C_hub(varindex(i,j,0,'w',Nodes)) = 1;
                hub_or_not = g(i) * g(j);
                cplex.addRows(0, C_hub, Demand(i,j) * hub_or_not,sprintf('TransferHub%d_%d',i,j));
            end
        end
    
    %   Max aircraft capacity
        for i = 1:Nodes
            for j = 1:Nodes
                C_cap = zeros(1,DV);
                C_cap(varindex(i,j,0,'x',Nodes)) = 1;
                if g(j) == 0 || g(i) == 0
                    if i~=j
                        for m = 1:Nodes
                            C_cap(varindex(m,j,0,'w',Nodes)) = 1 - g(i);
                            C_cap(varindex(i,m,0,'w',Nodes)) = 1 - g(j);
                        end
                    end
                end
                for k = 1:k_ac
                    C_cap(varindex(i,j,k,'z',Nodes)) = -(ACData(2,k)*LF);
                end
                cplex.addRows(-inf, C_cap,0,sprintf('ACcapacity%d_%d',i,j));
            end
        end

        
        %Only aircraft with a long enough range can fly between two cities
        for k=1:k_ac
            for i=1:Nodes
                for j=1:Nodes
                    C_range=zeros(1,DV);
                    distance=arclen(i,j,Airport_data);
                    C_range(varindex(i,j,k,'z',Nodes))=1;
                    if distance>ACData(4,k)
                        cplex.addRows(0, C_range, 0, sprintf('range%d_%d_%d',i,j,k));
                    end
                end
            end
        end
        
       
        
       %Constraint for the runway length wich should be long enough. 
%         for k=1:k_ac
%             for i=1:Nodes
%                 C_runway=zeros(1,DV);
%                 if ACData(5,k)>Airport_data(3,i)
%                     for j=1:Nodes
%                         C_runway(varindex(i,j,k,'z',Nodes))=1;
%                     end
%                     cplex.addRows(0, C_runway, 0, sprintf('runway%d_%d',i,k));
%                 end
%             end
%         end
        
        % Slots contraint
%         for i=1:Nodes
%             C_slots=zeros(1,DV);
%             for j=1:Nodes
%                 for k=1:k_ac
%                    C_slots(varindex(i,j,k,'z',Nodes))=1;
%                 end
%             end
%             cplex.addRows(0, C_slots, Airport_data(4,i), sprintf('slots%d',i));
%         end
        
        %flow inside of airport should be equal to flow outside of airport
        for k=1:k_ac
            for i=1:Nodes
                C_flow=zeros(1,DV);
                for j=1:Nodes
                    if i~=j
                        C_flow(varindex(i,j,k,'z',Nodes))=1;
                        C_flow(varindex(j,i,k,'z',Nodes))=-1;
                    end
                end
                cplex.addRows(0, C_flow, 0, sprintf('flow%d_%d',i, k));
            end
        end
        
            
     %%  Execute model
        %cplex.Param.mip.limits.nodes.Cur    = 1e+8;         %max number of nodes to be visited (kind of max iterations)
        cplex.Param.timelimit.Cur           = 250;         %max time in seconds
        cplex.Param.mip.tolerances.mipgap.Cur   = 0.009;
        
        
%   Run CPLEX
        cplex.solve();
        cplex.writeModel([model '.lp']);
    
     %%  Postprocessing
%   Store direct results
    status                      =   cplex.Solution.status;
    sol.profit      =   cplex.Solution.objval - fixed_cost;
    sol.values      =   cplex.Solution.x;
    sol.PassengerDirect (:,:)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,'x', Nodes):varindex(Nodes, Nodes, k, 'x', Nodes)), Nodes, Nodes))';
    sol.PassengerIndirect (:,:)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,'w', Nodes):varindex(Nodes, Nodes, k, 'w', Nodes)), Nodes, Nodes))';
    for k = 1:k_ac
        sol.Flow (:,:,k)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,'z', Nodes):varindex(Nodes, Nodes, k, 'z', Nodes)), Nodes, Nodes))';
    end
    
    worldmap([35 65],[-15 30])
    land = shaperead('landareas.shp', 'UseGeoCoords', true);
    geoshow(land, 'FaceColor', [0.6 0.6 0.6])
    color = ['r','g','b','c','m'];
    hold on
    
    % Count the numbers of slots used
    slots=zeros(Nodes,1);
%   Write output
    fprintf('\n-----------------------------------------------------------------\n');
    fprintf ('Objective function value:          %10.1f  \n', sol.profit);
    fprintf ('\n') 
    fprintf ('Link From   To         AC1    AC2   AC3    Total (Revenue per Seat)  (Revenue) (ASK)  (RASK) (CASK) (RPK) (Profit) (Yield) (ALF) (BELF)\n');
    NL      =   0;
    alf_array = [];
    belf_array = [];
    profit_array = [];
    for i = 1:Nodes
        for j = 1:Nodes
            if sol.Flow(i,j,1)+sol.Flow(i,j,2)+sol.Flow(i,j,3)>0
                yield=obj(varindex(i,j,1,'x',Nodes))*sol.PassengerDirect(i,j);
                pass_transfer = 0;
                if g(i) == 0 && g(j) ~= 0
                    %sol.PassengerIndirect(i,j)
                    for l = 1:Nodes
                        %for m = 1:Nodes
                            if l == 3
                                continue;
                            end
                            yield = yield + obj(varindex(l,j,1,'w',Nodes))*sol.PassengerIndirect(l,j)*(arclen(3,j,Airport_data)/(arclen(3,j,Airport_data)+arclen(l,3,Airport_data)));
                            pass_transfer = pass_transfer + sol.PassengerIndirect(l,j);
                        %end
                    end
                elseif g(j) == 0 && g(i) ~= 0
                    %sol.PassengerIndirect(i,j)
                    %for l = 1:Nodes
                        for m = 1:Nodes
                            if m == 3
                                continue;
                            end
                            yield = yield + obj(varindex(i,m,1,'w',Nodes))*sol.PassengerIndirect(i,m)*(arclen(i,3,Airport_data)/(arclen(3,m,Airport_data)+arclen(i,3,Airport_data)));
                            pass_transfer = pass_transfer + sol.PassengerIndirect(i,m);
                        end
                    %end
                end
                profit = yield;
                for k= 1:k_ac
                    profit = profit + obj(varindex(i,j,k,'z',Nodes))*sol.Flow (i,j,k);
                end
                slots(i,1)=slots(i,1)+sol.Flow(i,j,1)+sol.Flow(i,j,2)+sol.Flow(i,j,3);
                NL      = NL + 1;
                revenue = obj(varindex(i,j,1,'x',Nodes))*(sol.Flow (i,j,1)*ACData(2,1) + sol.Flow (i,j,2)*ACData(2,2) + sol.Flow (i,j,3)*ACData(2,3));
                ask = arclen(i,j,Airport_data)*(sol.Flow (i,j,1)*ACData(2,1) + sol.Flow (i,j,2)*ACData(2,2) + sol.Flow (i,j,3)*ACData(2,3));
                rask = revenue*ask;
                cost_tot = revenue-profit;
                cask = cost_tot/ask;
                rpk = arclen(i,j,Airport_data)*(sol.PassengerDirect(i,j)+pass_transfer);
                %profit = rpk*yield - ask*cask;
                profit_array = [profit_array,profit];
                alf = (sol.PassengerDirect(i,j)+pass_transfer)/(sol.Flow (i,j,1)*ACData(2,1) + sol.Flow (i,j,2)*ACData(2,2) + sol.Flow (i,j,3)*ACData(2,3));
                alf_array = [alf_array,alf];
                belf = (cost_tot/revenue)*alf;
                belf_array = [belf_array,belf];
                fprintf (' %2d  %s   %s  %5d  %5d %5d %6d  %5d  %5d 5%d %5d %5d %5d %5d %5d %5d %5d\n', NL, Airport_name{i}, ...
                            Airport_name{j}, sol.Flow (i,j,1), sol.Flow (i,j,2), ...
                            sol.Flow (i,j,3), sol.Flow (i,j,1)+sol.Flow (i,j,2)+sol.Flow(i,j,3), ...
                            obj(varindex(i,j,1,'x',Nodes)), revenue, ask, rask, cask, rpk, profit, yield, alf, belf);
                for k=1:k_ac
                    if sol.Flow(i,j,k)>0
                        h = geoshow([Airport_data(1,i);Airport_data(1,j)],...
                                [Airport_data(2,i);Airport_data(2,j)]);
                        h.Marker = '*';
                        h.Color = color(k);
                        h.LineWidth = sol.Flow(i,j,k);
                    end
                end
%                 x_loc=(Airport_data(1,i)+Airport_data(1,j))/2;
%                 y_loc=(Airport_data(2,i)+Airport_data(2,j))/2;
%                 Capacity=ACData(2,1:k_ac)*[sol.values(varindex(i,j,1,'z', Nodes)),sol.values(varindex(i,j,2,'z', Nodes)),sol.values(varindex(i,j,3,'z', Nodes))]';
%                 textm(x_loc,y_loc,['(',num2str(Capacity), ')'])
                textm(Airport_data(1,i), Airport_data(2,i),Airport_name(i))

            end
        end
    end
    
    anlf = sum(alf_array)/NL
    nbelf = sum(belf_array)/NL
    tot_profit = sum(profit_array)- fixed_cost
    
    
    fprintf('\n------------------------Slots-------------------------------------\n');
    fprintf ('Used Available \n');
    for i=1:Nodes
        fprintf (' %2d       %5d     \n', slots(i,1), Airport_data(4,i));
    end
   
    fprintf('\n------------------------Cost per AC-------------------------------------\n');
    
    for k=1:k_ac
        cost=0;
        for i=1:Nodes
            for j=1:Nodes
                cost=cost+obj(varindex(i,j,k,'z', Nodes))*sol.values(varindex(i,j,k,'z', Nodes));
            end
        end
        fprintf ('Cost AC type:     %d  \n', k);
        fprintf ('            :     %10.1f  \n', cost);
    end
    
        
end
function out = varindex(i,j,k,letter,nodes)
    if letter == 'x'
        out=(i-1)*nodes+j;
    elseif letter == 'w'
        out=(i-1)*nodes+j+nodes^2;
    elseif letter == 'z'
        out=(i-1)*nodes+j+2*nodes^2+(k-1)*nodes^2;
    end   
        % Function given the variable index for each DV (i,j,k) the letter
        % denotes wheter you would like to have the variable x,w or z. 
end

%{
function out = arclen(airport_i,airport_j,Airport_data)
    spheroid = wgs84Ellipsoid;
    spheroid.SemimajorAxis = spheroid.SemimajorAxis;
    spheroid.SemiminorAxis = spheroid.SemiminorAxis;
    %[arclen, azimuth] = distance(Airport_data(1:2,(1:end-1)), Airport_data(1:2,(2:end)), spheroid);
    out = (distance(Airport_data(1,airport_i), Airport_data(2,airport_i),Airport_data(1,airport_j), Airport_data(2,airport_j),spheroid))/1000;
end
%}

function out = arclen(airport_i,airport_j,Airport_data)
    delta_sigma = 2*asin(sqrt((sin(deg2rad((Airport_data(1,airport_i)-Airport_data(1,airport_j))/2)))^2+cos(deg2rad(Airport_data(1,airport_i)))*cos(deg2rad(Airport_data(1,airport_j)))*(sin(deg2rad((Airport_data(2,airport_i)-Airport_data(2,airport_j))/2)))^2));
    %[arclen, azimuth] = distance(Airport_data(1:2,(1:end-1)), Airport_data(1:2,(2:end)), spheroid);
    out = 6371*delta_sigma;
end

    