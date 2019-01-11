%Group 8:
%Koen den Hertog - 4457803
%Paul van Kessel - 4453182
%Benyamin de Leeuw - 4471121


itterations=5;
tot_passengers=zeros(itterations,2);
for iter=1:itterations
    close all
    demand ()
    close all
    tot_passengers=Airlineplanning(iter,tot_passengers);
end
hold off
figure(20)
bar([1 2 3 4 5],tot_passengers)




function tot_passengers = Airlineplanning (itterations,tot_passengers)
%%  Initialization
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
    savepath
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\examples\src\matlab');
    savepath
    close all
    warning('off','MATLAB:lang:badlyScopedReturnValue')
    warning('off','MATLAB:xlswrite:NoCOMServer')

    %%  Determine input
%   Select input file and sheet
    filn        =   [pwd '/AE4423_Datasheets.xlsx'];
    filn2       =   [pwd '/Group8_results.xlsx'];
    
    Demand(:,:,1)      =   xlsread(filn2,'new_demands', 'A1:X24'); % High season
    Demand(:,:,2)      =   xlsread(filn2,'new_demands', 'A25:X48'); % Low season

    Airport_data=   xlsread(filn,'Group 8', 'C6:Z9');
    [~,Airport_name] =   xlsread(filn,'Group 8', 'C5:Z5'); 
    
    fleet       =   xlsread(filn,'Group 8', 'B12:F12');

 
    ACData      =   xlsread(filn,'Group 8', 'B116:F124');
    
    Nodes      =    length(Demand(1,:,1));
    
    %%  Initiate CPLEX model
%   Create model
        model                   =   'Airline_planning';
        cplex                   =   Cplex(model);
        cplex.Model.sense       =   'maximize';

%       Number of aircraft types
        k_ac                       =   5;
        
%       Number of seasons
        season=2;
        
%       Load factor
        LF = 0.75;
        LF_US = 0.85;
        
%       Variable g which determines if the airport is a hub or not
        g=ones(Nodes,1);
        g(3,1)=0; % G is zero for the hub
        
%       Utilisation time in hours
        time_used=10*7;
        
        %Fuel price
        Fuel_price              = 1.6; %USD/gallon
        Max_US_flow             = 7500;
        
%   Decision variables
    
    %%  Objective Function
        DV_Xnodes                  =   Nodes*Nodes*season; % number of nodes for direct passenger
        DV_Wnodes                  =   Nodes*Nodes*season; % number of nodes for transfer passenger
        DV_Znodes                  =   Nodes*Nodes*k_ac*season;% number of nodes for flights with ac type
        DV_Knodes                  =   k_ac*2; %number of nodes for leasing additional ac or stopping contracts

        %   Decision variables
        DV=DV_Xnodes+DV_Wnodes+DV_Znodes+DV_Knodes;
        
        
        obj                     =   ones(DV,1);
        lb                      =   zeros(DV, 1);                                 % Lower bounds
        ub                      =   inf(DV, 1);                                   % Upper bounds
        ctype                   =   char(ones(1, (DV)) * ('I'));                  % Variable types 'C'=continuous; 'I'=integer; 'B'=binary
        
        
        l = 1;                                      % Array with DV names  (OPTIONAL, BUT HELPS READING THE .lp FILE)
        for p=1:season %For both high and low season 
            for i =1:Nodes % objective function values for direct passengers
                for j = 1:Nodes                    % of the x_{ij}^k variables
                    if i>20 || j>20
                        obj(l,1)      = 0.05*(arclen(i,j,Airport_data));
                    else
                        obj(l,1)      = (5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043)*(arclen(i,j,Airport_data));
                    end
                    NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') ',' num2str(p,'%02d') '   '];
                    l = l + 1;
                end
            end
        end
        
        for p=1:season %For both high and low season 
            for i =1:Nodes % objective function values for trasfer passenger
                for j = 1:Nodes                    % of the x_{ij}^k variables
                    if i>20 || j>20
                        obj(l,1)      = 0.05*(arclen(i,j,Airport_data));
                    else
                        obj(l,1)      = (5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043)*(arclen(i,j,Airport_data));
                    end
                    NameDV (l,:)  = ['W_' num2str(i,'%02d') ',' num2str(j,'%02d') ',' num2str(p,'%02d') '   '];
                    l = l + 1;
                end
            end
        end
        
        for p=1:season %For both high and low season  
            for k=1:k_ac % Objective function values for flights per ac type
                for i =1:Nodes
                    for j = 1:Nodes                    % of the x_{ij}^k variables
                        if g(i)==0 || g(j)==0 %Economies of scale, for flights arriving or departing at hub
                           obj(l,1)    =   -0.7*(ACData(8,k)*arclen(i,j,Airport_data)/ACData(1,k)+ACData(9,k)*Fuel_price/(1.5)*arclen(i,j,Airport_data)+ACData(7,k));
                        else
                            obj(l,1)    =   -(ACData(8,k)*arclen(i,j,Airport_data)/ACData(1,k)+ACData(9,k)*Fuel_price/(1.5)*arclen(i,j,Airport_data)+ACData(7,k));
                        end
                        NameDV (l,:)  = ['z_' num2str(i,'%02d') ',' num2str(j,'%02d') ',' num2str(k, '%02d') ',' num2str(p,'%02d') ];
                        l = l + 1;
                    end
                end
            end
        end
        
        for k=1:k_ac
            obj(l,1) = 2*(-8000+ACData(6,k));
            NameDV(l,:) = ['K_stop__' num2str(k,'%02d') '   '];
            l=l+1;
        end
        
        for k=1:k_ac
            obj(l,1) = 2*(-2000-ACData(6,k));
            NameDV(l,:) = ['K_extra_' num2str(k,'%02d') '   '];
            l=l+1;
        end
            
       
        % cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
        cplex.addCols(obj, [], lb, ub, ctype, NameDV);
        cplex.writeModel([model '.lp']);
        
    %%  Fixed cost calculation
        fixed_cost=0;
        for k=1:k_ac
            fixed_cost=fixed_cost+2*fleet(k)*ACData(6,k);
        end
        
    %%  Constraints
        % Constraints need to be evaluated for both the high and low demand
        % season
        utilisation_time=zeros(Nodes,Nodes,k_ac,season);
        for p=1:season
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
                        C_time_ac(varindex(i,j,k,p,'z',Nodes))=distance/ACData(1,k)+TAT;
                        utilisation_time(i,j,k,p)=distance/ACData(1,k)+TAT;
                    end
                end
                C_time_ac(varindex(1,1,k,p,'s',Nodes))=time_used;
                C_time_ac(varindex(1,1,k,p,'e',Nodes))=-time_used;
                rightvariable=time_used*fleet(k);
                cplex.addRows(-inf, C_time_ac, rightvariable, sprintf('Timeusedac%d_%d',k,p));
            end
        %   Passengers not more than demand
            for i = 1:Nodes
                for j = 1:Nodes
                    C_dem = zeros(1,DV);
                    C_dem(varindex(i,j,0,p,'x',Nodes)) = 1;
                    C_dem(varindex(i,j,0,p,'w',Nodes)) = 1;
                    cplex.addRows(0, C_dem, Demand(i,j,p), sprintf('DemandConstraint%d_%d_%d',i,j,p));
                end
            end
        %   No transfer if one the airports is hub 
            for i = 1:Nodes
                for j = 1:Nodes
                    C_hub     =   zeros(1, DV);
                    C_hub(varindex(i,j,0,p,'w',Nodes)) = 1;
                    hub_or_not = g(i) * g(j);
                    cplex.addRows(0, C_hub, Demand(i,j,p) * hub_or_not,sprintf('TransferHub%d_%d_%d',i,j,p));
                end
            end

        %   Max aircraft capacity
            for i = 1:Nodes
                for j = 1:Nodes
                    C_cap = zeros(1,DV);
                    C_cap(varindex(i,j,0,p,'x',Nodes)) = 1;
                    if g(j) == 0 || g(i) == 0
                        if i~=j
                            for m = 1:Nodes
                                C_cap(varindex(m,j,0,p,'w',Nodes)) = 1 - g(i);
                                C_cap(varindex(i,m,0,p,'w',Nodes)) = 1 - g(j);
                            end
                        end
                    end
                    for k = 1:k_ac
                        if i>20 || j>20
                            C_cap(varindex(i,j,k,p,'z',Nodes)) = -ACData(2,k)*LF_US;
                        else
                            C_cap(varindex(i,j,k,p,'z',Nodes)) = -ACData(2,k)*LF;
                        end
                    end
                    cplex.addRows(-inf, C_cap,0,sprintf('ACcapacity%d_%d_%d',i,j,p));
                end
            end


            %Only aircraft with a long enough range can fly between two cities
            for k=1:k_ac
                for i=1:Nodes
                    for j=1:Nodes
                        C_range=zeros(1,DV);
                        distance=arclen(i,j,Airport_data);
                        C_range(varindex(i,j,k,p,'z',Nodes))=1;
                        if distance>ACData(4,k)
                            cplex.addRows(0, C_range, 0, sprintf('range%d_%d_%d_%d',i,j,k,p));
                        end
                    end
                end
            end

           %Constraint for the runway length wich should be long enough. 
            for k=1:k_ac
                for i=1:Nodes
                    C_runway=zeros(1,DV);
                    if ACData(5,k)>Airport_data(3,i)
                        for j=1:Nodes
                            C_runway(varindex(i,j,k,p,'z',Nodes))=1;
                        end
                        cplex.addRows(0, C_runway, 0, sprintf('runway%d_%d_%d',i,k,p));
                    end
                end
            end

            % Slots contraint
            for j=1:Nodes
                C_slots=zeros(1,DV);
                for i=1:Nodes
                    for k=1:k_ac
                       C_slots(varindex(i,j,k,p,'z',Nodes))=1;
                    end
                end
                cplex.addRows(0, C_slots, Airport_data(4,j), sprintf('slots%d_%d',j,p));
            end

            %flow inside of airport should be equal to flow outside of airport
            for k=1:k_ac
                for i=1:Nodes
                    C_flow=zeros(1,DV);
                    for j=1:Nodes
                        if i~=j
                            C_flow(varindex(i,j,k,p,'z',Nodes))=1;
                            C_flow(varindex(j,i,k,p,'z',Nodes))=-1;
                        end
                    end
                    cplex.addRows(0, C_flow, 0, sprintf('flow%d_%d_%d',i,k,p));
                end
            end

            %Max amount of passengers to the US
            C_US_capacity=zeros(1,DV);
            for i=1:Nodes
                for j=21:Nodes %Only flights with destination to the US
                    C_US_capacity(varindex(i,j,1,p,'x',Nodes))=1; %direct passengers from hub
                    C_US_capacity(varindex(i,j,1,p,'w',Nodes))=1; %Transfer passengers from europe
                end
            end
            cplex.addRows(0, C_US_capacity,Max_US_flow, sprintf('Max_capacity_EU_to_US%d',p)); % Total passenger should be lower than maximum allowed

            %Max amounr of passengers from the US to the EU
            C_EU_capacity=zeros(1,DV);
            for i=21:Nodes
                for j=1:Nodes %Only flights with destination to the US
                    C_EU_capacity(varindex(i,j,1,p,'x',Nodes))=1; %direct passengers from hub
                    C_EU_capacity(varindex(i,j,1,p,'w',Nodes))=1; %Transfer passengers from europe
                end
            end
            cplex.addRows(0, C_EU_capacity,Max_US_flow, sprintf('Max_capacity_US_to_EU%d',p)); % Total passenger should be lower than maximum allowed


            % No direct flights between european cities and US, besides hub,
            % and no inner US flights
            C_no_hub_US=zeros(1,DV);
            for i=1:Nodes %Only consider flights departing from european flights
                for j=21:Nodes %Only codering flights arriving in the US
                    if g(i)==1 % IF the departing city is not the hub
                        C_no_hub_US(varindex(i,j,1,p,'x',Nodes))=1; % No flight from europe city to US
                        C_no_hub_US(varindex(j,i,1,p,'x',Nodes))=1; % No flight from US to europe city
                    end
                end
            end
            cplex.addRows(0, C_no_hub_US, 0, sprintf('No_hub_US%d',p));

            % No flights within europe for ac type 4 and 5
            C_ac4_ac5=zeros(1,DV);
            for k=4:k_ac
                for i=1:Nodes-4
                    for j=1:Nodes-4
                        C_ac4_ac5(varindex(i,j,k,p,'z',Nodes))=1;
                    end
                end
            end
            cplex.addRows(0,C_ac4_ac5, 0, sprintf('No_europe_flight_ac4_ac5_%d',p));
        
        end
        
        %Do not stop more contracts than there are ac available, for ac
        %type
        for k=1:k_ac
            C_stop_contract=zeros(1,DV);
            C_stop_contract(varindex(1,1,k,p,'s',Nodes))=1;
            cplex.addRows(0, C_stop_contract, fleet(k), sprintf('Max_stop_contract%d',k));
        end
            
        
     %%  Execute model
        %cplex.Param.mip.limits.nodes.Cur    = 1e+5;         %max number of nodes to be visited (kind of max iterations)
        cplex.Param.timelimit.Cur           = 500;         %max time in seconds
        cplex.Param.mip.tolerances.mipgap.Cur   = 0.009;
        
        
%   Run CPLEX
        cplex.writeModel([model '.lp']);
        cplex.solve();
    
     %%  Postprocessing
%   Store direct results
    status                      =   cplex.Solution.status;
    if status == 101 || status == 102 || status == 105  %http://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.0/ilog.odms.cplex.help/refcallablelibrary/macros/Solution_status_codes.html
        sol.profit      =   cplex.Solution.objval-fixed_cost;
        sol.values      =   cplex.Solution.x;
        for p = 1:season
            for k = 1:k_ac
                sol.Flow (:,:,k+(p-1)*k_ac)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,p,'z', Nodes):varindex(Nodes, Nodes, k, p, 'z', Nodes)), Nodes, Nodes))';
            end
        end
    end
    % Count the numbers of slots used
    slots=zeros(Nodes,2);
%   Write output
    new_fleet=zeros(k_ac,1);
    % The new fleet after considering buying additional and selling
    for k=1:k_ac
            new_fleet(k,1)=fleet(k)-cplex.Solution.x(varindex(1,1,k,1,'s',Nodes),1)+cplex.Solution.x(varindex(1,1,k,1,'e',Nodes),1);
    end
    
    %Profit high season
    profit_high=-fixed_cost*0.5;
    for i=varindex(1,1,1,1,'x',Nodes):varindex(Nodes,Nodes,k_ac,1,'x',Nodes)
        if isnan(obj(i,1))==0
            profit_high=profit_high+obj(i,1)*sol.values(i);
        end
    end
    for i=varindex(1,1,1,1,'w',Nodes):varindex(Nodes,Nodes,k_ac,1,'w',Nodes)
        if isnan(obj(i,1))==0
            profit_high=profit_high+obj(i,1)*sol.values(i);
        end
    end
    for i=varindex(1,1,1,1,'z',Nodes):varindex(Nodes,Nodes,k_ac,1,'z',Nodes)
        if isnan(obj(i,1))==0
            profit_high=profit_high+obj(i,1)*sol.values(i);
        end
    end
   
   
    %Profit low season
    profit_low=-fixed_cost*0.5;
    for i=varindex(1,1,1,2,'x',Nodes):varindex(Nodes,Nodes,k_ac,2,'x',Nodes)
        if isnan(obj(i,1))==0
            profit_low=profit_low+obj(i,1)*sol.values(i);
        end
    end
    for i=varindex(1,1,1,2,'w',Nodes):varindex(Nodes,Nodes,k_ac,2,'w',Nodes)
        if isnan(obj(i,1))==0
            profit_low=profit_low+obj(i,1)*sol.values(i);
        end
    end
    for i=varindex(1,1,1,2,'z',Nodes):varindex(Nodes,Nodes,k_ac,2,'z',Nodes)
        if isnan(obj(i,1))==0
            profit_low=profit_low+obj(i,1)*sol.values(i);
        end
    end
    
    % Adding cost for stopping/new constract
    for i=1:k_ac
        profit_high=profit_high+0.5*obj(varindex(1,1,i,2,'e',Nodes),1)*sol.values(varindex(1,1,i,2,'e',Nodes))+...
        0.5*obj(varindex(1,1,i,2,'s',Nodes),1)*sol.values(varindex(1,1,i,2,'s',Nodes));
        profit_low=profit_low+0.5*obj(varindex(1,1,i,2,'e',Nodes),1)*sol.values(varindex(1,1,i,2,'e',Nodes))+...
        0.5*obj(varindex(1,1,i,2,'s',Nodes),1)*sol.values(varindex(1,1,i,2,'s',Nodes));
    end
    
    fprintf('\n-------------------------Original fleet---------------------------\n');
    fprintf('AC type 1: (1) %d \n',new_fleet(1,1));
    fprintf('AC type 2: (1) %d \n',new_fleet(2,1));
    fprintf('AC type 3: (1) %d \n',new_fleet(3,1));
    fprintf('AC type 4: (0) %d \n',new_fleet(4,1));
    fprintf('AC type 5: (0) %d \n',new_fleet(5,1));
    fprintf('\n-------------------------Original fleet---------------------------\n');
    fprintf ('Objective function value:          %10.1f  \n', sol.profit);
    utilisation=zeros(1,k_ac,season);
    for p=1:season
        sol.PassengerDirect (:,:)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,p,'x', Nodes):varindex(Nodes, Nodes, k,p, 'x', Nodes)), Nodes, Nodes))';
        sol.PassengerIndirect (:,:)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,p,'w', Nodes):varindex(Nodes, Nodes, k,p, 'w', Nodes)), Nodes, Nodes))';
        tot_passengers(itterations,p)=sum(sum(sol.PassengerDirect (:,:)))+sum(sum(sol.PassengerIndirect (:,:)));
        figure(p);
            worldmap([20 45],[-120 -70]) % For the US
%         worldmap([35 65],[-15 30]) % For the EU
        land = shaperead('landareas.shp', 'UseGeoCoords', true);
        geoshow(land, 'FaceColor', [0.6 0.6 0.6])
        color = ['r','g','b','c','m'];
        hold on
        if p==1
            fprintf('\n-------------------Network operated high season------------------\n');
            fprintf ('Profit:          %10.1f  \n', profit_high);
        else 
            fprintf('\n-------------------Network operated low season------------------\n');
            fprintf ('Profit:          %10.1f  \n', profit_low);
        end
        fprintf ('\n') 
        fprintf ('Link From   To         AC1    AC2   AC3   AC4   AC5    Total(indirect) (Revenue per seat)  (Revenue) (ASK)  (RASK) (CASK) (RPK) (Profit) (Yield) (ALF) (BELF)\n');
        NL      =   0;
        NL_EU      =   0;
        NL_US      =   0;
        alf_array_EU = [];
        alf_array_US = [];
        belf_array_EU = [];
        belf_array_US = [];
        profit_array = [];
        for i = 1:Nodes
            for j = 1:Nodes
                if sum(sol.Flow(i,j,1+(p-1)*k_ac:p*k_ac))>0
                    slots(i,p)=slots(i,p)+sum(sol.Flow(i,j,1+(p-1)*k_ac:p*k_ac));
                    NL      = NL + 1;
                    yield=obj(varindex(i,j,1,p,'x',Nodes))*sol.PassengerDirect(i,j);
                    pass_transfer = 0;
                    if g(i) == 0 && g(j) ~= 0
                        %sol.PassengerIndirect(i,j)
                        for l = 1:Nodes
                            %for m = 1:Nodes
                                if l == 3 || l == j
                                    continue;
                                end
                                yield = yield + obj(varindex(l,j,1,p,'w',Nodes))*sol.PassengerIndirect(l,j)*(arclen(3,j,Airport_data)/(arclen(3,j,Airport_data)+arclen(l,3,Airport_data)));
                                pass_transfer = pass_transfer + sol.PassengerIndirect(l,j);
                            %end
                        end
                    elseif g(j) == 0 && g(i) ~= 0
                        %sol.PassengerIndirect(i,j)
                        %for l = 1:Nodes
                            for m = 1:Nodes
                                if m == 3 || m == i
                                    continue;
                                end
                                yield = yield + obj(varindex(i,m,1,p,'w',Nodes))*sol.PassengerIndirect(i,m)*(arclen(i,3,Airport_data)/(arclen(3,m,Airport_data)+arclen(i,3,Airport_data)));
                                pass_transfer = pass_transfer + sol.PassengerIndirect(i,m);
                            end
                        %end
                    end
                    profit = yield;
                    for k= 1:k_ac
                        profit=profit+obj(varindex(i,j,k,p,'z',Nodes))*sol.Flow(i,j,k+(p-1)*k_ac);
                    end
                    profit_array = [profit_array,profit];
                    for k= 1:k_ac
                        utilisation(1,k,p)=utilisation(1,k,p)+sol.Flow(i,j,k+(p-1)*k_ac)*utilisation_time(i,j,k,p);
                    end
                    revenue = obj(varindex(i,j,1,p,'x',Nodes))*(sol.Flow (i,j,1+(p-1)*k_ac)*ACData(2,1) + sol.Flow (i,j,2+(p-1)*k_ac)*ACData(2,2) + sol.Flow (i,j,3+(p-1)*k_ac)*ACData(2,3) + sol.Flow (i,j,4+(p-1)*k_ac)*ACData(2,4) + sol.Flow (i,j,5+(p-1)*k_ac)*ACData(2,5));
                    ask = arclen(i,j,Airport_data)*(sol.Flow (i,j,1+(p-1)*k_ac)*ACData(2,1) + sol.Flow (i,j,2+(p-1)*k_ac)*ACData(2,2) + sol.Flow (i,j,3+(p-1)*k_ac)*ACData(2,3) + sol.Flow (i,j,4+(p-1)*k_ac)*ACData(2,4) + sol.Flow (i,j,5+(p-1)*k_ac)*ACData(2,5));
                    rask = revenue/ask;
                    cost_tot = revenue-profit;
                    cask = cost_tot/ask;
                    rpk = arclen(i,j,Airport_data)*(sol.PassengerDirect(i,j)+pass_transfer);
                    alf = (sol.PassengerDirect(i,j)+pass_transfer)/(sol.Flow (i,j,1+(p-1)*k_ac)*ACData(2,1) + sol.Flow (i,j,2+(p-1)*k_ac)*ACData(2,2) + sol.Flow (i,j,3+(p-1)*k_ac)*ACData(2,3) + sol.Flow (i,j,4+(p-1)*k_ac)*ACData(2,4) + sol.Flow (i,j,5+(p-1)*k_ac)*ACData(2,5));
                    belf = (cost_tot/revenue)*alf;
                    if i > 20 || j > 20
                        alf_array_US = [alf_array_US,alf];
                        belf_array_US = [belf_array_US,belf];
                        NL_US      = NL_US + 1;
                    else
                        alf_array_EU = [alf_array_EU,alf];
                        belf_array_EU = [belf_array_EU,belf];
                        NL_EU      = NL_EU + 1;
                    end
                    indirect = pass_transfer /(pass_transfer + sol.PassengerDirect(i,j));
                    fprintf (' %2d  %s   %s  %5d  %5d %5d %5d %5d  %6d %5d %5d  %5d 5%d %5d %5d %5d %5d %5d %5d %5d\n', NL, Airport_name{i}, ...
                                Airport_name{j}, sol.Flow (i,j,1+(p-1)*k_ac), sol.Flow (i,j,2+(p-1)*k_ac), ...
                                sol.Flow (i,j,3+(p-1)*k_ac), sol.Flow (i,j,4+(p-1)*k_ac), sol.Flow(i,j,5+(p-1)*k_ac), ...
                                sol.Flow (i,j,1+(p-1)*k_ac)+sol.Flow (i,j,2+(p-1)*k_ac)+sol.Flow(i,j,3+(p-1)*k_ac)+...
                                sol.Flow (i,j,4+(p-1)*k_ac)+sol.Flow(i,j,5+(p-1)*k_ac), indirect, ...
                                obj(varindex(i,j,1,p,'x',Nodes)), revenue, ask, rask, cask, rpk, profit, yield, alf, belf);
                    for k=1:k_ac
%                         if sol.Flow(i,j,(p-1)*k_ac+k)>0 && i<=20 && j<=20  % For the EU
                        if sol.Flow(i,j,(p-1)*k_ac+k)>0   % For the US    
                            h = geoshow([Airport_data(1,i);Airport_data(1,j)],...
                                    [Airport_data(2,i);Airport_data(2,j)]);
                            h.Marker = '*';
                            h.Color = color(k);
                            h.LineWidth = sol.Flow(i,j,(p-1)*k_ac+k);
                        end
                
%                 x_loc=(Airport_data(1,i)+Airport_data(1,j))/2;
%                 y_loc=(Airport_data(2,i)+Airport_data(2,j))/2;
%                 Capacity=ACData(2,1:k_ac)*[sol.values(varindex(i,j,1,'z', Nodes)),sol.values(varindex(i,j,2,'z', Nodes)),sol.values(varindex(i,j,3,'z', Nodes))]';
%                 textm(x_loc,y_loc,['(',num2str(Capacity), ')'])
                    textm(Airport_data(1,i), Airport_data(2,i),Airport_name(i))
                    end
                end
            end
        end
        figure(15+p);
        worldmap([35 65],[-15 30]) % For the EU

        land = shaperead('landareas.shp', 'UseGeoCoords', true);
        geoshow(land, 'FaceColor', [0.6 0.6 0.6])
        color = ['r','b','g','c','m'];
        hold on
        list=[10 1 4 24];
        for loop=1:size(list,2)
            i=list(loop);
            for j=1:Nodes
                if sol.PassengerIndirect(i,j)>0 || sol.PassengerDirect(i,j)>0
                    h = geoshow([Airport_data(1,i);Airport_data(1,j)],...
                        [Airport_data(2,i);Airport_data(2,j)]);
                    h.Marker = '*';
                    h.Color = color(loop);
                    h.LineWidth = ceil((sol.PassengerIndirect(i,j)+sol.PassengerDirect(i,j))/40);
                end
            end  
        end
        hold off
    end
    
    
    
    extra_fixed_cost=0;
    for k=1:k_ac
        extra_fixed_cost=extra_fixed_cost-sol.values(varindex(1,1,k,p,'s', Nodes))*obj(varindex(1,1,k,p,'s', Nodes),1)+-sol.values(varindex(1,1,k,p,'e', Nodes))*obj(varindex(1,1,k,p,'e', Nodes),1);
    end
    
    anlf_EU = sum(alf_array_EU)/NL_EU
    nbelf_EU = sum(belf_array_EU)/NL_EU
    anlf_US = sum(alf_array_US)/NL_US
    nbelf_US = sum(belf_array_US)/NL_US
    tot_profit = sum(profit_array)- 0.5*fixed_cost - 0.5*extra_fixed_cost
    
    fprintf('\n------------------------High Season Slots-------------------------------------\n');
    fprintf ('Used Available \n');
    for i=1:Nodes
        fprintf (' %2d       %5d     \n', slots(i,1), Airport_data(4,i));
    end
   
    fprintf('\n------------------------Low Season Slots-------------------------------------\n');
    fprintf ('Used Available \n');
    for i=1:Nodes
        fprintf (' %2d       %5d     \n', slots(i,2), Airport_data(4,i));
    end
    
        fprintf('\n------------------------Utilisation High Season-------------------------------------\n');
    fprintf ('                 Used Available \n');
    for i=1:k_ac
        fprintf ('Aircraft type %1d %2d       %5d     \n',i, utilisation(1,i,1), time_used*new_fleet(i,1));
    end
    
        fprintf('\n------------------------Utilisation Low Season-------------------------------------\n');
    fprintf ('                 Used Available \n');
    for i=1:k_ac
        fprintf ('Aircraft type %1d %2d       %5d     \n',i, utilisation(1,i,2), time_used*new_fleet(i,1));
    end
    
    
    fprintf('\n------------------------Cost per AC High Season-------------------------------------\n');
    for k=1:k_ac
        cost=0;
        for i=1:Nodes
            for j=1:Nodes
                cost=cost+obj(varindex(i,j,k,1,'z', Nodes))*sol.values(varindex(i,j,k,1,'z', Nodes));
            end
        end
        fprintf ('Cost AC type %d:  %10.1f  \n',k, -cost);
    end
    
        fprintf('\n------------------------Cost per AC Low Season-------------------------------------\n');
    for k=1:k_ac
        cost=0;
        for i=1:Nodes
            for j=1:Nodes
                cost=cost+obj(varindex(i,j,k,2,'z', Nodes))*sol.values(varindex(i,j,k,2,'z', Nodes));
            end
        end
        fprintf ('Cost AC type %d:  %10.1f  \n',k, -cost);
    end
    %% Write frequency to excel file in the group8-data tab

    for i=1:Nodes
        for j=1:Nodes
            sol.Flow(i,j,season*k_ac+1)=sum(sol.Flow(i,j,1:k_ac));
        end
    end
    for i=1:Nodes
        for j=1:Nodes
            sol.Flow(i,j,season*k_ac+2)=sum(sol.Flow(i,j,k_ac+1:season*k_ac));
        end
    end
    
    xlswrite(filn2,sol.Flow(:,:,season*k_ac+1),'Group8-data','A1:X24')
    xlswrite(filn2,sol.Flow(:,:,season*k_ac+2),'Group8-data','A25:X48')
end


function out = varindex(i,j,k,season,letter,nodes)
    if letter == 'x'
        if season==1
            out=(i-1)*nodes+j;
        else
            out=(i-1)*nodes+j+nodes^2;
        end
    elseif letter == 'w'
        if season==1
            out=(i-1)*nodes+j+2*nodes^2;
        else
            out=(i-1)*nodes+j+3*nodes^2;
        end
    elseif letter == 'z'
        if season==1
            out=(i-1)*nodes+j+4*nodes^2+(k-1)*nodes^2;
        else
            out=(i-1)*nodes+j+9*nodes^2+(k-1)*nodes^2;
        end
    elseif letter == 's'% gives index for K-stop variable
        out=4*nodes^2+10*nodes^2+k;
    elseif letter == 'e'% gives index for K-extra variable
        out=4*nodes^2+10*nodes^2+k+5;
    end   
%         Function given the variable index for each DV (i,j,k) the letter
%         denotes wheter you would like to have the variable x,w,z or k. 

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

function demand ()
    %%  Determine inputsni
%   Select input file and sheet
    filn        =   [pwd '/AE4423_Datasheets.xlsx'];
    filn2       =   [pwd '/Group8_results.xlsx'];
    
    frequencies(:,:,1) = xlsread(filn2,'Group8-data','A1:X24');
    frequencies(:,:,2) = xlsread(filn2,'Group8-data','A25:X48');
    frequencies_c = xlsread(filn,'Group 8','C89:Z112');
    
    demand_low = xlsread(filn,'Group 8','C63:Z86');
    demand_high = xlsread(filn,'Group 8','C37:Z60');
    
    a = 1.0;
    b = 1.7;
    
    %% Market Share
    for t=1:2
        demand_low_new = zeros(length(frequencies(:,:,t)),length(frequencies(:,:,t)'));
        demand_high_new = zeros(length(frequencies(:,:,t)),length(frequencies(:,:,t)'));
        for i = 1:length(frequencies(:,:,t))
            for j = 1:length(frequencies(:,:,t)')
                freq_d = frequencies(i,j,t);
                freq_i = min(frequencies(3,j,t),frequencies(i,3,t));
                freq_c = frequencies_c(i,j);
                ms = (freq_d^a + freq_i^b)/(freq_d^a + freq_i^b + freq_c^a + 1e-15);
                d_lo = ms*demand_low(i,j);
                d_hi = ms*demand_high(i,j);
                demand_low_new(i,j) = d_lo;
                demand_high_new(i,j) = d_hi;
            end
        end
        if t==1
            demand_new = round(demand_high_new);
            xlswrite(filn2,demand_new,'new_demands','A1:X24')
        elseif t==2
            demand_new = round(demand_low_new);
            xlswrite(filn2,demand_new,'new_demands','A25:X48')
        end
    end
end

    