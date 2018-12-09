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
        DV_Xnodes                  =   Nodes*Nodes;
        DV_Wnodes                  =   Nodes*Nodes;
        DV_Znodes                  =   Nodes*Nodes*k_ac; 

        %   Decision variables
        DV=DV_Xnodes+DV_Wnodes+DV_Znodes;
        
        
        obj                     =   ones(DV,1);
        lb                      =   zeros(DV, 1);                                 % Lower bounds
        ub                      =   inf(DV, 1);                                   % Upper bounds
        ctype                   =   char(ones(1, (DV)) * ('I'));                  % Variable types 'C'=continuous; 'I'=integer; 'B'=binary
        
        
        l = 1;                                      % Array with DV names  (OPTIONAL, BUT HELPS READING THE .lp FILE)
        for i =1:Nodes
            for j = 1:Nodes                    % of the x_{ij}^k variables
                obj(l,1)      = (5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043)*(arclen(i,j,Airport_data));
                NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '   '];
                l = l + 1;
            end
        end
        
        for i =1:Nodes
            for j = 1:Nodes                    % of the x_{ij}^k variables
                obj(l,1)      = (5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043)*(arclen(i,j,Airport_data));
                NameDV (l,:)  = ['W_' num2str(i,'%02d') ',' num2str(j,'%02d') '   '];
                l = l + 1;
            end
        end
        
        for k=1:k_ac
            for i =1:Nodes
                for j = 1:Nodes                    % of the x_{ij}^k variables
                    if g(i)==0 || g(j)==3
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
                    TAT=ACData(3,k)/60;
                    if g(j)==0
                        TAT=max([1 TAT*2]);
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
                    for m = 1:Nodes
                        C_cap(varindex(m,j,0,'w',Nodes)) = 1 - g(i);
                        C_cap(varindex(i,m,0,'w',Nodes)) = 1 - g(j);
                    end
                end
                for k = 1:k_ac
                    C_cap(varindex(i,j,k,'z',Nodes)) = -ACData(2,k)*LF;
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
        for k=1:k_ac
            for i=1:Nodes
                C_runway=zeros(1,DV);
                if ACData(5,k)>Airport_data(3,i)
                    for j=1:Nodes
                        C_runway(varindex(i,j,k,'z',Nodes))=1;
                    end
                    cplex.addRows(0, C_runway, 0, sprintf('runway%d_%d',i,k));
                end
            end
        end
        
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
                    C_flow(varindex(i,j,k,'z',Nodes))=1;
                    C_flow(varindex(j,i,k,'z',Nodes))=-1;
                end
                cplex.addRows(0, C_flow, 0, sprintf('flow%d_%d',i, k));
            end
        end
        
        fixed_cost=0
        for k=1:k_ac
            fixed_cost=fixed_cost+fleet(k)*ACData(6,k);
        end
            
     %%  Execute model
        cplex.Param.mip.limits.nodes.Cur    = 1e+8;         %max number of nodes to be visited (kind of max iterations)
        cplex.Param.timelimit.Cur           = 8*3600;         %max time in seconds
        
%   Run CPLEX
        cplex.solve();
        cplex.writeModel([model '.lp']);
    
     %%  Postprocessing
%   Store direct results
    status                      =   cplex.Solution.status;
    if status == 101 || status == 102 || status == 105  %http://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.0/ilog.odms.cplex.help/refcallablelibrary/macros/Solution_status_codes.html
        sol.profit      =   cplex.Solution.objval-fixed_cost;
        sol.values      =   cplex.Solution.x;
        for k = 1:k_ac
            sol.Flow (:,:,k)   =   round(reshape(cplex.Solution.x(varindex(1,1,k,'z', Nodes):varindex(Nodes, Nodes, k, 'z', Nodes)), Nodes, Nodes))';
        end
    end
    % Count the numbers of slots used
    slots=zeros(Nodes,1);
%   Write output
    fprintf('\n-----------------------------------------------------------------\n');
    fprintf ('Objective function value:          %10.1f  \n', sol.profit);
    fprintf ('\n') 
    fprintf ('Link From   To         AC1    AC2   AC3    Total (Demand) \n');
    NL      =   0;
    for i = 1:Nodes
        for j = 1:Nodes
            if sol.Flow(i,j,1)+sol.Flow(i,j,2)+sol.Flow(i,j,3)>0
                slots(i,1)=slots(i,1)+sol.Flow(i,j,1)+sol.Flow(i,j,2)+sol.Flow(i,j,3);
                NL      = NL + 1;
                fprintf (' %2d  %s   %s  %5d  %5d %5d %6d  (%5d) \n', NL, Airport_name{i}, ...
                            Airport_name{j}, sol.Flow (i,j,1), sol.Flow (i,j,2), ...
                            sol.Flow (i,j,3), sol.Flow (i,j,1)+sol.Flow (i,j,2)+sol.Flow(i,j,3), ...
                            Demand(i,j));
            end
        end
    end
    
    fprintf('\n------------------------Slots-------------------------------------\n');
    fprintf ('Used Available \n');
    for i=1:Nodes
        fprintf (' %2d       %5d     \n', slots(i,1), Airport_data(4,i));
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


function out = arclen(airport_i,airport_j,Airport_data)
    spheroid = wgs84Ellipsoid;
    spheroid.SemimajorAxis = spheroid.SemimajorAxis;
    spheroid.SemiminorAxis = spheroid.SemiminorAxis;
    %[arclen, azimuth] = distance(Airport_data(1:2,(1:end-1)), Airport_data(1:2,(2:end)), spheroid);
    out = (distance(Airport_data(1,airport_i), Airport_data(2,airport_i),Airport_data(1,airport_j), Airport_data(2,airport_j),spheroid))/1000;
end

    