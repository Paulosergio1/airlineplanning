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
    
    fleet       =   xlsread(filn,'Group 8', 'B12:F12');
    
    Population  =   xlsread(filn,'General', 'B4:C27');
    GDOP        =   xlsread(filn,'General', 'F4:G27');
 
    ACData      =   xlsread(filn,'Group 8', 'B116:F124');
    
    Nodes      =   length(Demand(1,:));

    
    %%  Initiate CPLEX model
%   Create model
        model                   =   'Airline_planning';
        cplex                   =   Cplex(model);
        cplex.Model.sense       =   'maximize';

%       Number of aircraft types
        k_ac                       =   3;
        
%       Variable g which determines if the airport is a hub or not
        g=ones(Nodes,1);
        g(3,1)=0;
        
%       Utilisation time in hours
        time_used=10;
        
        %Fuel price
        Fuel_price              = 1.42; %USD/gallon
        
%   Decision variables
    
    %%  Objective Function
        DV_Xnodes                  =   Nodes*Nodes;
        DV_Wnodes                  =   Nodes*Nodes;
        DV_Znodes                  =   Nodes*Nodes*3; 

        %   Decision variables
        DV=DV_Xnodes+DV_Wnodes+DV_Znodes;
        
        
        obj                     =   ones(DV,1);
        lb                      =   zeros(DV, 1);                                 % Lower bounds
        ub                      =   inf(DV, 1);                                   % Upper bounds
        ctype                   =   char(ones(1, (DV)) * ('I'));                  % Variable types 'C'=continuous; 'I'=integer; 'B'=binary
        
        
        l = 1;                                      % Array with DV names  (OPTIONAL, BUT HELPS READING THE .lp FILE)
        for i =1:Nodes
            for j = 1:Nodes                    % of the x_{ij}^k variables
                obj(l,1)      = 5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043;
                NameDV (l,:)  = ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '   '];
                l = l + 1;
            end
        end
        
        for i =1:Nodes
            for j = 1:Nodes                    % of the x_{ij}^k variables
                obj(l,1)      = 5.9*(arclen(i,j,Airport_data))^(-0.76)+0.043;
                NameDV (l,:)  = ['W_' num2str(i,'%02d') ',' num2str(j,'%02d') '   '];
                l = l + 1;
            end
        end
        
        for k=1:k_ac
            for i =1:Nodes
                for j = 1:Nodes                    % of the x_{ij}^k variables
                    obj(l,1)    =   -(ACData(8,k)*arclen(i,j,Airport_data)/ACData(1,k)+ACData(9,k)*Fuel_price/(1.5)*arclen(i,j,Airport_data));
                    NameDV (l,:)  = ['z_' num2str(i,'%02d') ',' num2str(j,'%02d') ',' num2str(k, '%02d') ];
                    l = l + 1;
                end
            end
        end
       
        % cplex.addCols(obj,A,lb,ub,ctype,name)  http://www-01.ibm.com/support/knowledgecenter/#!/SSSA5P_12.2.0/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX1213.html
        cplex.addCols(obj, [], lb, ub, ctype, NameDV);
        cplex.writeModel([model '.lp']);
        
        
    %%  Constraints
    %   Flow conservation at the nodes          
%         for i = 1:Nodes
%             for k = 1:Classes
%                 C1      =   zeros(1, DV);    %Setting coefficient matrix with zeros
%                 for j = 1:Nodes
%                     C1(varindex(i,j,k))   =    1;              %Link getting IN the node
%                     C1(varindex(j,i,k))   =   -1;              %Link getting OUT the node
%                 end
%                 cplex.addRows(Flow(i,k), C1, Flow(i,k), sprintf('FlowBalanceNode%d_%d',i,k));
%             end
%         end
%         
%     %   Capacity per class in each link
%         for i = 1:Nodes;
%             for j = 1:Nodes;
%                 C2      =   zeros(1, DV);       %Setting coefficient matrix with zeros
%                 for k = 1:Classes;
%                     C2(varindex(i,j,k))   =   1;      %Only the X_{i,j} (for both k) for the {i,j} pair under consideration
%                 end
%                 cplex.addRows(0, C2, Cap(i,j),sprintf('CapacityLink%d_%d_%d',i,j,k));
%             end
%         end
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
                    rightvariable=time_used*fleet(k);
                end
            end
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
        
        
     %%  Execute model
        cplex.Param.mip.limits.nodes.Cur    = 1e+8;         %max number of nodes to be visited (kind of max iterations)
        cplex.Param.timelimit.Cur           = 3600;         %max time in seconds
        
%   Run CPLEX
        cplex.solve();
        cplex.writeModel([model '.lp']);
    
     %%  Postprocessing
%   Store direct results
    status                      =   cplex.Solution.status;
    if status == 101 || status == 102 || status == 105  %http://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.0/ilog.odms.cplex.help/refcallablelibrary/macros/Solution_status_codes.html
        sol.profit      =   cplex.Solution.objval;
        for k = 1:Classes
            sol.Flow (:,:,k)   =   round(reshape(cplex.Solution.x(varindex(1,1,k):varindex(Nodes, Nodes, k)), Nodes, Nodes))';
        end
    end
%   Write output
    fprintf('\n-----------------------------------------------------------------\n');
    fprintf ('Objective function value:          %10.1f  \n', sol.profit);
    fprintf ('\n') 
    fprintf ('Link     From     To    Flow_Y   Flow_J   Total  (  Cap)    Cost \n');
    NL      =   0;
    for i = 1:Nodes
        for j = 1:Nodes
            if Cost(i,j)<10000
                NL      = NL + 1;
                if sol.Flow(i,j,1)+sol.Flow(i,j,2)>0
                    fprintf (' %2d \t  %s  \t  %s \t  %5d  %5d   %6d  (%5d)   %6d \n', NL, Airport{i}, ...
                                Airport{j}, sol.Flow (i,j,1), sol.Flow (i,j,2), ...
                                sol.Flow (i,j,1)+sol.Flow (i,j,2), Cap(i,j), ...
                                Cost(i,j)*(sol.Flow (i,j,1)+sol.Flow (i,j,2)));
                end
            end
        end
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

    