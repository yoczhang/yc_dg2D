function interfacemeshInfo=getCoupledEdgeInfo(StokesmeshInfo, S_interfaceEdgeIndex, DarcymeshInfo, D_interfaceEdgeIndex, whichcoord)
%
%   %--------------------------- Noting !!! ----------------------------------------
%       This function is designed to get the interface information of
%       StokesDarcy problem for the first time. Actually, I think this function can
%       be applied to other problems, but here we just reserve the
%       Stokes(simplified as: S), and Darcy(simplified as: D) notations,
%       because these notations can make the program easy to read.
%   %----------------------------------------------------------------------------------
%
%   input:
%       StokesmeshInfo, mesh structure-type, the Stokes mesh information.
%       S_interfaceEdgeIndex, the interface edge information belong to Stokes mesh.
%       DarcymeshInfo, mesh structure-type, the Darcy mesh information.
%       D_interfaceEdgeIndex, the interface edge information belong to Darcy mesh.
%       whichcoord, can be valued as 'xcoord', 1, or 'ycoord', 2.
%                           'xcoord' and 1 stands for  in the following sequence will be based on the x-coordinate,
%                           'ycoord' and 2 stands for  in the following sequence will be based on the y-coordinate.
%
%   output:
%       interfacemeshInfo, the coupled interface information.
%
%
%
%   YcZhang 30/8/2017
%
%   Last modified 30/8/2017
%

if strcmpi(whichcoord,'xcoord') 
    whichcoord = 1;
elseif strcmpi(whichcoord,'ycoord') 
    whichcoord = 2;
end


%% the COUPLED interface edges setting
nE_interface_S = length(S_interfaceEdgeIndex);
nE_interface_D = length(D_interfaceEdgeIndex);

% nodes and edges
node_S = StokesmeshInfo.node;
node_D = DarcymeshInfo.node;
    
edge_S = StokesmeshInfo.edge;
edge_D = DarcymeshInfo.edge;
    
interfaceNodes_S = zeros(2*nE_interface_S,3);
interfaceNodes_D = zeros(2*nE_interface_D,3);

% 
for neS = 1:nE_interface_S
    node_index = edge_S(S_interfaceEdgeIndex(neS),:);
        
    interfaceNodes_S(2*neS-1,1:2) = node_S(node_index(1),:);
    interfaceNodes_S(2*neS-1,3) = S_interfaceEdgeIndex(neS);
        
    interfaceNodes_S(2*neS,1:2) = node_S(node_index(2),:);
    interfaceNodes_S(2*neS,3) = S_interfaceEdgeIndex(neS);
end % for neS
interfaceNodes_S = sortrows(interfaceNodes_S,whichcoord); 
    %> re-sort the first column(i.e. the xcoord), or the second column(i.e. the ycoord).
    
for neD = 1:nE_interface_D
    node_index = edge_D(D_interfaceEdgeIndex(neD),:);
        
    interfaceNodes_D(2*neD-1,1:2) = node_D(node_index(1),:);
    interfaceNodes_D(2*neD-1,3) = D_interfaceEdgeIndex(neD);
        
    interfaceNodes_D(2*neD,1:2) = node_D(node_index(2),:);
    interfaceNodes_D(2*neD,3) = D_interfaceEdgeIndex(neD);
end % for neS
interfaceNodes_D = sortrows(interfaceNodes_D,whichcoord); 
    %> re-sort the first column(i.e. the xcoord), or the second column(i.e. the ycoord).
    
% couple the interface 
interfaceNodes_couple = [interfaceNodes_S; interfaceNodes_D];
interfaceNodes_couple = sortrows(interfaceNodes_couple,whichcoord);
    %> re-sort the first column(i.e. the xcoord), or the second column(i.e. the ycoord).
    
    

%% to get the interface edges belong to which element in Stokes and Darcy domain.
interfaceNodes_couple_onecoord = unique(interfaceNodes_couple(:,whichcoord));
Ninterface = length(interfaceNodes_couple_onecoord)-1; %>1, interfacemeshInfo.Ninterface

interface2StokesElem = getInterfaceInfo(whichcoord, StokesmeshInfo, interfaceNodes_couple_onecoord, interfaceNodes_S, Ninterface);
interface2DarcyElem = getInterfaceInfo(whichcoord, DarcymeshInfo, interfaceNodes_couple_onecoord, interfaceNodes_D, Ninterface);

interfaceEdge = zeros(Ninterface,2);
interfaceEdge(1:Ninterface,1) = 1:Ninterface; interfaceEdge(1:Ninterface,2) = 2:Ninterface+1; 


%% construct the interfacemeshInfor structure
interfacemeshInfo.Ninterface = Ninterface;
interfacemeshInfo.interface2StokesElem = interface2StokesElem; 
    %> i.e., [n, e] = interface2StokesElem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
    %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 
interfacemeshInfo.interface2DarcyElem = interface2DarcyElem;
    %> i.e., [n, e] = interface2DarcyElem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
    %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 
interfacemeshInfo.interfaceNode = unique(interfaceNodes_couple(:,1:2),'rows','stable');
interfacemeshInfo.interfaceEdge = interfaceEdge;

end % function



%
%--------------------------------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------------%
%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function interface2Elem = getInterfaceInfo(whichcoord, meshInfo, interfaceNodes_couple_coord, interfaceNodes, Ninterface)
%
%
%
%
%   YcZhang 29/8/2017
%

interface2Elem = zeros(Ninterface,2); 
    %> i.e., [n, e] = interface2Elem(ii,:), n standsfor the ii-the interface belongs to n-th elem, 
    %> and e standsfor the ii-th interface is the e-th edge of n-th elem. 

for ne = 1:Ninterface
    coord_node1 = interfaceNodes_couple_coord(ne);
    coord_node2 = interfaceNodes_couple_coord(ne+1);
        
    %% the Stokes 
    % node 1
    if ismember(coord_node1, interfaceNodes(:,whichcoord))
        % and notice that, if A is a matrix, c is a value, the
        % v=ismember(c, A) is a scalar, but v=ismember(A,c) is a matrix have the same size of A.
            
        node1coupledEdge = ismember(interfaceNodes(:,whichcoord),coord_node1);
        node1coupledEdge = interfaceNodes(node1coupledEdge,3);
    else
        if whichcoord == 1
            temp1 = [coord_node1, 0, 0; interfaceNodes(:,1:3)];
        elseif whichcoord == 2
            temp1 = [0, coord_node1, 0; interfaceNodes(:,1:3)];
        end
        temp1 = sortrows(temp1,whichcoord);
        temp2 = 1:size(temp1,1);
            
        node_index = temp2( ismember(temp1(:,whichcoord),coord_node1) );
        if (node_index==1) || (node_index==Ninterface+2)
            disp(['Error in node_index, node1, Ninterface = ',num2str(Ninterface)])
            disp(['ne = ',num2str(ne),'; node_index =',num2str(node_index)])
            break;
        end
        upper_node = node_index-1;
        lower_node = node_index+1;
            
        coord_upper = temp1(upper_node,1);
        coord_lower = temp1(lower_node,1);
            
        upper_edges = temp1( ismember(temp1(:,whichcoord),coord_upper), 3);
        lower_edges = temp1( ismember(temp1(:,whichcoord),coord_lower), 3);
            
        node1coupledEdge = intersect(upper_edges, lower_edges);
    end % if ismember
        
    % node 2
    if ismember(coord_node2, interfaceNodes(:,whichcoord))
        % and notice that, if A is a matrix, c is a value, the
        % v=ismember(c, A) is a scalar, but v=ismember(A,c) is a matrix have the same size of A.
            
        node2coupledEdge = ismember(interfaceNodes(:,whichcoord),coord_node2);
        node2coupledEdge = interfaceNodes(node2coupledEdge,3);
    else
        if whichcoord == 1
            temp1 = [coord_node2, 0, 0; interfaceNodes(:,1:3)];
        elseif whichcoord == 2
            temp1 = [0, coord_node2, 0; interfaceNodes(:,1:3)];
        end
        temp1 = sortrows(temp1,whichcoord);
        temp2 = 1:size(temp1,1);
        
        node_index = temp2( ismember(temp1(:,whichcoord),coord_node2) );
        if (node_index==1) || (node_index==Ninterface+2)
            disp(['Error in node_index, node2, Ninterface = ',num2str(Ninterface)])
            disp(['ne = ',num2str(ne),'; node_index =',num2str(node_index)])
            break;
        end
        upper_node = node_index-1;
        lower_node = node_index+1;
            
        coord_upper = temp1(upper_node,1);
        coord_lower = temp1(lower_node,1);
           
        upper_edges = temp1( ismember(temp1(:,whichcoord),coord_upper), 3);
        lower_edges = temp1( ismember(temp1(:,whichcoord),coord_lower), 3);
            
        node2coupledEdge = intersect(upper_edges, lower_edges);
         
    end % if ismember
      
    % get the edge (node1, node2) belong to which edge in Stokes domain.
    coupledE = intersect(node1coupledEdge,node2coupledEdge);
    
    % to get elem
    elem1 = meshInfo.edge2elem(coupledE,1);
    elem2 = meshInfo.edge2elem(coupledE,2);
    local_edge1 = meshInfo.edge2elem(coupledE,3);
    local_edge2 = meshInfo.edge2elem(coupledE,4);
    
    interface2Elem(ne,1) = max(elem1,elem2); 
        %> because the interface is also the boundary, so according edge2elem, 
        %> for the boundary edge, there has 0 between elem1 and elem2.
        
    if interface2Elem(ne,1) == elem1
        interface2Elem(ne,2) = local_edge1;
    else
        interface2Elem(ne,2) = local_edge2;
    end
        
end % for ne

end % function 