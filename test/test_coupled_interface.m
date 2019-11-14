function interfacemeshInfo=test_coupled_interface(StokesmeshInfo, DarcymeshInfo, which_coord)
%
%
%   input:
%       which_coord, 'xcoord' stands for  in the following sequence will be based on the x-coordinate,
%                           'ycoord'stands for  in the following sequence will be based on the y-coordinate.
%
%
%
%
%

% a = 0:2:12;
% a = [a'*0.1 a'*0+0.001];
% 
% b = 0:3:12;
% b = [b'*0.1 b'*0+0.001];
% 
% c = [a;b];
% 
% [coupled, i, j] = unique(c,'rows');

if strcmpi(which_coord,'xcoord') 
    whichcoord = 1;
elseif strcmpi(which_coord,'ycoord') 
    whichcoord = 2;
end


%> in pde case2, the interface is setted by y==1.
    
%% bdEdge setting
% Stokes domain bdEdge setting
S_interEdgeIndex = StokesmeshInfo.interEdgeIndex; % [Ninter x 1]
    
S_DirichletEdgeIndex = StokesmeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
S_interfaceEdgeIndex = S_DirichletEdgeIndex( StokesmeshInfo.baryEdge(S_DirichletEdgeIndex,2)==1 );
S_DirichletEdgeIndex = setdiff(S_DirichletEdgeIndex, S_interfaceEdgeIndex);
    
StokesmeshInfo.interEdgeIndex = S_interEdgeIndex;
StokesmeshInfo.DirichletEdgeIndex = S_DirichletEdgeIndex;
StokesmeshInfo.interfaceEdgeIndex = S_interfaceEdgeIndex;
    
% Darcy domain bdEdge setting
D_interEdgeIndex = DarcymeshInfo.interEdgeIndex; % [Ninter x 1]
    
D_bdEdgeIndex = DarcymeshInfo.bdEdgeIndex; % here, all the bdEdge is set as Diri edge, [Ndir x 1] 
D_interfaceEdgeIndex = D_bdEdgeIndex( DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)==1 );
D_DirichletEdgeIndex = D_bdEdgeIndex( DarcymeshInfo.baryEdge(D_bdEdgeIndex,2)==0 );
D_NeumannEdgeIndex = setdiff(D_bdEdgeIndex, union(D_interfaceEdgeIndex,D_DirichletEdgeIndex));
    
DarcymeshInfo.interEdgeIndex = D_interEdgeIndex;
DarcymeshInfo.interfaceEdgeIndex =D_interfaceEdgeIndex;
DarcymeshInfo.DirichletEdgeIndex = D_DirichletEdgeIndex;
DarcymeshInfo.NeumannEdgeIndex = D_NeumannEdgeIndex;
    

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

%
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




% function my_first_edition_just_to_Stokes_domain
% 
% for ne = 1:Ninterface
%     coord_node1 = interfaceNodes_couple_coord(ne);
%     coord_node2 = interfaceNodes_couple_coord(ne+1);
%         
%     %% the Stokes 
%     % node 1
%     if ismember(coord_node1, interfaceNodes_S(:,coord))
%         % and notice that, if A is a matrix, c is a value, the
%         % v=ismember(c, A) is a scalar, but v=ismember(A,c) is a matrix have the same size of A.
%             
%         coupledE_S1 = ismember(interfaceNodes_S(:,coord),coord_node1);
%         coupledE_S1 = interfaceNodes_S(coupledE_S1,3);
%     else
%         if coord == 1
%             temp1 = [coord_node1, 0, 0; interfaceNodes_S(:,1:3)];
%         elseif coord == 2
%             temp1 = [0, coord_node1, 0; interfaceNodes_S(:,1:3)];
%         end
%         temp1 = sortrows(temp1,coord);
%         temp2 = 1:size(temp1,1);
%             
%         node_index = temp2( ismember(temp1(:,coord),coord_node1) );
%         if (node_index==1) || (node_index==Ninterface+2)
%             disp('Error in node_index')
%             break;
%         end
%         upper_node = node_index-1;
%         lower_node = node_index+1;
%             
%         coord_upper = temp1(upper_node,1);
%         coord_lower = temp1(lower_node,1);
%             
%         upper_edges = temp1( ismember(temp1(:,coord),coord_upper), 3);
%         lower_edges = temp1( ismember(temp1(:,coord),coord_lower), 3);
%             
%         coupledE_S1 = intersect(upper_edges, lower_edges);
%     end % if ismember
%         
%     % node 2
%     if ismember(coord_node2, interfaceNodes_S(:,coord))
%         % and notice that, if A is a matrix, c is a value, the
%         % v=ismember(c, A) is a scalar, but v=ismember(A,c) is a matrix have the same size of A.
%             
%         coupledE_S2 = ismember(interfaceNodes_S(:,coord),coord_node2);
%         coupledE_S2 = interfaceNodes_S(coupledE_S2,3);
%     else
%         if coord == 1
%             temp1 = [coord_node2, 0, 0; interfaceNodes_S(:,1:3)];
%         elseif coord == 2
%             temp1 = [0, coord_node2, 0; interfaceNodes_S(:,1:3)];
%         end
%         temp1 = sortrows(temp1,coord);
%         temp2 = 1:size(temp1,1);
%         
%         node_index = temp2( ismember(temp1(:,coord),coord_node2) );
%         if (node_index==1) || (node_index==Ninterface+2)
%             disp('Error in node_index')
%             break;
%         end
%         upper_node = node_index-1;
%         lower_node = node_index+1;
%             
%         coord_upper = temp1(upper_node,1);
%         coord_lower = temp1(lower_node,1);
%            
%         upper_edges = temp1( ismember(temp1(:,coord),coord_upper), 3);
%         lower_edges = temp1( ismember(temp1(:,coord),coord_lower), 3);
%             
%         coupledE_S2 = intersect(upper_edges, lower_edges);
%          
%     end % if ismember
%       
%     % get the edge (node1, node2) belong to which edge in Stokes domain.
%     interface2StokesElem(ne) = intersect(coupledE_S1,coupledE_S2);
% 
% end % for ne
% 
% end % function
