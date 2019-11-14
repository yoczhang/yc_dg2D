function interfacemeshInfo=test_coupled_interface_1(StokesmeshInfo, DarcymeshInfo, which_coord)
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
    

interfacemeshInfo=getCoupledEdgeInfo(StokesmeshInfo, S_interfaceEdgeIndex, DarcymeshInfo, D_interfaceEdgeIndex, whichcoord);


end % function


