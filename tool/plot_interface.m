function plot_interface(StokesmeshInfo, DarcymeshInfo)

interfaceEdges_S = StokesmeshInfo.interfaceEdgeIndex;
interfaceEdges_D = DarcymeshInfo.interfaceEdgeIndex;

% for ii = 1:length(interfaceEdges_S)
%     point1 = StokesmeshInfo.node(StokesmeshInfo.edge(ii,1),:);
%     point2 = StokesmeshInfo.node(StokesmeshInfo.edge(ii,2),:);
%     hold on
%     plot(point1(1),point1(2),point2(1),point2(2),'r');
%     
% end % for ii
figure
hold on
S_point1_index = StokesmeshInfo.edge(interfaceEdges_S,1);
S_point2_index = StokesmeshInfo.edge(interfaceEdges_S,2);
S_point = union(S_point1_index,S_point2_index);
S_point_xcoord = StokesmeshInfo.node(S_point,1);
S_point_ycoord = StokesmeshInfo.node(S_point,2);
plot(S_point_xcoord,S_point_ycoord,'r-o');
for ii = 1:length(S_point)
    str_name = ['(',num2str(S_point_xcoord(ii)),',',num2str(S_point_ycoord(ii)),')'];
    str_name_1 = num2str(S_point_xcoord(ii));
    text(S_point_xcoord(ii),S_point_ycoord(ii)+S_point_ycoord(ii)/2,str_name_1');
end

% for ii = 1:length(interfaceEdges_D)
%     point1 = DarcymeshInfo.node(DarcymeshInfo.edge(ii,1),:);
%     point2 = DarcymeshInfo.node(DarcymeshInfo.edge(ii,2),:);
%     hold on
%     plot(point1(1),point1(2),point2(1),point2(2),'b');
%     
% end % for ii

D_point1_index = DarcymeshInfo.edge(interfaceEdges_D,1);
D_point2_index = DarcymeshInfo.edge(interfaceEdges_D,2);
D_point = union(D_point1_index,D_point2_index);
D_point_xcoord = DarcymeshInfo.node(D_point,1);
D_point_ycoord = DarcymeshInfo.node(D_point,2);
plot(D_point_xcoord,D_point_ycoord,'b-*');
for ii = 1:length(D_point)
    str_name = ['(',num2str(D_point_xcoord(ii)),',',num2str(D_point_ycoord(ii)),')'];
    str_name_1 = num2str(D_point_xcoord(ii));
    text(D_point_xcoord(ii),D_point_ycoord(ii)-D_point_ycoord(ii)/2,str_name_1');
    
end

end % function