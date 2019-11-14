function plotUhAlong_fixedY_coord(Uh, y_coord, meshInfo, degreek, varargin)
%
%   Plot the Uh on the given fixed-Y-coordinate line.
%
%   input:
%       Uh, the numerical solution.
%       y_coord, the fixed-Y-coordinate.
%       meshInfo, the mesh information.
%       degreek, the degree of polynominal.
%
%
%   YcZhang 2/10/2017
%
%   Last modified 2/10/2017
%
%

if length(varargin)==1
    using_local_bases = 0;
else
    using_local_bases = 1;
end 


Nbases = nchoosek(degreek+2,2);
node = meshInfo.node;
elem = meshInfo.elem;

xx = node(:,1); yy = node(:,2);

%% level-func
fixed_y_coord = y_coord;
ff = yy - fixed_y_coord;
levelFunc.ff = ff;
[cutElem,cutNode] = getEnrichedElem(levelFunc, node, elem, 0);

% if there have some ff == 0, stands for: 
%   Level-set function: ff, is exactly zero at some nodes,
%   then we need to choose one side of the elems.
if isempty(find(sign(ff)==0,1)) == 0
    barycutElem_y = meshInfo.baryElem(cutElem,2);
    cutElem = cutElem(barycutElem_y <= y_coord);
end 

hold on
Xpoints = [];
Uhpoints = [];
for CurrElem = 1:length(cutElem)
    CurrCutElem = cutElem(CurrElem);
    
    uh_local = Uh((CurrCutElem-1)*Nbases+1:CurrCutElem*Nbases);

    CurrCutNode = cutNode(cutNode(:,3)==CurrCutElem,1:2);
    
    x_min = min(CurrCutNode(:,1));
    x_max = max(CurrCutNode(:,1));
    
    x_points = (x_min:(x_max-x_min)/1000:x_max)';
    y_points = fixed_y_coord*ones(length(x_points),1);
    
    if using_local_bases
        elem_xT = meshInfo.centroidElem(CurrCutElem,1);
        elem_yT = meshInfo.centroidElem(CurrCutElem,2);
        elem_hT = meshInfo.hElem(CurrCutElem);
    else
        elem_xT = 0;
        elem_yT = 0;
        elem_hT = 1;
    end 
    
    [Pb, ~, ~] = ...
        localBases2D(elem_xT, elem_yT, elem_hT, x_points, y_points, 1);
        %> elePtrialPb, elePtrialPbx, elePtrialPby, [Npoints x NTg1bases_trial]
        
    if Pb ~= 0
        uh_points = Pb*uh_local;
    else
        continue
    end
    
%     uh_points = zeros(length(x_points),1);
%     for ii = 1:length(x_points)
%         uh_points(ii) = fe_solution(x_points(ii),y_points(ii),uh_local,0,0,1);
%     end % for
    
    Xpoints_temp = [Xpoints; x_points];
    Uhpoints_temp = [Uhpoints; uh_points];
    
    Xpoints = Xpoints_temp; Uhpoints = Uhpoints_temp;
end % for

X_Uh = [Xpoints, Uhpoints];
X_Uh = sortrows(X_Uh); % rerange the x-coord from min to max.

plot(X_Uh(:,1),X_Uh(:,2));
grid on
% axis([0 2 0 1])

end % function