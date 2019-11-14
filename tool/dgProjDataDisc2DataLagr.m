function [LagrX, LagrY, dataLagr] = dgProjDataDisc2DataLagr(dataDisc, meshInfo, degreek)
%
%
%   We let  NLpoints denote the number of Lagrange-Points on Elem,
%               Nelems denote the number of the elements of Th,
%               Nbases denote the number of LOCAL bases on each element.
%
%   input:
%       dataDisc, [Nbases*Nelems x 1], the DG solution on Elems.
%       elem, [Nelems x 3], to get the number of elems.
%       degreek, polynomial degree.
%
%   output:
%       LagrX, [NLpoints x Nelems], map the Lx on REF elem TO PHY elem.
%       LagrY, [NLpoints x Nelems], map the Ly on REF elem TO PHY elem.
%       dataLagr, [NLpoints x Nelems], converts the representation matrix in the HDG/modal 
%           basis to the respective representation matrix in a Lagrange/nodal basis.
%
%
%   YcZhang 24/9/2017
%
%   Last modified 24/9/2017
%

Nelems = meshInfo.Nelems;
Nbases = nchoosek(degreek+2,2);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% In matlab I only can plot P1, so I comment the switch statement, 
% and just let Lx = [0,1,0]; Ly=[0,0,1]; 
% switch degreek
%     case 0, Lx = 1/3; Ly = 1/3;
%     case 1, Lx = [0,1,0]; Ly=[0,0,1];
%     otherwise, Lx = [0,1,0,1/2,0,1/2]; Ly=[0,0,1,1/2,1/2,0];
% %     otherwise, Lx = [0,1,0,2/3,1/3,0,0,1/3,2/3,1/3]; Ly = [0,0,1,1/3,2/3,2/3,1/3,0,0,1/3];
% end

Lx = [0,1,0]; Ly=[0,0,1]; 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

dataDisc = reshape(dataDisc, Nbases, Nelems); % [Nbases x Nelems]

% get the Lagr-Points on ALL phy elems, let Nelem denote the num of elems,
% and In order to facilitate the calculation we need transpose of LagrX, LagrY.
% Also note xhat, yhat is the COLUMN vectors, but mapRef2Phy required the
% input are ROW vectors.
LagrX = mapRef2Phy(1,Lx,Ly); % [Nelems x NLpoints]
LagrY = mapRef2Phy(2,Lx,Ly); % [Nelems x NLpoints]
LagrX = LagrX'; LagrY = LagrY'; % [NLpoints x Nelems]

[Pb,~,~] = hdgRefBasesOnQuad2D(Lx', Ly', degreek); % [NLpoints x Nbases]

dataLagr = Pb*dataDisc; % [NLpoints x Nelems]


end % function