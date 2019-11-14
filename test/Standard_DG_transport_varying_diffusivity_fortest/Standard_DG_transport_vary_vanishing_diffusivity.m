function Standard_DG_transport_vary_vanishing_diffusivity

% 05/23/2016

clc
close all
clearvars

tic 
basis_type=1;  N_loc=(basis_type+1)*(basis_type+2)/2;
kappa=1; sigma_F=1;
kappa_tidle=1; sigma_alpha_tidle=1;
delta_t=1e-3;T=1;
theta_time=0;

%% Mesh information
% 512 elements for DGT, 2048 elements for new_DGT
% Firstly, load data
file=fopen('DGM.txt'); points=fscanf(file,'%f'); fclose(file); nppp=size(points); npp=nppp(1)/2; DGM=zeros(npp,2);
for i=1:npp
    DGM(i,1)=points(2*i-1);DGM(i,2)=points(2*i);
end
clear points nppp npp
file=fopen('DGT.txt'); connecticy=fscanf(file,'%f'); fclose(file); neee=size(connecticy); nee=neee(1)/3; DGT=zeros(nee,3);
for i=1:nee
    DGT(i,1)=connecticy(3*i-2);DGT(i,2)=connecticy(3*i-1);DGT(i,3)=connecticy(3*i);
end
clear connecticy neee nee
file=fopen('eipsilon.txt'); temp=fscanf(file,'%f'); fclose(file); n_temp=size(temp); eipsilon=zeros(n_temp(1),1);
for i=1:n_temp
    eipsilon(i,1)=temp(i);
end
DGM=DGM'; DGT=DGT'; eipsilon=eipsilon';

%----------------------- yc test ---------------------------

%-- mesh 1
% [node,elem] = squaremesh([0,2,0,1],1/20);
% bdFlag = setboundary(node,elem,'Dirichlet');
% for k = 1:1
% 	[node,elem,bdFlag] = uniformbisect(node,elem,bdFlag);
% end

%-- mesh 2
% h_x = 1/20; h_y = h_x;
% h=[h_x,h_y];
% left=0;right=2;bottom=0;top=1;
% [M,T]=generate_quad_P_T(left,right,bottom,top,h,1);
% node = M'; elem = T';

%-- mesh 3
% [node, elem] = get_BR_paper_mesh(1/10);

%-- mesh 4
% [node, elem] = triMeshOnPolygonDomain([0,2,2,0],[0,0,1,1],1/20);
load('Delaunaymesh_80times40_[0_2]_[0_1]');

% save delaunaymesh node elem

%-- get meshInfo
meshInfo = polyMeshAuxStructure(node, elem);

%-- get the enriched elems
enrichedElem = getEnrichedElemForExamples(node, elem, 5);

eipsilon = ones(size(elem,1),1);
eipsilon(enrichedElem) = 1e-3;
eipsilon=eipsilon';

[DGM,DGT]=generate_DGM_DGT_DGE_DG_edge_flag(node',elem');
%-------------------------------------------------------------



[DGE,DG_edge_flag,F_i,F_in,F_out,F_HP_i]=generate_DGE_DG_edge_flag(DGM,DGT,eipsilon);

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(12);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(8);

% If necessary,
% [new_DGM,new_DGT,new_eipsilon]=generate_refined_mesh_eipsilon(DGM,DGT,eipsilon,1);

%% Assemble mass and stiffness matrix

[M,A1,A2]=generate_stiffness_matrix_local_DG(DGM,DGT,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,eipsilon,basis_type,N_loc);
A3=generate_stiffness_matrix_local_DG_edge('integrand_1_on_edges',DGT,DGE,DG_edge_flag,Gauss_coefficient_reference_1D,Gauss_point_reference_1D...
                                                  ,F_i,basis_type,N_loc,eipsilon,kappa,sigma_F);
A4=generate_stiffness_matrix_from_convection_term(DGE,DGT,DG_edge_flag,basis_type,N_loc,F_i,F_HP_i,F_out,'integrand_2_on_edges'...
                                                           ,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
% A5=generate_stiffness_matrix_diffusion_upwind_scheme_nu_part(DGT,DGE,DG_edge_flag,basis_type,N_loc,F_HP_i,Gauss_coefficient_reference_1D...
%                                                             ,Gauss_point_reference_1D,'integrand_2_on_edges','integrand_3_on_edges'...
%                                                             ,'integrand_4_on_edges',eipsilon,kappa_tidle);                                                       
A6=generate_stiffness_matrix_diffusion_upwind_scheme_alpha_part('integrand_1_on_edges',DGT,DGE,DG_edge_flag,Gauss_coefficient_reference_1D...
                                                               ,Gauss_point_reference_1D...
                                                               ,F_HP_i,basis_type,N_loc,eipsilon,kappa_tidle,sigma_alpha_tidle);

%% Assemble_load_vector
b=generate_load_vector_inflow_boundary_edges(DGT,DGE,basis_type,N_loc,F_in,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);


%% Begin time loop

left_matrix=M+(1-theta_time)*delta_t*(A1-A2+A3+A4+A6);
right_matrix=M-theta_time*delta_t*(A1-A2+A3+A4+A6);
clear M A1 A2 A3 A4 A5 A6

r_0=zeros(size(DGT,2)*N_loc,1);
for i=1:T/delta_t
    right_term=right_matrix*r_0+b*delta_t;
    
    r=left_matrix\right_term;
    
    if mod(i,100)==0
        uh=zeros(3*size(DGT,2),1);
        for n=1:size(DGT,2)
            uh_local=r((n-1)*N_loc+1:n*N_loc,1);
            uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_type);
            uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_type);
            uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_type);
        end
%         figure(i/100);
%         trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);
%         hold on
%         trisurf(DGT',DGM(1,:)',DGM(2,:)',zeros(3*size(DGT,2),1));
%         shading interp, xlabel('x'), ylabel('y'), colorbar , axis equal, axis off, view(2)
    end
    
    r_0=r;
end

figure
plotUhAlong_fixedY_coord(r_0, 0.545, meshInfo, 1, 0)


% %% Figure
uh=zeros(3*size(DGT,2),1);
for n=1:size(DGT,2)
    uh_local=r((n-1)*N_loc+1:n*N_loc,1);
    uh(3*n-2,1)=fe_solution(DGM(1,DGT(1,n)),DGM(2,DGT(1,n)),uh_local,0,0,basis_type);
    uh(3*n-1,1)=fe_solution(DGM(1,DGT(2,n)),DGM(2,DGT(2,n)),uh_local,0,0,basis_type);
    uh(3*n,1)=fe_solution(DGM(1,DGT(3,n)),DGM(2,DGT(3,n)),uh_local,0,0,basis_type);
end

figure(123456789);
trisurf(DGT',DGM(1,:)',DGM(2,:)',uh);
shading interp, xlabel('x'), ylabel('y'), colorbar , axis equal, axis off, view(2)

time=toc
end % function




%% ------------------------------------------ sub function -------------------------------------------------- %%
%--------------------------------------------------------------------------------------------------------------------%
%%>>-- Begin sub function 1 ---------------------------------------------------------
function enrichedElem = getEnrichedElemAdapTransport(node, elem, testcase)

Nnodes = size(node,1);
figure(2002)
patchPlotMesh(node, elem);
switch testcase
    case 1
        %% get sub-enriched elem and plot
        % get the level-set func1
        aa = 1; bb = 0; cc = -0.5;
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)

        % get the level-set func2
        aa = 1; bb = 0; cc = -1.5;
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0,0.7];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)

        % get the level-set func3
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.6).^2 );
        ff = -(Dist - 0.1);
        x_domain=[0,1]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3)

        % get the level-set func4
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.4).^2 );
        ff = -(Dist - 0.1);
        x_domain=[1,2]; y_domain=[0,1];

        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem4 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem4)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem2; cutElem3; cutElem4]);
        
    case 2
        %% get sub-enriched elem and plot
        % get the level-set func1
        aa = 1; bb = 0; cc = -0.5;
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        
        [cutElem1,~] = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)

        % get the level-set func2
        aa = 1; bb = 0; cc = -1.5;
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem2 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem2]);
        
    case 3
        % get the level-set func1
        xx = node(:,1); yy = node(:,2);
        Dist = sqrt( (xx-1).^2 + (yy-0.5).^2 );
        ff = -(Dist - 0.3);

        levelFunc.ff = ff;

        cutElem1 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)
        
        %% get the all-enriched elem
        enrichedElem = unique(cutElem1);
        
    case 4
        mesh_h = 1./(sqrt(Nnodes)-1);
        %% get sub-enriched elem and plot
        % get the level-set func1
        aa = 1; bb = 0; cc = -(0.5-mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem1,~] = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1)
        
        % get the level-set func1_1
        aa = 1; bb = 0; cc = -(0.5-2*mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem1_1,~] = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem1_1)
        
        % get the level-set func2
        aa = 1; bb = 0; cc = -(0.5+mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem2,~] = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2)
        
        % get the level-set func2_1
        aa = 1; bb = 0; cc = -(0.5+2*mesh_h);
        %aa = 1; bb = -1; cc = 0.1;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;
        
        [cutElem2_1,~] = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem2_1)

        % get the level-set func3
        aa = 1; bb = 0; cc = -(1.5-mesh_h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem3 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem3)
        
        % get the level-set func4
        aa = 1; bb = 0; cc = -(1.5+mesh_h);
        %aa = 1; bb = -1; cc = -0.8;
        ff = levelSetFunc(aa,bb,cc,node);
        x_domain=[0,2]; y_domain=[0.3,1];

        levelFunc.aa = aa;
        levelFunc.bb = bb;
        levelFunc.cc = cc;
        levelFunc.ff = ff;
        levelFunc.x_domain = x_domain;
        levelFunc.y_domain = y_domain;

        cutElem4 = getEnrichedElem(levelFunc, node, elem, 1);
        patchPlotEnrichment(node, elem, cutElem4)
        
        %% get the all-enriched elem
        enrichedElem = unique([cutElem1; cutElem1_1; cutElem2; cutElem2_1; cutElem3; cutElem4]);
        

end % switch

end % function
%%<<-- End sub function 1 ---------------------------------------------------------------


function plotUhOnLine(Uh, y_coord, meshInfo, degreek)
%
%
%   plot the Uh on the given line
%

Nbases = nchoosek(degreek+2,2);
node = meshInfo.node;
elem = meshInfo.elem;

xx = node(:,1); yy = node(:,2);

%% level-func
aa_y = y_coord;
ff = yy - aa_y;
levelFunc.ff = ff;
[cutElem,cutNode] = getEnrichedElem(levelFunc, node, elem, 0);

barycutElem_y = meshInfo.baryElem(cutElem,2);
cutElem = cutElem(barycutElem_y <= y_coord);

hold on
X_points = [];
Uh_points = [];
for CurrElem = 1:length(cutElem)
    CurrCutElem = cutElem(CurrElem);
    
    uh_local = Uh((CurrCutElem-1)*Nbases+1:CurrCutElem*Nbases);

    CurrCutNode = cutNode(cutNode(:,3)==CurrCutElem,1:2);
    
    x_min = min(CurrCutNode(:,1));
    x_max = max(CurrCutNode(:,1));
    
    x_points = (x_min:(x_max-x_min)/100:x_max)';
    y_points = aa_y*ones(length(x_points),1);
    
    uh_points = zeros(length(x_points),1);
    for ii = 1:length(x_points)
        uh_points(ii) = fe_solution(x_points(ii),y_points(ii),uh_local,0,0,1);
    end % for
    
    X_points = [X_points; x_points];
    Uh_points = [Uh_points; uh_points];
    
end % for

X_Uh = [X_points, Uh_points];
X_Uh = sortrows(X_Uh);

plot(X_Uh(:,1),X_Uh(:,2));
grid on
% axis([0 2 0 1])

end % function

