% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotDeformation(DISP)
% This function plots the deformed geometry.

global CONNEC CRACK DOMAIN NODES PLOT PSI XYZ

nXElem  = DOMAIN(1);                                                        % Number of elements in the x-direction
nYElem  = DOMAIN(2);                                                        % Number of elements in the y-direction
XYZ2    = XYZ;                                                              % Store undeformed XYZ information
dispMax = max(max(DISP));                                                   % The maximum displacement
dispMin = min(min(DISP));                                                   % The minimum displacement
globMax = max(abs(dispMax),abs(dispMin));                                   % The maximum value of displacement

if nnz(PLOT(3,6)) == 0
    exp   = abs(fix(log10(globMax)))-0.5;                                   % Magnitude of scaling for plot
    scale = 10^exp;                                                         % Scaling factor    
else
    scale = PLOT(3,6);                                                      % User-defined scaling factor
end

% Calculate the displaced XYZ information
nNode = 1;
for iYNode = 1:(nYElem+1)
    for iXNode = 1:(nXElem+1)
        XYZ2(nNode,2) = XYZ2(nNode,2)+scale*DISP(2*nNode-1);                % Find the displacement in the x-direction
        XYZ2(nNode,3) = XYZ2(nNode,3)+scale*DISP(2*nNode);                  % Find the displacement in the y-direction
        nNode = nNode+1;
    end
end

% Plot the deformed and undeformed geometries
xf = zeros(4,nXElem*nYElem); yf = xf;
for iElem = 1:(nXElem*nYElem)
    NN = CONNEC(iElem,2:5);                                                 % Nodes for current element
    xf(:,iElem) = XYZ2(NN',2);                                              % Final x-coordinates of nodes
    yf(:,iElem) = XYZ2(NN',3);                                              % Final y-coordinates of nodes
end

figure; hold on;
patch(xf,yf,[0.8 0.8 0.8])
title('Deformed Geometry')
xlabel('X'); ylabel('Y'); axis off; axis equal;

% Plot the node numbers
if PLOT(3,2) == 1
    for i = 1:(nXElem+1)*(nYElem+1)
        text(XYZ2(i,2),XYZ2(i,3),num2str(i));
    end
end

% Plot the element numbers
if PLOT(3,3) == 1
    for i = 1:nXElem*nYElem
        NN = CONNEC(i,2:5);
        XN = XYZ2(NN,2);
        YN = XYZ2(NN,3);
        text(mean(XN),mean(YN),num2str(i));
    end
end

% Plot the enriched nodes
if PLOT(3,4) == 1
    for iNode = 1:size(NODES,1)
        if NODES(iNode,2) ~= 0
            plot(XYZ2(iNode,2),XYZ2(iNode,3),'bo')
        elseif NODES(iNode,4) ~= 0
            plot(XYZ2(iNode,2),XYZ2(iNode,3),'bs')
        end
    end
end

% Plot the deformed crack
if PLOT(3,5) == 1
    if isempty(CRACK) == 0
        cplt = NaN(size(NODES,1),2); index = 1;
        for iElem = 1:(nXElem*nYElem)
            N1  = CONNEC(iElem,2);
            N2  = CONNEC(iElem,3);
            N3  = CONNEC(iElem,4);
            N4  = CONNEC(iElem,5);
            NN  = NODES([N1 N2 N3 N4]',:);
            HEN = length(setdiff(NN(:,2),0));
            CTN = length(setdiff(NN(:,4),0));
            NEN = HEN+CTN;
            if NEN == 4
                X1 = XYZ2(N1,2); X2 = XYZ2(N2,2); X3 = XYZ2(N3,2); X4 = XYZ2(N4,2);
                Y1 = XYZ2(N1,3); Y2 = XYZ2(N2,3); Y3 = XYZ2(N3,3); Y4 = XYZ2(N4,3);
                xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];
                psi    = [PSI(N1) PSI(N2) PSI(N3) PSI(N4)];
                corner = [1 2 3 4 1];
                node   = [-1 -1;1 -1;1 1;-1 1];
                U = [DISP(2*N1-1) DISP(2*N1) DISP(2*N2-1) DISP(2*N2)...
                    DISP(2*N3-1) DISP(2*N3) DISP(2*N4-1) DISP(2*N4)];
                for i = 1:4
                    n1 = corner(i);
                    n2 = corner(i+1);
                    if psi(n1)*psi(n2) < 0
                        r    = psi(n1)/(psi(n1)-psi(n2));
                        pnt  = (1-r)*node(n1,:)+r*node(n2,:);
                        xi   = pnt(1); eta = pnt(2);
                        N    = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...
                            (1+xi)*(1+eta);(1-xi)*(1+eta)];
                        xpnt = dot(N,xyz(1:4,1)');
                        ypnt = dot(N,xyz(1:4,2)');
                        upnt = N(1)*U(1)+N(2)*U(3)+N(3)*U(5)+N(4)*U(7);
                        vpnt = N(1)*U(2)+N(2)*U(4)+N(3)*U(6)+N(4)*U(8);
                        cplt(index,1) = xpnt+upnt;
                        cplt(index,2) = ypnt+vpnt;
                        index = index+1;
                    end
                end
            end
        end
        cplt(isnan(cplt(:,1)),:) = [];
        cplt = sortrows(cplt);
        npt  = size(cplt,1);
        plot(cplt([1 npt],1),cplt([1 npt],2),'w','LineWidth',1)
    end
end

hold off