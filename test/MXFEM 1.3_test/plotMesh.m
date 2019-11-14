% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotMesh
% This function plots the initial geometry.

global CONNEC CRACK DOMAIN INC NODES PLOT VOID XYZ

nXElem = DOMAIN(1);                                                         % Number of elements in the x-direction
nYElem = DOMAIN(2);                                                         % Number of elements in the y-direction
lXElem = DOMAIN(3);

% Plot the deformed and undeformed geometries
xi = zeros(4,nXElem*nYElem); yi = xi;
for iElem = 1:(nXElem*nYElem)
    NN = CONNEC(iElem,2:5);                                                 % Nodes for current element
    xi(:,iElem) = XYZ(NN',2);                                               % Initial x-coordinates of nodes
    yi(:,iElem) = XYZ(NN',3);                                               % Initial y-coordinates of nodes
end

figure; hold on;
patch(xi,yi,[0.8 0.8 0.8])
title('Finite Element Mesh')
xlabel('X'); ylabel('Y'); axis off; axis equal; 

% Plot the node numbers
if PLOT(2,2) == 1
    for i = 1:(nXElem+1)*(nYElem+1)
        text(XYZ(i,2),XYZ(i,3),num2str(i));
    end
end

% Plot the element numbers
if PLOT(2,3) == 1
    for i = 1:nXElem*nYElem
        NN = CONNEC(i,2:5);
        XN = XYZ(NN,2);
        YN = XYZ(NN,3);
        text(mean(XN),mean(YN),num2str(i));
    end
end

% Plot the enriched nodes
if PLOT(2,4) == 1
    for iNode = 1:size(NODES,1)
        if NODES(iNode,2) ~= 0
            plot(XYZ(iNode,2),XYZ(iNode,3),'bo')
        elseif NODES(iNode,4) ~= 0
            plot(XYZ(iNode,2),XYZ(iNode,3),'bs')
        end
        if NODES(iNode,28) ~= 0
            plot(XYZ(iNode,2),XYZ(iNode,3),'rx')
        end
        if NODES(iNode,30) ~= 0
            plot(XYZ(iNode,2),XYZ(iNode,3),'ko')
        end
    end
end

% Plot the discontinuities
if PLOT(2,5) == 1
    nPt = size(CRACK,1);
    for iPt = 2:nPt
        x = [CRACK(iPt-1,1) CRACK(iPt,1)];
        y = [CRACK(iPt-1,2) CRACK(iPt,2)];
        plot(x,y,'w','LineWidth',2)
    end
    
    if isempty(INC) == 0
        if size(INC,2) == 3
            nINC = size(INC,1);
            xMax = nXElem*lXElem;
            yMax = nYElem*lXElem;
            for i = 1:nINC
                xc = INC(i,1); yc = INC(i,2); rc = INC(i,3);
                theta = (0:256)*2*pi/256;
                xp = rc*cos(theta)+xc;
                yp = rc*sin(theta)+yc;
                for j = 1:length(xp)
                    if (xp(j) > xMax) || (xp(j) < 0)
                        xp(j) = NaN; yp(j) = NaN;
                    elseif (yp(j) > yMax) || (yp(j) < 0)
                        xp(j) = NaN; yp(j) = NaN;
                    end
                end
                plot(xp,yp,'w','LineWidth',2)
            end
        elseif size(INC,2) == 4
            plot([INC(1) INC(3)],[INC(2) INC(4)],'w','LineWidth',2)            
        end
    end
    
    if isempty(VOID) == 0
        nVoid = size(VOID,1);
        xMax = nXElem*lXElem;
        yMax = nYElem*lXElem;
        for i = 1:nVoid
            xc = VOID(i,1); yc = VOID(i,2); rc = VOID(i,3);
            theta = (0:256)*2*pi/256;
            xp = rc*cos(theta)+xc;
            yp = rc*sin(theta)+yc;
            for j = 1:length(xp)
                if (xp(j) > xMax) || (xp(j) < 0)
                    xp(j) = NaN; yp(j) = NaN;
                elseif (yp(j) > yMax) || (yp(j) < 0)
                    xp(j) = NaN; yp(j) = NaN;
                end
            end
            plot(xp,yp,'w','LineWidth',2)
        end
    end
end

% Plot the search radius for the J-domain
if PLOT(2,6) == 1
    nPt = size(CRACK,1);
    if nPt ~= 0
        xCT = CRACK(nPt,1);
        yCT = CRACK(nPt,2);
        r = 4*lXElem;
        theta = linspace(0,360);
        for i = 1:length(theta)
            x = r*cosd(theta);
            y = r*sind(theta);
        end
        x = xCT+x;
        y = yCT+y;
        plot(x,y,'y');
    end
end
hold off