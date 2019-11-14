% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotLevelSet
% This function plots the PHI, PSI and ZETA level set functions.

global CHI CRACK DOMAIN INC PHI PLOT PSI VOID XYZ ZETA
 
nXElem = DOMAIN(1);                                                         % Number of elements in the x-direction
nYElem = DOMAIN(2);                                                         % Number of elements in the y-direction
lXElem = DOMAIN(3);                                                         % Length of elements in the x-direction

X = zeros((nXElem+1),(nYElem+1));
Y = X; Z = X; Zp = X;

% Plot the level set function representing the crack
if (PLOT(1,2) == 1) || (PLOT(1,3) == 1)
    if PLOT(1,2) == 1
        nNode = 1;
        for iYNode = 1:(nYElem+1)
            for iXNode = 1:(nXElem+1)
                X(iXNode,iYNode)  = XYZ(nNode,2);
                Y(iXNode,iYNode)  = XYZ(nNode,3);
                if PHI(nNode) == 0
                    Zp(iXNode,iYNode) = NaN;
                else
                    Zp(iXNode,iYNode) = PHI(nNode);
                end
                nNode = nNode+1;
            end
        end
        figure; hold on; contourf(X,Y,Zp,40,'LineStyle','none');
        axis equal; axis([0 nXElem*lXElem 0 nYElem*lXElem]); 
        set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
        colorbar('horiz'); title('\phi Level Set');
        
        if PLOT(1,6) == 1
            nPt = size(CRACK,1);
            for iPt = 2:nPt
                x = [CRACK(iPt-1,1) CRACK(iPt,1)];
                y = [CRACK(iPt-1,2) CRACK(iPt,2)];
                plot(x,y,'w','LineWidth',2)
            end
        end
        
        hold off;
    end
    
    if PLOT(1,3) == 1
        nNode = 1;
        for iYNode = 1:(nYElem+1)
            for iXNode = 1:(nXElem+1)
                X(iXNode,iYNode)  = XYZ(nNode,2);
                Y(iXNode,iYNode)  = XYZ(nNode,3);
                if PSI(nNode) == 0
                    Zp(iXNode,iYNode) = NaN;
                else
                    Zp(iXNode,iYNode) = PSI(nNode);
                end
                nNode = nNode+1;
            end
        end
        figure; hold on; contourf(X,Y,Zp,40,'LineStyle','none');
        
        if PLOT(1,6) == 1
            nPt = size(CRACK,1);
            for iPt = 2:nPt
                x = [CRACK(iPt-1,1) CRACK(iPt,1)];
                y = [CRACK(iPt-1,2) CRACK(iPt,2)];
                plot(x,y,'w','LineWidth',2)
            end
        end
        
        axis equal; axis([0 nXElem*lXElem 0 nYElem*lXElem]);
        set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
        colorbar('horiz'); title('\psi Level Set');
    end
end

% Plot the level set function representing the inclusions
if PLOT(1,4) == 1
    figure; hold on;
    nNode = 1;
    for iYNode = 1:(nYElem+1)
        for iXNode = 1:(nXElem+1)
            X(iXNode,iYNode) = XYZ(nNode,2);
            Y(iXNode,iYNode) = XYZ(nNode,3);
            Z(iXNode,iYNode) = ZETA(nNode);
            nNode = nNode+1;
        end
    end
    contourf(X,Y,Z,40,'LineStyle','none'); 
    axis equal; axis([0 nXElem*lXElem 0 nYElem*lXElem]);
    set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
    colorbar('horiz'); title('\zeta Level Set');
    
    if PLOT(1,6) == 1
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
    
    hold off;
end

% Plot the level set function representing the voids
if PLOT(1,5) == 1
    figure; hold on;
    nNode = 1;
    for iYNode = 1:(nYElem+1)
        for iXNode = 1:(nXElem+1)
            X(iXNode,iYNode) = XYZ(nNode,2);
            Y(iXNode,iYNode) = XYZ(nNode,3);
            Z(iXNode,iYNode) = CHI(nNode);
            nNode = nNode+1;
        end
    end
    contourf(X,Y,Z,40,'LineStyle','none'); 
    axis equal; axis([0 nXElem*lXElem 0 nYElem*lXElem]);
    set(gca,'XTick',[],'YTick',[],'XColor','w','YColor','w')
    colorbar('horiz'); title('\chi Level Set');
    
    if PLOT(1,6) == 1
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
    
    hold off;
end