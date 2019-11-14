% Original version of subDomain code written by Nguyen Vinh Phu (2006)
% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [gp,gw,J] = subDomain(npt,chi,psi,zeta,xyz,ctflag,jiflag,CCS)
% This function subdivides enriched elements and determines the guass 
% points and weights to be used in the integration during the assembly of 
% the stiffness matrix.

global CRACK
ntip  = []; xCT = []; yCT = [];

if ctflag == 1
    m   = size(CRACK,1);                                                    % Determine number of data points defining crack
    
    xCTip = [CRACK(m,1) CRACK(1,1)];                                        % X-coordinate of crack tip
    yCTip = [CRACK(m,2) CRACK(1,2)];                                        % Y-coordinate of crack tip        

    xMn = min(xyz(:,1));    xMx = max(xyz(:,1));    xMd = (xMn+xMx)/2;
    yMn = min(xyz(:,2));    yMx = max(xyz(:,2));    yMd = (yMn+yMx)/2;
    for i = 1:length(xCTip)
        if (xCTip(i) >= xMn) && (xCTip(i) <= xMx)
            if (yCTip(i) >= yMn) && (yCTip(i) <= yMx)
                xCT = xCTip(i); yCT = yCTip(i); break
            end
        end
    end
    
    if isempty(xCT) == 1
        for i = 1:length(xCTip)
            d(i) = sqrt((xCTip(i)-xMd)^2+(yCTip(i)-yMd)^2);
        end
        dmin = min(d);
        ind  = find(d == dmin);
        xCT = xCTip(ind); yCT = yCTip(ind);
    end
    
    xyz = [xyz;xCT yCT];
    
    % Map crack tip from real space to parametric space
    le   = xyz(2,1)-xyz(1,1);
    xm   = (xyz(2,1)+xyz(1,1))/2;
    ym   = (xyz(3,2)+xyz(2,2))/2;
    xd   = xCT-xm;
    yd   = yCT-ym;
    xi   = 2*xd/le;
    eta  = 2*yd/le;
    ntip = [xi eta];
end

corner = [1 2 3 4 1];
node   = [-1 -1;1 -1;1 1;-1 1];

% Add crack tip to critical nodes for triangularization
node = [node;ntip];

% Loop through crack level set function
if isempty(psi) == 0
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
            xyz  = [xyz;xpnt ypnt];
            node = [node;pnt];
        end
    end
end

% Loop through inclusion level set function
if isempty(zeta) == 0
    for i = 1:4
        n1 = corner(i);
        n2 = corner(i+1);
        if zeta(n1)*zeta(n2) < 0
            r    = zeta(n1)/(zeta(n1)-zeta(n2));
            pnt  = (1-r)*node(n1,:)+r*node(n2,:);
            xi   = pnt(1); eta = pnt(2);
            N    = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...
                        (1+xi)*(1+eta);(1-xi)*(1+eta)];
            xpnt = dot(N,xyz(1:4,1)');
            ypnt = dot(N,xyz(1:4,2)');
            xyz  = [xyz;xpnt ypnt];
            node = [node;pnt];
        end
    end
end

% Loop through void level set function
if isempty(chi) == 0
    for i = 1:4
        n1 = corner(i);
        n2 = corner(i+1);
        if chi(n1)*chi(n2) < 0
            r    = chi(n1)/(chi(n1)-chi(n2));
            pnt  = (1-r)*node(n1,:)+r*node(n2,:);
            xi   = pnt(1); eta = pnt(2);
            N    = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...
                        (1+xi)*(1+eta);(1-xi)*(1+eta)];
            xpnt = dot(N,xyz(1:4,1)');
            ypnt = dot(N,xyz(1:4,2)');
            xyz  = [xyz;xpnt ypnt];
            node = [node;pnt];
        end
    end
end

% Find the triangles defining the subdomains
warning off MATLAB:delaunayn:DuplicateDataPoints
tri = delaunay(node(:,1),node(:,2),{'Qt','Qbb','Qc','Qz','Q3'});

% Loop over subtriangles to get quadrature points and weights
[q,w] = gauss(npt,'TRI');

% Convert to crack tip coordinate system
if jiflag == 1
    m   = size(CRACK,1);                                                    % Number of segments defining the crack
    mj  = size(xyz,1);                                                      % Determine number of data points defining crack
    xCT = CRACK(m,1);                                                       % X-coordinate of crack tip
    yCT = CRACK(m,2);                                                       % Y-coordinate of crack tip
    xyz(:,1) = xyz(:,1)-xCT;                                                % Translate to crack tip as origin
    xyz(:,2) = xyz(:,2)-yCT;                                                % Translate to crack tip as origin
    for i = 1:mj
        xyzp(i,:) = CCS*xyz(i,:)';                                          % Rotate to crack tip coordinate system
    end
    xyz = xyzp;
end

pt = 1;
for e = 1:size(tri,1)
    coord = node(tri(e,:),:);
    xyzl  = xyz(tri(e,:),:);
    for i = 1:length(w)
        xi = q(i,1); eta = q(i,2);
        N  = [1-xi-eta;xi;eta];
        gp(pt,:) = N'*coord;
        gw(pt,1) = w(i)/2;
        J(pt,:)  = [-xyzl(1,1)+xyzl(2,1) -xyzl(1,2)+xyzl(2,2)...
                    -xyzl(1,1)+xyzl(3,1) -xyzl(1,2)+xyzl(3,2)];
        pt = pt+1;
    end
end