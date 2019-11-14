close all
clearvars
clc

%% example 1
doplot = 1;
initH = 0.03;
interfaceFix_x = (0:initH:2)';
interfaceFix_y = sqrt(4-interfaceFix_x.^2);
interfaceFix = [interfaceFix_x, interfaceFix_y];

%%
%--- paras
squarDomain1 = [0,0;3,3];
%---
fd1 = @(p) dintersect( drectangle(p,0,3,0,3), ddiff(dcircle(p,0,0,3),dcircle(p,0,0,2)) );
fh1 = @(p) 0.005*ones(size(p,1),1);
pfix1 = [0,0;3,0;3,3;0,3;interfaceFix];
[p1,t1] = distmesh( fd1, fh1, initH, squarDomain1, pfix1 );

if doplot
    patch( 'vertices', p1, 'faces', t1, 'facecolor', [.9, .9, .9] )
    axis tight
    axis equal
end

%%
%--- paras
squarDomain2 = [0,0;2,2];
%---
fd2 = @(p) dintersect( drectangle(p,0,2,0,2), ddiff(dcircle(p,0,0,2),dcircle(p,0,0,1)) );
fh2 = @(p) 0.005*ones(size(p,1),1);
pfix2 = [0,0;2,0;2,2;0,2; interfaceFix];
%- some fixed points for Fractures
fixF1x = (0.3:initH:1)';
% fixF1y = 1.1*ones(size(fixF1x))+0.02;
fixF1y = 0.1*fixF1x + 1.25;
fixF1 = [fixF1x,fixF1y];
pfix2 = [pfix2;fixF1];

fixF2x = (0.2:initH:0.6)';
fixF2y = -0.3*fixF2x + 1.8;
fixF2 = [fixF2x,fixF2y];
pfix2 = [pfix2;fixF2];

%- the intersect Fracture
fixF3x = (1.15:initH:1.5)';
fixF3y = -0.5*fixF3x + 1.;
fixF3 = [fixF3x,fixF3y];
pfix2 = [pfix2;fixF3];

fixF4y = (0.5:initH:0.9)';
fixF4x = 1.3*ones(size(fixF4y));
fixF4 = [fixF4x,fixF4y];
pfix2 = [pfix2;fixF4];
[p2,t2] = distmesh( fd2, fh2, initH, squarDomain2, pfix2 );

if doplot
    patch( 'vertices', p2, 'faces', t2, 'facecolor', [.9, .9, .9] )
    hold on
    plot(interfaceFix_x,interfaceFix_y,'.r')
    plot(fixF1x,fixF1y,'-r')
    plot(fixF2x,fixF2y,'-r')
    plot(fixF3x,fixF3y,'-r')
    plot(fixF4x,fixF4y,'-r')
    axis tight
    axis equal
end

%%



