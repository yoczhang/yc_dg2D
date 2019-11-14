clear
clc
close all

fix_x = -0.6:0.05:0.6;
fix_y = -fix_x;
fix_x=fix_x';
fix_y=fix_y';
fixP=[fix_x,fix_y];
fixP=[fixP;[-1,-1;-1,1;1,-1;1,1]];

fd = @(p) -min(min(min(1+p(:,2),1-p(:,2)),1+p(:,1)),1-p(:,1));
fh = @(p) ones(size(p,1),1);
[p,t] = distmesh( fd, fh, 0.2, [-1,-1;1,1], fixP );
patch( 'vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] )