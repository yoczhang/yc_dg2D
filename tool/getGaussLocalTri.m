function [phyGpoints, phyGweights] = getGaussLocalTri(coordTri, formulaGauss2D)
%
%   %-------------------------------------------------
%       Do not modify this function.
%       At the begining, this function is just a sub function in every body-integration m-function,
%       but here we just separating it from each body-integration m-function.   
%   %-------------------------------------------------
%
%   input:
%       coordTri, [3 x 2], each row is the coordinates of the veritices of the triangle.
%       formulaGauss2D, [Npoints x 2], the ref-coordinates of Gauss-points.
%
%   output:
%       phyGpoints, [Npoints x 2], the physical (or local) Gauss coordinates of Gauss-points.
%       phyGweights, [Npoints x 1], the physical (or local) Gauss weights of Gauss-points.
%
%
%   YcZhang 14/9/2017
%
%   Last modified 14/9/2017
%

x1=coordTri(1,1);
y1=coordTri(1,2);
x2=coordTri(2,1);
y2=coordTri(2,2);
x3=coordTri(3,1);
y3=coordTri(3,2);
JacobiTri=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

phyGweights = JacobiTri * formulaGauss2D(:,3);
phyGpoints(:,1)=x1+(x2-x1)*formulaGauss2D(:,1)+(x3-x1)*formulaGauss2D(:,2);
phyGpoints(:,2)=y1+(y2-y1)*formulaGauss2D(:,1)+(y3-y1)*formulaGauss2D(:,2);
end % function getGaussLocalTri