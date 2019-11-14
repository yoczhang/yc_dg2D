function [Pb, Pbx, Pby] = localBases2D(xT, yT, hT, x, y, k)
%
%   In this function we compute the values of 2D local bases at physical local Gauss-Points.
%
%   We let Nbases denote the number of bases, Npoints denote the number of
%   Gauss-Points.
%
%   input:
%       (xT,yT), the centroid(xing xin) (x-coord, y-coord) of element.
%       hT, the diameter of the element.
%       x, the x-coordinates of physical local Gauss-Points, size: 
%                               a row vector, [Npoints x 1].
%       y, the y-coordinates of physical local Gauss-Points, size: 
%                               a row vector, [Npoints x 1].
%       k, the degree of the polynomial, size: a scalar.
%
%   output:
%       Pb, the value of bases at given coordinates, size: 
%                               a matrix, [Npoints x Nbases].
%       Pbx, the derivative to x of the bases, size:
%                               a matrix, [Npoints x Nbases].
%       Pby, the derivative to y of the bases, size:
%                               a matrix, [Npoints x Nbases].
%              
%
%   YcZhang 4/8/2017
%   
%   Last modified 4/8/2017
%

useLocalBasis = 0;
if ~useLocalBasis
    xT = 0;
    yT = 0;
    hT = 1;
end

[r, c] = size(x);
if r < c % we need to guarantee the input x is a column vector.
    x = x'; 
end
[r, c] = size(y);
if r < c % we need to guarantee the input y is a column vector.
    y = y';
end

vone = ones(length(x),1);
if k == 0
    Pb = 1*vone;
    Pbx = 0*vone;
    Pby = 0*vone;
elseif k == 1
    Pb = [1*vone, -(xT - x)/hT, -(yT - y)/hT];
    Pbx = [0*vone, 1/hT*vone, 0*vone];
    Pby = [0*vone, 0*vone, 1/hT*vone];
elseif k == 2
    Pb = [1*vone, -(xT - x)/hT, -(yT - y)/hT, (xT - x).^2/hT.^2, ((xT - x).*(yT - y))/hT.^2, (yT - y).^2/hT.^2];
    Pbx = [0*vone, 1/hT*vone, 0*vone, -(2.*xT - 2.*x)/hT.^2, -(yT - y)/hT.^2, 0*vone];
    Pby = [0*vone, 0*vone, 1/hT*vone, 0*vone, -(xT - x)/hT.^2, -(2.*yT - 2.*y)/hT.^2];
elseif k == 3
    Pb = [ 1.*vone, -(xT - x)/hT, -(yT - y)/hT, (xT - x).^2/hT.^2, ((xT - x).*(yT - y))/hT.^2, (yT - y).^2/hT.^2, -(xT - x).^3/hT.^3, -((xT - x).^2.*(yT - y))/hT.^3, -((xT - x).*(yT - y).^2)/hT.^3, -(yT - y).^3/hT.^3];
    Pbx = [ 0.*vone, 1.*vone/hT, 0.*vone, -(2.*xT - 2.*x)/hT.^2, -(yT - y)/hT.^2, 0.*vone, (3.*(xT - x).^2)/hT.^3, ((2.*xT - 2.*x).*(yT - y))/hT.^3, (yT - y).^2/hT.^3, 0.*vone];
    Pby = [ 0.*vone, 0.*vone, 1.*vone/hT, 0.*vone, -(xT - x)/hT.^2, -(2.*yT - 2.*y)/hT.^2, 0.*vone, (xT - x).^2/hT.^3, ((2.*yT - 2.*y).*(xT - x))/hT.^3, (3.*(yT - y).^2)/hT.^3];
elseif k == 4
    Pb = [ 1.*vone, -(xT - x)/hT, -(yT - y)/hT, (xT - x).^2/hT.^2, ((xT - x).*(yT - y))/hT.^2, (yT - y).^2/hT.^2, -(xT - x).^3/hT.^3, -((xT - x).^2.*(yT - y))/hT.^3, -((xT - x).*(yT - y).^2)/hT.^3, -(yT - y).^3/hT.^3, (xT - x).^4/hT.^4, ((xT - x).^3.*(yT - y))/hT.^4, ((xT - x).^2.*(yT - y).^2)/hT.^4, ((xT - x).*(yT - y).^3)/hT.^4, (yT - y).^4/hT.^4];
    Pbx = [ 0.*vone, 1.*vone/hT, 0.*vone, -(2.*xT - 2.*x)/hT.^2, -(yT - y)/hT.^2, 0.*vone, (3.*(xT - x).^2)/hT.^3, ((2.*xT - 2.*x).*(yT - y))/hT.^3, (yT - y).^2/hT.^3, 0.*vone, -(4.*(xT - x).^3)/hT.^4, -(3.*(xT - x).^2.*(yT - y))/hT.^4, -((2.*xT - 2.*x).*(yT - y).^2)/hT.^4, -(yT - y).^3/hT.^4, 0.*vone];
    Pby = [ 0.*vone, 0.*vone, 1.*vone/hT, 0.*vone, -(xT - x)/hT.^2, -(2.*yT - 2.*y)/hT.^2, 0.*vone, (xT - x).^2/hT.^3, ((2.*yT - 2.*y).*(xT - x))/hT.^3, (3.*(yT - y).^2)/hT.^3, 0.*vone, -(xT - x).^3/hT.^4, -((2.*yT - 2.*y).*(xT - x).^2)/hT.^4, -(3.*(xT - x).*(yT - y).^2)/hT.^4, -(4.*(yT - y).^3)/hT.^4];
elseif k == 5
    Pb = [ 1.*vone, -(xT - x)/hT, -(yT - y)/hT, (xT - x).^2/hT.^2, ((xT - x).*(yT - y))/hT.^2, (yT - y).^2/hT.^2, -(xT - x).^3/hT.^3, -((xT - x).^2.*(yT - y))/hT.^3, -((xT - x).*(yT - y).^2)/hT.^3, -(yT - y).^3/hT.^3, (xT - x).^4/hT.^4, ((xT - x).^3.*(yT - y))/hT.^4, ((xT - x).^2.*(yT - y).^2)/hT.^4, ((xT - x).*(yT - y).^3)/hT.^4, (yT - y).^4/hT.^4, -(xT - x).^5/hT.^5, -((xT - x).^4.*(yT - y))/hT.^5, -((xT - x).^3.*(yT - y).^2)/hT.^5, -((xT - x).^2.*(yT - y).^3)/hT.^5, -((xT - x).*(yT - y).^4)/hT.^5, -(yT - y).^5/hT.^5];
    Pbx = [ 0.*vone, 1.*vone/hT, 0.*vone, -(2.*xT - 2.*x)/hT.^2, -(yT - y)/hT.^2, 0.*vone, (3.*(xT - x).^2)/hT.^3, ((2.*xT - 2.*x).*(yT - y))/hT.^3, (yT - y).^2/hT.^3, 0.*vone, -(4.*(xT - x).^3)/hT.^4, -(3.*(xT - x).^2.*(yT - y))/hT.^4, -((2.*xT - 2.*x).*(yT - y).^2)/hT.^4, -(yT - y).^3/hT.^4, 0.*vone, (5.*(xT - x).^4)/hT.^5, (4.*(xT - x).^3.*(yT - y))/hT.^5, (3.*(xT - x).^2.*(yT - y).^2)/hT.^5, ((2.*xT - 2.*x).*(yT - y).^3)/hT.^5, (yT - y).^4/hT.^5, 0.*vone];
    Pby = [ 0.*vone, 0.*vone, 1.*vone/hT, 0.*vone, -(xT - x)/hT.^2, -(2.*yT - 2.*y)/hT.^2, 0.*vone, (xT - x).^2/hT.^3, ((2.*yT - 2.*y).*(xT - x))/hT.^3, (3.*(yT - y).^2)/hT.^3, 0.*vone, -(xT - x).^3/hT.^4, -((2.*yT - 2.*y).*(xT - x).^2)/hT.^4, -(3.*(xT - x).*(yT - y).^2)/hT.^4, -(4.*(yT - y).^3)/hT.^4, 0.*vone, (xT - x).^4/hT.^5, ((2.*yT - 2.*y).*(xT - x).^3)/hT.^5, (3.*(xT - x).^2.*(yT - y).^2)/hT.^5, (4.*(xT - x).*(yT - y).^3)/hT.^5, (5.*(yT - y).^4)/hT.^5];
    
elseif k == 6
    Pb = [ 1.*vone, -(xT - x)/hT, -(yT - y)/hT, (xT - x).^2/hT.^2, ((xT - x).*(yT - y))/hT.^2, (yT - y).^2/hT.^2, -(xT - x).^3/hT.^3, -((xT - x).^2.*(yT - y))/hT.^3, -((xT - x).*(yT - y).^2)/hT.^3, -(yT - y).^3/hT.^3, (xT - x).^4/hT.^4, ((xT - x).^3.*(yT - y))/hT.^4, ((xT - x).^2.*(yT - y).^2)/hT.^4, ((xT - x).*(yT - y).^3)/hT.^4, (yT - y).^4/hT.^4, -(xT - x).^5/hT.^5, -((xT - x).^4.*(yT - y))/hT.^5, -((xT - x).^3.*(yT - y).^2)/hT.^5, -((xT - x).^2.*(yT - y).^3)/hT.^5, -((xT - x).*(yT - y).^4)/hT.^5, -(yT - y).^5/hT.^5, (xT - x).^6/hT.^6, ((xT - x).^5.*(yT - y))/hT.^6, ((xT - x).^4.*(yT - y).^2)/hT.^6, ((xT - x).^3.*(yT - y).^3)/hT.^6, ((xT - x).^2.*(yT - y).^4)/hT.^6, ((xT - x).*(yT - y).^5)/hT.^6, (yT - y).^6/hT.^6];
    Pbx = [ 0.*vone, 1.*vone/hT, 0.*vone, -(2.*xT - 2.*x)/hT.^2, -(yT - y)/hT.^2, 0.*vone, (3.*(xT - x).^2)/hT.^3, ((2.*xT - 2.*x).*(yT - y))/hT.^3, (yT - y).^2/hT.^3, 0.*vone, -(4.*(xT - x).^3)/hT.^4, -(3.*(xT - x).^2.*(yT - y))/hT.^4, -((2.*xT - 2.*x).*(yT - y).^2)/hT.^4, -(yT - y).^3/hT.^4, 0.*vone, (5.*(xT - x).^4)/hT.^5, (4.*(xT - x).^3.*(yT - y))/hT.^5, (3.*(xT - x).^2.*(yT - y).^2)/hT.^5, ((2.*xT - 2.*x).*(yT - y).^3)/hT.^5, (yT - y).^4/hT.^5, 0.*vone, -(6.*(xT - x).^5)/hT.^6, -(5.*(xT - x).^4.*(yT - y))/hT.^6, -(4.*(xT - x).^3.*(yT - y).^2)/hT.^6, -(3.*(xT - x).^2.*(yT - y).^3)/hT.^6, -((2.*xT - 2.*x).*(yT - y).^4)/hT.^6, -(yT - y).^5/hT.^6, 0.*vone];
    Pby = [ 0.*vone, 0.*vone, 1.*vone/hT, 0.*vone, -(xT - x)/hT.^2, -(2.*yT - 2.*y)/hT.^2, 0.*vone, (xT - x).^2/hT.^3, ((2.*yT - 2.*y).*(xT - x))/hT.^3, (3.*(yT - y).^2)/hT.^3, 0.*vone, -(xT - x).^3/hT.^4, -((2.*yT - 2.*y).*(xT - x).^2)/hT.^4, -(3.*(xT - x).*(yT - y).^2)/hT.^4, -(4.*(yT - y).^3)/hT.^4, 0.*vone, (xT - x).^4/hT.^5, ((2.*yT - 2.*y).*(xT - x).^3)/hT.^5, (3.*(xT - x).^2.*(yT - y).^2)/hT.^5, (4.*(xT - x).*(yT - y).^3)/hT.^5, (5.*(yT - y).^4)/hT.^5, 0.*vone, -(xT - x).^5/hT.^6, -((2.*yT - 2.*y).*(xT - x).^4)/hT.^6, -(3.*(xT - x).^3.*(yT - y).^2)/hT.^6, -(4.*(xT - x).^2.*(yT - y).^3)/hT.^6, -(5.*(xT - x).*(yT - y).^4)/hT.^6, -(6.*(yT - y).^5)/hT.^6];
    
elseif k == 7
    Pb = [ 1.*vone, -(xT - x)/hT, -(yT - y)/hT, (xT - x).^2/hT.^2, ((xT - x).*(yT - y))/hT.^2, (yT - y).^2/hT.^2, -(xT - x).^3/hT.^3, -((xT - x).^2.*(yT - y))/hT.^3, -((xT - x).*(yT - y).^2)/hT.^3, -(yT - y).^3/hT.^3, (xT - x).^4/hT.^4, ((xT - x).^3.*(yT - y))/hT.^4, ((xT - x).^2.*(yT - y).^2)/hT.^4, ((xT - x).*(yT - y).^3)/hT.^4, (yT - y).^4/hT.^4, -(xT - x).^5/hT.^5, -((xT - x).^4.*(yT - y))/hT.^5, -((xT - x).^3.*(yT - y).^2)/hT.^5, -((xT - x).^2.*(yT - y).^3)/hT.^5, -((xT - x).*(yT - y).^4)/hT.^5, -(yT - y).^5/hT.^5, (xT - x).^6/hT.^6, ((xT - x).^5.*(yT - y))/hT.^6, ((xT - x).^4.*(yT - y).^2)/hT.^6, ((xT - x).^3.*(yT - y).^3)/hT.^6, ((xT - x).^2.*(yT - y).^4)/hT.^6, ((xT - x).*(yT - y).^5)/hT.^6, (yT - y).^6/hT.^6, -(xT - x).^7/hT.^7, -((xT - x).^6.*(yT - y))/hT.^7, -((xT - x).^5.*(yT - y).^2)/hT.^7, -((xT - x).^4.*(yT - y).^3)/hT.^7, -((xT - x).^3.*(yT - y).^4)/hT.^7, -((xT - x).^2.*(yT - y).^5)/hT.^7, -((xT - x).*(yT - y).^6)/hT.^7, -(yT - y).^7/hT.^7];
    Pbx = [ 0.*vone, 1.*vone/hT, 0.*vone, -(2.*xT - 2.*x)/hT.^2, -(yT - y)/hT.^2, 0.*vone, (3.*(xT - x).^2)/hT.^3, ((2.*xT - 2.*x).*(yT - y))/hT.^3, (yT - y).^2/hT.^3, 0.*vone, -(4.*(xT - x).^3)/hT.^4, -(3.*(xT - x).^2.*(yT - y))/hT.^4, -((2.*xT - 2.*x).*(yT - y).^2)/hT.^4, -(yT - y).^3/hT.^4, 0.*vone, (5.*(xT - x).^4)/hT.^5, (4.*(xT - x).^3.*(yT - y))/hT.^5, (3.*(xT - x).^2.*(yT - y).^2)/hT.^5, ((2.*xT - 2.*x).*(yT - y).^3)/hT.^5, (yT - y).^4/hT.^5, 0.*vone, -(6.*(xT - x).^5)/hT.^6, -(5.*(xT - x).^4.*(yT - y))/hT.^6, -(4.*(xT - x).^3.*(yT - y).^2)/hT.^6, -(3.*(xT - x).^2.*(yT - y).^3)/hT.^6, -((2.*xT - 2.*x).*(yT - y).^4)/hT.^6, -(yT - y).^5/hT.^6, 0.*vone, (7.*(xT - x).^6)/hT.^7, (6.*(xT - x).^5.*(yT - y))/hT.^7, (5.*(xT - x).^4.*(yT - y).^2)/hT.^7, (4.*(xT - x).^3.*(yT - y).^3)/hT.^7, (3.*(xT - x).^2.*(yT - y).^4)/hT.^7, ((2.*xT - 2.*x).*(yT - y).^5)/hT.^7, (yT - y).^6/hT.^7, 0.*vone];
    Pby = [ 0.*vone, 0.*vone, 1.*vone/hT, 0.*vone, -(xT - x)/hT.^2, -(2.*yT - 2.*y)/hT.^2, 0.*vone, (xT - x).^2/hT.^3, ((2.*yT - 2.*y).*(xT - x))/hT.^3, (3.*(yT - y).^2)/hT.^3, 0.*vone, -(xT - x).^3/hT.^4, -((2.*yT - 2.*y).*(xT - x).^2)/hT.^4, -(3.*(xT - x).*(yT - y).^2)/hT.^4, -(4.*(yT - y).^3)/hT.^4, 0.*vone, (xT - x).^4/hT.^5, ((2.*yT - 2.*y).*(xT - x).^3)/hT.^5, (3.*(xT - x).^2.*(yT - y).^2)/hT.^5, (4.*(xT - x).*(yT - y).^3)/hT.^5, (5.*(yT - y).^4)/hT.^5, 0.*vone, -(xT - x).^5/hT.^6, -((2.*yT - 2.*y).*(xT - x).^4)/hT.^6, -(3.*(xT - x).^3.*(yT - y).^2)/hT.^6, -(4.*(xT - x).^2.*(yT - y).^3)/hT.^6, -(5.*(xT - x).*(yT - y).^4)/hT.^6, -(6.*(yT - y).^5)/hT.^6, 0.*vone, (xT - x).^6/hT.^7, ((2.*yT - 2.*y).*(xT - x).^5)/hT.^7, (3.*(xT - x).^4.*(yT - y).^2)/hT.^7, (4.*(xT - x).^3.*(yT - y).^3)/hT.^7, (5.*(xT - x).^2.*(yT - y).^4)/hT.^7, (6.*(xT - x).*(yT - y).^5)/hT.^7, (7.*(yT - y).^6)/hT.^7];
    
end








end % function
