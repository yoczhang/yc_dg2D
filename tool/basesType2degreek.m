function degreek = basesType2degreek(elemType)
%
%   input:
%       elemType, a character string variable, i.e 'P1', 'P2', ...
%
%   output:
%       degreek, the polynomial degree.
%
%
%
%   YcZhang 14/5/2017
%
%   Last modified 14/5/2017
%
if strcmpi(elemType,'P0') 
    degreek = 0;
elseif strcmpi(elemType,'P1') 
    degreek = 1;
elseif  strcmpi(elemType,'P2') 
    degreek = 2;
elseif strcmpi(elemType,'P3')
    degreek = 3;
elseif strcmpi(elemType,'P4')
    degreek = 4;
elseif strcmpi(elemType,'P5')
    degreek = 5;
elseif strcmpi(elemType,'P6')
    degreek = 6;
else % set default degreek
    degreek = 2;
end

end % function