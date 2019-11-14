function r = settingElemFlagAdapTransport(meshInfo, n)
%
%   In this function we want to get the flag of each element, considering
%   there may have the different vector_u in different element.
%
%   input:
%       n, the element index.
%
%   output:
%       r, the flag.
%
%
%   YcZhang 23/9/2017
%
%   Last modified 23/9/2017
%
%

%% Case 1, the input velocity is given by (-1,0).
% r = 1;


%% Case 2, the input velocity is from Stokes-Darcy.
%> Stokes-Darcy domain, y==0 is the interface,
%> 1 stands for Stokes domain,
%> -1 stands for Darcy domain.

if meshInfo.baryElem(n,2) > 0 
    r = 1;
else
    r =-1;
end

end % function