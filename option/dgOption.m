function option = dgOption(option)
%
%   DG method default option.
%
%
%   YcZhang 12/8/2017
%
%   Last modified 15/8/2017
%

if ~isfield(option,'basesType_trial')
    option.elemType = 'P2';    
end

if ~isfield(option,'basesType_test')
    option.elemType = 'P2';    
end

if ~isfield(option,'basesType_u')
    option.elemType = 'P2';    
end

if ~isfield(option,'basesType_p')
    option.elemType = 'P1';    
end

if ~isfield(option,'faceType')
    option.faceType = 'P1';    
end

if ~isfield(option,'export')
    option.export = 0; % 0 stands for do not export system, solvers.
end

if ~isfield(option,'maxIt') % to define the maximum circulation, to generate the Rate.
    option.maxIt = 4;    
end

if ~isfield(option,'maxN') % the maximum number of nodes
    option.maxN = 2e5;    
end

if ~isfield(option,'L0') % to uniformrefine the mesh and generate an initial mesh
    option.L0 = 0;    
end

if ~isfield(option,'refType')
    option.refType = 'red';
end

if ~isfield(option,'printlevel')
    option.printlevel = 1;    
end

if ~isfield(option,'plotflag')
    option.plotflag = 1;    
end

if ~isfield(option,'rateflag')
    option.rateflag = 1;    
end

if ~isfield(option,'dispflag')
    option.dispflag = 1;    
end

if ~isfield(option,'tol')
    option.tol = 1e-8;    
end

end % function