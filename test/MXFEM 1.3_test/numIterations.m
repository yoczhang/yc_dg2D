% Written By: Matthew Jon Pais, University of Florida (2011)
% Website: www.matthewpais.com
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function iter = numIterations

global GROW

if     isempty(GROW) == 1, iter = 1;                % No crack growth
elseif  length(GROW) == 2, iter = GROW(1);          % Incremental crack growth
elseif  length(GROW) == 4, iter = GROW(1)/GROW(2);  % Classical Paris model crack growth
elseif  length(GROW) == 7, iter = GROW(1)/GROW(2);  % Modified Paris model crack growth
end