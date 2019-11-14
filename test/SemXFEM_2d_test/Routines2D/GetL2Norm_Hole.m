% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [L2Norm, EnergyNorm] = GetL2Norm_Hole(ElemNum, NodeNum, Mesh, ShapeFctData, ...
    Sol, aa, mu1, lambda1, EE1, nu1, kappa1)

% Compute L2-norm and energy norm according to the plate with hole.

% reset(cla), reset(clf), hold on

uuL2Norm = 0; vvL2Norm = 0; EnergyNorm = 0;
uuFact   = 0; vvFact   = 0; EnergyFact = 0;

CC = [lambda1+2*mu1  lambda1  0; lambda1  lambda1+2*mu1  0; 0  0  4*mu1];

% Loop over elements.
for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);

    xxInt = ShapeFctData(CurrElem).xxInt;
    yyInt = ShapeFctData(CurrElem).yyInt;
    wwInt = ShapeFctData(CurrElem).wwInt;
    ffInt = ShapeFctData(CurrElem).ffInt;
    NN = ShapeFctData(CurrElem).N ; NNx = ShapeFctData(CurrElem).Nx ; NNy = ShapeFctData(CurrElem).Ny ;
    MM = ShapeFctData(CurrElem).M ; MMx = ShapeFctData(CurrElem).Mx ; MMy = ShapeFctData(CurrElem).My ;
    nQ = ShapeFctData(CurrElem).nQ;

    % Loop over integration points.
    for i = 1 : nQ 

        if ffInt(i) < 0
            continue
        end

        N  =  NN(:, i); M  =  MM(:, i);
        Nx = NNx(:, i); Mx = MMx(:, i);
        Ny = NNy(:, i); My = MMy(:, i);
        uuAp = N' * Sol(Nodes        ) + M' * Sol(Nodes+2*NodeNum);
        vvAp = N' * Sol(Nodes+NodeNum) + M' * Sol(Nodes+3*NodeNum);
        epsXXApprox = ... % dudx
             Nx' * Sol(Nodes           ) + ...
             Mx' * Sol(Nodes+2*NodeNum);
        epsYYApprox = ... % dvdy
             Ny' * Sol(Nodes+  NodeNum) + ...
             My' * Sol(Nodes+3*NodeNum);
        epsXYApprox = 0.5 * (... % 0.5*(dudy+dvdx)
             Ny' * Sol(Nodes          ) +  Nx' * Sol(Nodes+  NodeNum) + ...
             My' * Sol(Nodes+2*NodeNum) +  Mx' * Sol(Nodes+3*NodeNum) );

        if ffInt(i) > 0
            [uuEx, vvEx, Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
                ExactSol_Hole(xxInt(i), yyInt(i), aa, kappa1, lambda1, mu1);
        else
            uuEx = 0;
            vvEx = 0;
        end

        EpsAp = [epsXXApprox-Eps11; epsYYApprox-Eps22; epsXYApprox-Eps12];
        EpsEx = [Eps11; Eps22; Eps12];

        uuError = (uuEx - uuAp) * (uuEx - uuAp);
        vvError = (vvEx - vvAp) * (vvEx - vvAp);

        uuL2Norm = uuL2Norm + wwInt(i) * uuError;
        vvL2Norm = vvL2Norm + wwInt(i) * vvError;

        uuFact = uuFact + wwInt(i) * uuEx * uuEx;
        vvFact = vvFact + wwInt(i) * vvEx * vvEx;

        EnergyNorm = EnergyNorm + wwInt(i) * (EpsAp' * CC * EpsAp);
        EnergyFact = EnergyFact + wwInt(i) * (EpsEx' * CC * EpsEx);

        %plot3(xxInt(i), yyInt(i), vvAp, 'g*')
        %plot3(xxInt(i), yyInt(i), vvEx , 'b*')

    end

    if (CurrElem/100) == round(CurrElem/100)
        disp([num2str(CurrElem) ' / ' num2str(ElemNum)])
    end

end

uuL2Norm = sqrt( uuL2Norm ) / sqrt( uuFact );
vvL2Norm = sqrt( vvL2Norm ) / sqrt( vvFact );
L2Norm = uuL2Norm + vvL2Norm;

EnergyNorm = sqrt( EnergyNorm ) / sqrt( EnergyFact );
