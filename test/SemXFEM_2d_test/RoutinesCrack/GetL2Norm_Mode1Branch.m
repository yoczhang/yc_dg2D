% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [L2Norm, EnergyNorm] = GetL2Norm_Mode1Branch(ElemNum, NodeNum, Mesh, ...
    ShapeFctData, Sol, k1, kappa, mu, lambda)

% Compute L2-norm and energy norm according to the crack-mode 1 test case.

% reset(cla), reset(clf), hold on

uuL2Norm = 0; vvL2Norm = 0; EnergyNorm = 0;
uuFact   = 0; vvFact   = 0; EnergyFact = 0;

CC = [lambda+2*mu  lambda  0; lambda  lambda+2*mu  0; 0  0  4*mu];

% Loop over elements.
for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);

    xxInt = ShapeFctData(CurrElem).xxInt;
    yyInt = ShapeFctData(CurrElem).yyInt;
    wwInt = ShapeFctData(CurrElem).wwInt;
    ffInt = ShapeFctData(CurrElem).ffInt;
    NN  = ShapeFctData(CurrElem).N ; NNx  = ShapeFctData(CurrElem).Nx ; NNy  = ShapeFctData(CurrElem).Ny ;
    MM  = ShapeFctData(CurrElem).M ; MMx  = ShapeFctData(CurrElem).Mx ; MMy  = ShapeFctData(CurrElem).My ;
    FF1 = ShapeFctData(CurrElem).F1; FF1x = ShapeFctData(CurrElem).F1x; FF1y = ShapeFctData(CurrElem).F1y;
    FF2 = ShapeFctData(CurrElem).F2; FF2x = ShapeFctData(CurrElem).F2x; FF2y = ShapeFctData(CurrElem).F2y;
    FF3 = ShapeFctData(CurrElem).F3; FF3x = ShapeFctData(CurrElem).F3x; FF3y = ShapeFctData(CurrElem).F3y;
    FF4 = ShapeFctData(CurrElem).F4; FF4x = ShapeFctData(CurrElem).F4x; FF4y = ShapeFctData(CurrElem).F4y;
    nQ = ShapeFctData(CurrElem).nQ;

    % Loop over integration points.
    for i = 1 : nQ 

        N = NN(:, i); M = MM(:, i); F1 = FF1(:, i); F2 = FF2(:, i); F3 = FF3(:, i); F4 = FF4(:, i);
        Nx = NNx(:, i); Mx = MMx(:, i); F1x = FF1x(:, i); F2x = FF2x(:, i); F3x = FF3x(:, i); F4x = FF4x(:, i);
        Ny = NNy(:, i); My = MMy(:, i); F1y = FF1y(:, i); F2y = FF2y(:, i); F3y = FF3y(:, i); F4y = FF4y(:, i);
        uuAp = ...
             N' * Sol(Nodes           ) + ...
             M' * Sol(Nodes+ 2*NodeNum) + ...
            F1' * Sol(Nodes+ 4*NodeNum) + ...
            F2' * Sol(Nodes+ 5*NodeNum) + ...
            F3' * Sol(Nodes+ 6*NodeNum) + ...
            F4' * Sol(Nodes+ 7*NodeNum);
        vvAp = ...
             N' * Sol(Nodes+   NodeNum) + ...
             M' * Sol(Nodes+ 3*NodeNum) + ...
            F1' * Sol(Nodes+ 8*NodeNum) + ...
            F2' * Sol(Nodes+ 9*NodeNum) + ...
            F3' * Sol(Nodes+10*NodeNum) + ...
            F4' * Sol(Nodes+11*NodeNum);
        epsXXAp = ... % dudx
             Nx' * Sol(Nodes           ) + ...
             Mx' * Sol(Nodes+ 2*NodeNum) + ...
            F1x' * Sol(Nodes+ 4*NodeNum) + ...
            F2x' * Sol(Nodes+ 5*NodeNum) + ...
            F3x' * Sol(Nodes+ 6*NodeNum) + ...
            F4x' * Sol(Nodes+ 7*NodeNum);
        epsYYAp = ... % dvdy
             Ny' * Sol(Nodes+   NodeNum) + ...
             My' * Sol(Nodes+ 3*NodeNum) + ...
            F1y' * Sol(Nodes+ 8*NodeNum) + ...
            F2y' * Sol(Nodes+ 9*NodeNum) + ...
            F3y' * Sol(Nodes+10*NodeNum) + ...
            F4y' * Sol(Nodes+11*NodeNum);
        epsXYAp = 0.5 * (... % 0.5*(dudy+dvdx)
             Ny' * Sol(Nodes           ) +  Nx' * Sol(Nodes+   NodeNum) + ...
             My' * Sol(Nodes+ 2*NodeNum) +  Mx' * Sol(Nodes+ 3*NodeNum) + ...
            F1y' * Sol(Nodes+ 4*NodeNum) + F1x' * Sol(Nodes+ 8*NodeNum) + ...
            F2y' * Sol(Nodes+ 5*NodeNum) + F2x' * Sol(Nodes+ 9*NodeNum) + ...
            F3y' * Sol(Nodes+ 6*NodeNum) + F3x' * Sol(Nodes+10*NodeNum) + ...
            F4y' * Sol(Nodes+ 7*NodeNum) + F4x' * Sol(Nodes+11*NodeNum) );

        [uuEx, vvEx, duudx, duudy, dvvdx, dvvdy, Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
            ExactSol_Mode1(xxInt(i), yyInt(i), k1, kappa, mu, lambda);

        EpsAp = [epsXXAp-Eps11; epsYYAp-Eps22; epsXYAp-Eps12];
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
