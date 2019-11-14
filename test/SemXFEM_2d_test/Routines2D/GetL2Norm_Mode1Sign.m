% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [L2Norm] = GetL2Norm_Mode1Sign(ElemNum, NodeNum, Mesh, ShapeFctData, ...
    Sol, k1, kappa, mu, lambda)

% Compute L2-norm according to the crack-mode 1 test case.

% reset(cla), reset(clf), hold on

uuL2Norm = 0; vvL2Norm = 0;
uuFact   = 0; vvFact   = 0;

% Loop over elements.
for CurrElem = 1 : ElemNum

    Nodes = Mesh(CurrElem, :);

    xxInt = ShapeFctData(CurrElem).xxInt;
    yyInt = ShapeFctData(CurrElem).yyInt;
    wwInt = ShapeFctData(CurrElem).wwInt;
    ffInt = ShapeFctData(CurrElem).ffInt;
    NN = ShapeFctData(CurrElem).N;
    MM = ShapeFctData(CurrElem).M;
    nQ = ShapeFctData(CurrElem).nQ;

    % Loop over integration points.
    for i = 1 : nQ 

        N = NN(:, i); 
        M = MM(:, i);
        uuAp = N' * Sol(Nodes        ) + M' * Sol(Nodes+2*NodeNum);
        vvAp = N' * Sol(Nodes+NodeNum) + M' * Sol(Nodes+3*NodeNum);

        [uuEx, vvEx, duudx, duudy, dvvdx, dvvdy, Eps11, Eps12, Eps22, Sigma11, Sigma12, Sigma22] = ...
            ExactSol_Mode1(xxInt(i), yyInt(i), k1, kappa, mu, lambda);
 
        uuError = (uuEx - uuAp) * (uuEx - uuAp);
        vvError = (vvEx - vvAp) * (vvEx - vvAp);

        uuL2Norm = uuL2Norm + wwInt(i) * uuError;
        vvL2Norm = vvL2Norm + wwInt(i) * vvError;

        uuFact = uuFact + wwInt(i) * uuEx * uuEx;
        vvFact = vvFact + wwInt(i) * vvEx * vvEx;

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
