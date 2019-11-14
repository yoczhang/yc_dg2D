% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [L2Norm] = GetL2Norm_RodBiMat(ElemNum, NodeNum, Mesh, ShapeFctData, ...
    Sol, EE1, EE2, LineLoadRight, TestCase, cc)

% Compute L2-norm according to the bi-material rod problem.

reset(cla), reset(clf), hold on

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

        if TestCase == 1
            if cc ~= -0.5
                cc
                error('The exact solution for TestCase 1 is only specified for a disc. at x = 0.5')
            end
            if ffInt(i) < 0
                uuEx = -(xxInt(i)-0.5).^2/(2*EE2)-(xxInt(i)-0.5)*(EE1-EE2)/(4*EE2*(EE2+EE1))+1/(4*(EE2+EE1));
                vvEx = 0;
            else    
                uuEx = -(xxInt(i)-0.5).^2/(2*EE1)-(xxInt(i)-0.5)*(EE1-EE2)/(4*EE1*(EE2+EE1))+1/(4*(EE2+EE1));
                vvEx = 0;
            end
        elseif TestCase == 2
            if ffInt(i) < 0
                uuEx = LineLoadRight/EE2*xxInt(i);
                vvEx = 0;
            else    
                uuEx = (LineLoadRight*EE2*(xxInt(i)-(-cc)) + LineLoadRight*EE1*(-cc)) / (EE1*EE2);;
                vvEx = 0;
            end
        else
            error('This value for TestCase is not defined!')
        end

        uuError = (uuEx - uuAp) * (uuEx - uuAp);
        vvError = (vvEx - vvAp) * (vvEx - vvAp);

        uuL2Norm = uuL2Norm + wwInt(i) * uuError;
        vvL2Norm = vvL2Norm + wwInt(i) * vvError;

        uuFact = uuFact + wwInt(i) * uuEx * uuEx;
        vvFact = vvFact + wwInt(i) * vvEx * vvEx;

        plot3(xxInt(i), yyInt(i), uuAp, 'g*')
        plot3(xxInt(i), yyInt(i), uuEx , 'b*')

    end

    if (CurrElem/100) == round(CurrElem/100)
        disp([num2str(CurrElem) ' / ' num2str(ElemNum)])
    end

end

uuL2Norm = sqrt( uuL2Norm ) / sqrt( uuFact );
vvL2Norm = 0;%sqrt( vvL2Norm ) / sqrt( vvFact ); % Note: sqrt( vvFact )=0 => vvL2Norm=Inf
L2Norm = uuL2Norm + vvL2Norm;
