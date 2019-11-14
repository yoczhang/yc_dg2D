% Copyright (c) 2010, Thomas-Peter Fries, RWTH Aachen University
function [k1, k2] = StressIntFactorRadial(Mesh, xx, yy, Sol, ShapeFctData, ...
    rr, k1, EE, nu, kappa, mu, lambda, NodeNum, ElemNum)

% Get stress intensity factors.

CasePlot = 0;

CC = [lambda+2*mu  lambda  0; lambda  lambda+2*mu  0; 0  0  4*mu];

if abs(lambda - EE*nu/(1-nu*nu)) < 1.e-12 % Plane stress.
    EEPrime = EE;
elseif abs(lambda - EE*nu/((1+nu)*(1-2*nu))) < 1.e-12 % Plane strain.
    EEPrime = EE / (1-nu^2);
else
    error('Internal error: Neither plane strain nor plane stress!')
end

k1_Mode1 = 1; % Do not change this!
k2_Mode2 = 1; % Do not change this!

% Get circular ring of elements around crack (DEFINED BY rr)
Count = 0;
Dist = sqrt(xx.*xx+yy.*yy);
ffHelp = Dist-rr;
for m = 1 : ElemNum
    Nodes = Mesh(m, :);
    SignVect = sign(ffHelp(Nodes));
    if isempty(find(SignVect==0)) == 0
        error('ffHelp is zero at a node!')
    end
    if min(SignVect)~=max(SignVect)
        Count = Count + 1;
        SpecialElems(Count, 1) = m;
    end
end

if CasePlot == 1
    % Draw special elements.
    PlotMesh(Mesh, xx, yy, ElemNum)
    for i = 1 : length(SpecialElems)
        CurrElem = SpecialElems(i);
        Nodes = Mesh(CurrElem, :);
        patch(xx(Nodes), yy(Nodes), 'r')
    end
    % Draw circle.
    alpha = linspace(0, 2*pi, 100);
    plot(rr*cos(alpha), rr*sin(alpha), 'b-')
end

% Loop over elements.
Sum_k1 = 0;
Sum_k2 = 0;
for m = 1 : length(SpecialElems)

    CurrElem = SpecialElems(m);
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

    % Construct ramp function q.
    PosInner = find(ffHelp(Nodes)<0);

    q = zeros(1, nQ); dqdx = zeros(1, nQ); dqdy = zeros(1, nQ);
    for i = 1 : length(PosInner)
        q    = q    + NN(PosInner(i), :);
        dqdx = dqdx + NNx(PosInner(i), :);
        dqdy = dqdy + NNy(PosInner(i), :);
    end

    % Plot ramp function
    if CasePlot == 1
        plot3(xxInt, yyInt, q, 'k*')
    end

    % Loop over integration points.
    for i = 1 : nQ

        N = NN(:, i); M = MM(:, i); F1 = FF1(:, i); F2 = FF2(:, i); F3 = FF3(:, i); F4 = FF4(:, i);
        Nx = NNx(:, i); Mx = MMx(:, i); F1x = FF1x(:, i); F2x = FF2x(:, i); F3x = FF3x(:, i); F4x = FF4x(:, i);
        Ny = NNy(:, i); My = MMy(:, i); F1y = FF1y(:, i); F2y = FF2y(:, i); F3y = FF3y(:, i); F4y = FF4y(:, i);

        % Get numerical data from approximation.
        uu_Ap = ...
             N' * Sol(Nodes           ) + ...
             M' * Sol(Nodes+ 2*NodeNum) + ...
            F1' * Sol(Nodes+ 4*NodeNum) + ...
            F2' * Sol(Nodes+ 5*NodeNum) + ...
            F3' * Sol(Nodes+ 6*NodeNum) + ...
            F4' * Sol(Nodes+ 7*NodeNum);
        vv_Ap = ...
             N' * Sol(Nodes+   NodeNum) + ...
             M' * Sol(Nodes+ 3*NodeNum) + ...
            F1' * Sol(Nodes+ 8*NodeNum) + ...
            F2' * Sol(Nodes+ 9*NodeNum) + ...
            F3' * Sol(Nodes+10*NodeNum) + ...
            F4' * Sol(Nodes+11*NodeNum);
        dudx_Ap = ...
             Nx' * Sol(Nodes           ) + ...
             Mx' * Sol(Nodes+ 2*NodeNum) + ...
            F1x' * Sol(Nodes+ 4*NodeNum) + ...
            F2x' * Sol(Nodes+ 5*NodeNum) + ...
            F3x' * Sol(Nodes+ 6*NodeNum) + ...
            F4x' * Sol(Nodes+ 7*NodeNum);
        dudy_Ap = ...
             Ny' * Sol(Nodes           ) + ...
             My' * Sol(Nodes+ 2*NodeNum) + ...
            F1y' * Sol(Nodes+ 4*NodeNum) + ...
            F2y' * Sol(Nodes+ 5*NodeNum) + ...
            F3y' * Sol(Nodes+ 6*NodeNum) + ...
            F4y' * Sol(Nodes+ 7*NodeNum);
        dvdx_Ap = ...
             Nx' * Sol(Nodes+   NodeNum) + ...
             Mx' * Sol(Nodes+ 3*NodeNum) + ...
            F1x' * Sol(Nodes+ 8*NodeNum) + ...
            F2x' * Sol(Nodes+ 9*NodeNum) + ...
            F3x' * Sol(Nodes+10*NodeNum) + ...
            F4x' * Sol(Nodes+11*NodeNum);
        dvdy_Ap = ...
             Ny' * Sol(Nodes+   NodeNum) + ...
             My' * Sol(Nodes+ 3*NodeNum) + ...
            F1y' * Sol(Nodes+ 8*NodeNum) + ...
            F2y' * Sol(Nodes+ 9*NodeNum) + ...
            F3y' * Sol(Nodes+10*NodeNum) + ...
            F4y' * Sol(Nodes+11*NodeNum);

        Eps11_Ap = dudx_Ap;
        Eps22_Ap = dvdy_Ap;
        Eps12_Ap = 0.5 * (dudy_Ap + dvdx_Ap);
        Sigma = [CC(1:2,1:3); CC(3,1:3)/2] * [Eps11_Ap; Eps22_Ap; Eps12_Ap]; % Modify C(3,3)-entry by 1/2.
        Sig11_Ap = Sigma(1);
        Sig22_Ap = Sigma(2);
        Sig12_Ap = Sigma(3);

        % Get exact values at integration points according to k1 = 1.
        [uu_Mode1, vv_Mode1, dudx_Mode1, dudy_Mode1, dvdx_Mode1, dvdy_Mode1, Eps11_Mode1, Eps12_Mode1, Eps22_Mode1, ...
            Sig11_Mode1, Sig12_Mode1, Sig22_Mode1] = ExactSol_Mode1(xxInt(i), yyInt(i), k1_Mode1, kappa, mu, lambda);        

        [uu_Mode2, vv_Mode2, dudx_Mode2, dudy_Mode2, dvdx_Mode2, dvdy_Mode2, Eps11_Mode2, Eps12_Mode2, Eps22_Mode2, ...
            Sig11_Mode2, Sig12_Mode2, Sig22_Mode2] = ExactSol_Mode2(xxInt(i), yyInt(i), k2_Mode2, kappa, mu, lambda);        

        % Compute W_12 = ([Eps11_Ap; Eps22_Ap; Eps12_Ap]' * CC * [Eps11_Mode1; Eps22_Mode1; Eps12_Mode1])
        W_12_Mode1 = (Sig11_Ap*Eps11_Mode1 + 2*Sig12_Ap*Eps12_Mode1 + Sig22_Ap*Eps22_Mode1);

        Sum_k1 = Sum_k1 + wwInt(i) * (...
            dqdx(i) * (Sig11_Ap*dudx_Mode1 + Sig12_Ap*dvdx_Mode1 + Sig11_Mode1*dudx_Ap + Sig12_Mode1*dvdx_Ap - W_12_Mode1 ) + ...
            dqdy(i) * (Sig12_Ap*dudx_Mode1 + Sig22_Ap*dvdx_Mode1 + Sig12_Mode1*dudx_Ap + Sig22_Mode1*dvdx_Ap) );

        % Compute W_12 = ([Eps11_Ap; Eps22_Ap; Eps12_Ap]' * CC * [Eps11_Mode2; Eps22_Mode2; Eps12_Mode2])
        W_12_Mode2 = (Sig11_Ap*Eps11_Mode2 + 2*Sig12_Ap*Eps12_Mode2 + Sig22_Ap*Eps22_Mode2);

        Sum_k2 = Sum_k2 + wwInt(i) * (...
            dqdx(i) * (Sig11_Ap*dudx_Mode2 + Sig12_Ap*dvdx_Mode2 + Sig11_Mode2*dudx_Ap + Sig12_Mode2*dvdx_Ap - W_12_Mode2 ) + ...
            dqdy(i) * (Sig12_Ap*dudx_Mode2 + Sig22_Ap*dvdx_Mode2 + Sig12_Mode2*dudx_Ap + Sig22_Mode2*dvdx_Ap) );

    end

end

k1 = EEPrime/2 * Sum_k1;
k2 = EEPrime/2 * Sum_k2;

if CasePlot == 1
    pause
end
