clearvars

syms u1P n1P v1P u1M n1M v1M u2P n2P v2P u2M n2M v2M

uPxnP = [u1P*n1P 1/2*(u1P*n2P+u2P*n1P);
    1/2*(u1P*n2P+u2P*n1P) u2P*n2P];

uMxnM = [u1M*n1M 1/2*(u1M*n2M+u2M*n1M);
    1/2*(u1M*n2M+u2M*n1M) u2M*n2M];

vPxnP = [v1P*n1P 1/2*(v1P*n2P+v2P*n1P);
    1/2*(v1P*n2P+v2P*n1P) v2P*n2P];

vMxnM = [v1M*n1M 1/2*(v1M*n2M+v2M*n1M);
    1/2*(v1M*n2M+v2M*n1M) v2M*n2M];

%------------------------------------------------
% n1P=-n1M; n2P=-n2M;
% v2P=0; v2M=0;
% r = uPxnP.*vPxnP + uPxnP.*vMxnM + uMxnM.*vPxnP + uMxnM.*vMxnM;
% r_11 = r(1,1); r_12 = r(1,2); r_21 = r(2,1); r_22 = r(2,2);

%------------------------------------------------
r = uPxnP.*vMxnM;

expand(sum(sum(r)))