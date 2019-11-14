clearvars
syms u1 u2 V1 V2 N1 N2 PxU1 PxU2 PyU1 PyU2
DU=1/2*[2*PxU1, PyU1+PxU2;
    PyU1+PxU2, 2*PyU2];
OV=[V1*N1, V1*N2;
    V2*N1, V2*N2];

r1 = DU.*OV;

r2=sum(sum(r1));

r3 = expand(r2)

r4 = latex(r3)