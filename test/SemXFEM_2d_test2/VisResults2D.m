function VisResults2D

% Visualize the results from the 2d test cases.

% Results obtained with XFEM2dRodBiMat_AbsEnr.m (TestCase=1)
hh1 = 2./[9 19 39 79];
Res1 = [5.26343704e-002 2.69755877e-002 1.37464753e-002 6.94621098e-003];

% Results obtained with XFEM2dRodBiMat_SignEnr_1dLM.m (TestCase=1)
hh2 = 2./[9 19 39 79];
Res2 = [2.87408468e-002 6.64317811e-003 1.59756341e-003 3.91763949e-004];

% Results obtained with XFEM2dHoleBiMat_AbsEnrMoes.m
hh3 = 2./[10 20 40 80];
Res3 = [2.54599525e-002 6.37563818e-003 1.70974043e-003 4.17820744e-004];

% Results obtained with XFEM2dHoleBiMat_SignEnr_1dLM.m
hh4 = 2./[10 20 40 80];
Res4 = [3.06109747e-002 1.39492684e-002 5.79551862e-003 8.16358970e-003];

% Results obtained with FEM2dVoid.m
hh5 = 2./[10 20 40 80];
Res5 = [4.47109331e-002 1.17942628e-002 2.95719689e-003 7.42623268e-004];

% Results obtained with XFEM2dVoid_HeaviEnr.m
hh6 = 2./[10 20 40 80];
Res6 = [4.47109331e-002 1.17942628e-002 2.95719689e-003 7.42623268e-004];

% Results obtained with XFEM2dCrack_SignEnr.m
hh7 = 2./[9 19 39 79];
Res7 = [5.51790038e-002 2.40729100e-002 1.12289285e-002 5.41042218e-003];

% Results obtained with XFEM2dCrack_BranchEnr.m
hh8 = 2./[9 19 39 59 99];
Res8 = [1.74139653e-002 5.20600733e-003 1.27397149e-003 5.74311e-004 2.06076e-004];

reset(cla), reset(clf), hold on
a(1) = plot(hh1, Res1, 'b-*');
a(2) = plot(hh2, Res2, 'k-v');
set(a, 'LineWidth', 1)
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
a = legend('XFEM2dRodBiMat_AbsEnr', 'XFEM2dRodBiMat_SignEnr_1dLM', 4);
set(a, 'Interpreter', 'none')
a = line([1.e-2 1.e0], [1.e-4 1.e0]);
set(a, 'LineStyle', ':')
title('2d test case (a), L_2-norms.')
disp('Paused, press button.')
pause

reset(cla), reset(clf), hold on
a(1) = plot(hh3, Res3, 'b-*');
a(2) = plot(hh4, Res4, 'k-v');
set(a, 'LineWidth', 1)
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
a = legend('XFEM2dHoleBiMat_AbsEnrMoes', 'XFEM2dHoleBiMat_SignEnr_1dLM', 4);
set(a, 'Interpreter', 'none')
a = line([1.e-2 1.e0], [1.e-4 1.e0]);
set(a, 'LineStyle', ':')
title('2d test case (b), L_2-norms.')
disp('Paused, press button.')
pause

reset(cla), reset(clf), hold on
a(1) = plot(hh5, Res5, 'b-*');
a(2) = plot(hh6, Res6, 'k-v');
set(a, 'LineWidth', 1)
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
a = legend('FEM2dVoid', 'XFEM2dVoid_HeaviEnr', 4);
set(a, 'Interpreter', 'none')
a = line([1.e-2 1.e0], [1.e-4 1.e0]);
set(a, 'LineStyle', ':')
title('2d test case (c), L_2-norms.')
disp('Paused, press button.')
pause

reset(cla), reset(clf), hold on
a(1) = plot(hh7, Res7, 'b-*');
a(2) = plot(hh8, Res8, 'k-v');
set(a, 'LineWidth', 1)
set(gca, 'XScale', 'log', 'YScale', 'log')
grid on
a = legend('XFEM2dCrack_SignEnr', 'XFEM2dCrack_BranchEnr', 4);
set(a, 'Interpreter', 'none')
a = line([1.e-2 1.e0], [1.e-4 1.e0]);
set(a, 'LineStyle', ':')
title('2d test case (d), L_2-norms.')
