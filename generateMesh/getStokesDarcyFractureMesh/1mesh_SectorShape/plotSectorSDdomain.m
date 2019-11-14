%--- subfunction3
function plotSectorSDdomain()
%
%
%
close all
clearvars
clc

initH = 0.001;

%%
%--- plot the whole domain
%- Darcy domain
x1_1 = (0:initH:sqrt(2)/2)';
y1_1 = sqrt(1-x1_1.^2);
y1_2 = (0:initH:sqrt(2)/2)';
x1_2 = sqrt(1-y1_2.^2);

figure
hold on
plot(x1_1, y1_1, '-k', 'LineWidth',1)
plot(x1_2, y1_2, '-k', 'LineWidth',1)

x1_3 = (1:initH:2)';
y1_3 = zeros(size(x1_3));
y1_4 = (1:initH:2)';
x1_4 = zeros(size(y1_4));
plot(x1_3, y1_3, '-k', 'LineWidth',1)
plot(x1_4, y1_4, '-k', 'LineWidth',1)

x1_5 = (0:initH:2*sqrt(2)/2)';
y1_5 = sqrt(4-x1_5.^2);
y1_6 = (0:initH:2*sqrt(2)/2)';
x1_6 = sqrt(4-y1_6.^2);
plot(x1_5, y1_5, '-k', 'LineWidth',1)
plot(x1_6, y1_6, '-k', 'LineWidth',1)

%- plot fractures
fixF1x = (0.3:initH:1)';
fixF1y = 0.1*fixF1x + 1.25;
text('Interpreter','latex','String', ...
    '$$ \gamma_1 $$', ...
    'Position', [1+0.15,0.1+1.27], ...
    'HorizontalAlignment','right','FontSize',16);  

fixF2x = (0.2:initH:0.6)';
fixF2y = -0.3*fixF2x + 1.8;
text('Interpreter','latex','String', ...
    '$$ \gamma_2 $$', ...
    'Position', [0.6+0.15,-0.3*0.6 + 1.8], ...
    'HorizontalAlignment','right','FontSize',16);  

fixF3x = (1.15:initH:1.5)';
fixF3y = -0.5*fixF3x + 1.;
text('Interpreter','latex','String', ...
    '$$ \gamma_3 $$', ...
    'Position', [1.5+0.15, -0.5*1.5 + 1], ...
    'HorizontalAlignment','right','FontSize',16); 

fixF4y = (0.5:initH:0.9)';
fixF4x = 1.3*ones(size(fixF4y));
text('Interpreter','latex','String', ...
    '$$ \gamma_4 $$', ...
    'Position', [1.3+0.08, 0.9+0.1], ...
    'HorizontalAlignment','right','FontSize',16); 

plot(fixF1x, fixF1y, '-k', 'LineWidth',1)
plot(fixF2x, fixF2y, '-k', 'LineWidth',1)
plot(fixF3x, fixF3y, '-k', 'LineWidth',1)
plot(fixF4x, fixF4y, '-k', 'LineWidth',1)

%%
%- Stokes domain
x2_1 = (0:initH:2*sqrt(2)/2)';
y2_1 = sqrt(4-x2_1.^2);
y2_2 = (0:initH:2*sqrt(2)/2)';
x2_2 = sqrt(4-y2_2.^2);
plot(x2_1, y2_1, '-k', 'LineWidth',1)
plot(x2_2, y2_2, '-k', 'LineWidth',1)

x2_3 = (2:initH:3)';
y2_3 = zeros(size(x2_3));
y2_4 = (2:initH:3)';
x2_4 = zeros(size(y2_4));
plot(x2_3, y2_3, '-k', 'LineWidth',1)
plot(x2_4, y2_4, '-k', 'LineWidth',1)

x2_5 = (0:initH:3*sqrt(2)/2)';
y2_5 = sqrt(9-x2_5.^2);
y2_6 = (0:initH:3*sqrt(2)/2)';
x2_6 = sqrt(9-y2_6.^2);
plot(x2_5, y2_5, '-k', 'LineWidth',1)
plot(x2_6, y2_6, '-k', 'LineWidth',1)

%% plot arrow
basex = -0.2;
arrow_1 = arrow( [basex,0],[basex,1], ...
    'TipAngle',12,'BaseAngle',30);
text('Interpreter','latex','String', ...
    '$$ r_1 $$', ...
    'Position', [basex+0.08,1+0.07], ...
    'HorizontalAlignment','right','FontSize',16);  

arrow_2 = arrow( [basex-0.06,0],[basex-0.06,2], ...
    'TipAngle',12,'BaseAngle',30);
text('Interpreter','latex','String', ...
    '$$ r_2 $$', ...
    'Position', [basex-0.06+0.08,2+0.07], ...
    'HorizontalAlignment','right','FontSize',16);  

arrow_3 = arrow( [basex-0.12,0],[basex-0.12,3], ...
    'TipAngle',12,'BaseAngle',30);
text('Interpreter','latex','String', ...
    '$$ r_3 $$', ...
    'Position', [basex-0.12+0.08,3+0.07], ...
    'HorizontalAlignment','right','FontSize',16);  


%% plot text
%- left
text('Interpreter','latex','String', ...
    '$$ \Gamma_{\partial D}^1 $$', ...
    'Position', [basex+0.19,1+0.5], ...
    'HorizontalAlignment','right','FontSize',16);  
text('Interpreter','latex','String', ...
    '$$ \Gamma_{\partial S}^1 $$', ...
    'Position', [basex+0.19,2+0.5], ...
    'HorizontalAlignment','right','FontSize',16); 

%- middle
basex1 = sqrt(2)/2;
text('Interpreter','latex','String', ...
    '$$ \Gamma_{\partial D}^2 $$', ...
    'Position', [basex1+0.07,basex1-0.13], ...
    'HorizontalAlignment','right','FontSize',16); 
text('Interpreter','latex','String', ...
    '$$ \Gamma_{SD} $$', ...
    'Position', [2*basex1+0.4,2*basex1-0.12], ...
    'HorizontalAlignment','right','FontSize',16);  
text('Interpreter','latex','String', ...
    '$$ \Gamma_{\partial S}^2 $$', ...
    'Position', [3*basex1+0.4,3*basex1-0.12], ...
    'HorizontalAlignment','right','FontSize',16);  

%- bottom
text('Interpreter','latex','String', ...
    '$$ \Gamma_{\partial D}^3 $$', ...
    'Position', [1+0.5, basex+0.09], ...
    'HorizontalAlignment','right','FontSize',16);  
text('Interpreter','latex','String', ...
    '$$ \Gamma_{\partial S}^3 $$', ...
    'Position', [2+0.5,basex+0.09], ...
    'HorizontalAlignment','right','FontSize',16); 


%- domain
text('Interpreter','latex','String', ...
    '$$ \Omega_D $$', ...
    'Position', [basex1+0.35,basex1+0.1], ...
    'HorizontalAlignment','right','FontSize',16); 
text('Interpreter','latex','String', ...
    '$$ \Omega_S $$', ...
    'Position', [2*basex1+0.42,2*basex1+0.32], ...
    'HorizontalAlignment','right','FontSize',16);  



%%
axis equal;
axis off;
print('SectorSDdomain.jpg','-djpeg', '-r750')
print('SectorSDdomain.eps','-depsc')

end % function