function My_symsDerivate()
clc;
clearvars;

%% --------------- log start ------------------
load('setpath_pwd.mat')
logFilepath = [setpath_pwd,'/symsCompute/'];
logFilename=[logFilepath,'sysmDerivate_log.txt'];
sgc_exist = exist(logFilename, 'file');
if sgc_exist
    delete(logFilename)
end

diary(logFilename);
% diary on; % begin diary

%% set syms vars
syms X1 X2 PI x y muu K k t Kappa xi

%% syms computing
haveTimeDerivate = false;
u_case = 1;
switch u_case
    case 1        
%         phi = @(t,x,y) (exp(y)-exp(-y))*sin(x)*exp(t);
%         u1 = @(t,x,y) (K/PI*sin(2*PI*y)*cos(x))*exp(t);
%         u2 = @(t,x,y) ((-2*K + K/PI^2*(sin(PI*y))^2)*sin(x))*exp(t);
%         p = @(t,x,y) 0*x;
        
%         phi = @(t,x,y) x*(1-x)*y*(1-y)*t;
%         u1 = @(t,x,y) (K/PI*sin(2*PI*y)*cos(x))*exp(t);
%         u2 = @(t,x,y) ((-2*K + K/PI^2*(sin(PI*y))^2)*sin(x))*exp(t);
%         p = @(t,x,y) 0*x;
        

%         phi = @(t,x,y) (y)*sin(PI*x)*cos(t);
%         u1 = @(t,x,y) -sin(PI*y)*cos(PI*x)*cos(t);
%         u2 = @(t,x,y) sin(PI*x)*cos(PI*y)*cos(t);
%         p = @(t,x,y) sin(PI*x)*cos(t);
        
        syms theta lambda_G
        phi = @(x,y) (1-theta)*exp(x) + (1-theta)*lambda_G/2*exp(x) + exp(y) + lambda_G/2;
        u1 = @(x,y) ( (1-theta)*exp(x)+1 )*(exp(y)+lambda_G);
        u2 = @(x,y) (1-theta)*exp(x+y) + exp(y);
        p = @(x,y) xi*(cos(PI) + sin(PI))*cos(PI*y);


%         ks = @(x,y) x^2*(x-1)^2*y^2*(y-1)^2;
%         phi = @(x,y) x+y;
%         %u1 = -diff(ks,y);
%         %u2 = diff(ks,x);
%         u1 = @(x,y) - 2*x^2*y*(x - 1)^2*(y - 1)^2 - x^2*y^2*(2*y - 2)*(x - 1)^2;
%         u2 = @(x,y) 2*x*y^2*(x - 1)^2*(y - 1)^2 + x^2*y^2*(2*x - 2)*(y - 1)^2;
%         p = @(x,y) x^5+y^5-1/3;

%         phi = @(x,y) 2*muu*x+(y-y^2+y^3/3-x*(x-1)*(y-1))/K;
%         u1 = @(x,y) y^2 -2*y+1; 
%         u2 = @(x,y) x^2-x;
%         p = @(x,y) 2*muu*(x+y-1)+1/(3*K);

%         phi = @(x,y) (exp(y)-exp(-y))*sin(x);
%         u1 = @(x,y) K/PI*sin(2*PI*y)*cos(x);
%         u2 = @(x,y) (-2*K + K/PI^2*(sin(PI*y))^2)*sin(x);
%         p = @(x,y) 0*x;

%         phi = @(x,y) (y)*sin(PI*x);
%         u1 = @(x,y) -sin(PI*y)*cos(PI*x);
%         u2 = @(x,y) sin(PI*x)*cos(PI*y);
%         p = @(x,y) sin(PI*x);

%         phi = @(x,y) x*y*(1-x/2)*(1-y/2)*exp(x+y);
%         u1 = @(x,y) K/PI*sin(2*PI*y)*cos(x);
%         u2 = @(x,y) (-2*K + K/PI^2*(sin(PI*y))^2)*sin(x);
%         p = @(x,y) 0*x;

       
        
        if haveTimeDerivate
            disp_case1_haveTime(phi,u1,u2,p)
        else
            disp_case1_noTime(phi,u1,u2,p)
        end
        
        %% judge the integration of pressure on whole domain is zero.
        %--- for example, domain: [0,1]x[0,1].
        %P = (2*x-1)*(2*y-1);
        %intValue = int(int(P,x,[0,1]), y,[0,1])
        
        
end % switch



%% -------------- log end ------------------------
diary off; % close the diary

%% write the bash file
bashFileName = [logFilepath,'bashFile/replaceStr_1.sh'];

%------------------------------------------------------------
if exist(bashFileName, 'file')
    fid_replaceStr=fopen(bashFileName,'w');

    replaceStrs = ['sed -i "s/=/ = /g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);
    
    replaceStrs = ['sed -i "s/0;/0*x; /g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);

    replaceStrs = ['sed -i "s/*/.*/g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);

    replaceStrs = ['sed -i "s/\^/.^/g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);

    replaceStrs = ['sed -i "s/PI/pi/g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);
    
    replaceStrs = ['sed -i "s/muu/mu/g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);

    replaceStrs = ['sed -i "s/x,y)/x,y) /g" ',logFilename];
    fprintf(fid_replaceStr,'%s\n',replaceStrs);
    
    fclose(fid_replaceStr);
    
    %----- change chmod
    changeChmod = ['sudo chmod 777 ',bashFileName];
    system(changeChmod);
end
%------------------------------------------------------------

%% run the sh file to relace: * --> .* ; ^ --> .^ and so on.
 if ismac
     temp = ['sed -i "" "s/-i/-i \"\"/g" ',bashFileName];
     system(temp);
 end

system(bashFileName);

end % function


%% ----------------------- sub function ---------------------------
%-----------------------------------------------------------------------------------
function disp_case1_haveTime(phi,u1,u2,p)
%
%
%   

syms X1 X2 PI x y muu K k t
%------ output phi ------
phit = diff(phi,t);
phix=diff(phi,x);
phiy=diff(phi,y);
phixy=diff(phix,y);
phixx=diff(phix,x);
phiyy=diff(phiy,y);
disp(['phi=',func2str(phi),';']);
disp(['phit=@(t,x,y)',char(phit),';']);
disp(['phix=@(t,x,y)',char(phix),';']);
disp(['phiy=@(t,x,y)',char(phiy),';']);
disp(['phixy=@(t,x,y)',char(phixy),';']);
disp(['phixx=@(t,x,y)',char(phixx),';']);
disp(['phiyy=@(t,x,y)',char(phiyy),';']);
disp(' ')
disp(' ')


%------ output u1 --------
u1t = diff(u1,t);
u1x=diff(u1,x);
u1y=diff(u1,y);
u1xy=diff(u1x,y);
u1xx = diff(u1x,x);
u1yy = diff(u1y,y);
disp(['u1=',func2str(u1),';']);
disp(['u1t=@(t,x,y)',char(u1t),';']);
disp(['u1x=@(t,x,y)',char(u1x),';']);
disp(['u1y=@(t,x,y)',char(u1y),';']);
disp(['u1xy=@(t,x,y)',char(u1xy),';']);
disp(['u1xx=@(t,x,y)',char(u1xx),';']);
disp(['u1yy=@(t,x,y)',char(u1yy),';']);
disp(' ')
disp(' ')

%------ output u2 --------
u2t = diff(u2,t);
u2x=diff(u2,x);
u2y=diff(u2,y);
u2xy=diff(u2x,y);
u2xx = diff(u2x,x);
u2yy = diff(u2y,y);
disp(['u2=',func2str(u2),';']);
disp(['u2t=@(t,x,y)',char(u2t),';']);
disp(['u2x=@(t,x,y)',char(u2x),';']);
disp(['u2y=@(t,x,y)',char(u2y),';']);
disp(['u2xy=@(t,x,y)',char(u2xy),';']);
disp(['u2xx=@(t,x,y)',char(u2xx),';']);
disp(['u2yy=@(t,x,y)',char(u2yy),';']);
disp(' ')
disp(' ')

%------ output p --------
px=diff(p,x);
py=diff(p,y);
disp(['p=',func2str(p),';']);
disp(['px=@(t,x,y)',char(px),';']);
disp(['py=@(t,x,y)',char(py),';']);
disp(' ')
disp(' ')


end 

function disp_case1_noTime(phi,u1,u2,p)
%
%
%   

syms X1 X2 PI x y muu K k t
%------ output phi ------
% phit = diff(phi,t);
phix=diff(phi,x);
phiy=diff(phi,y);
phixy=diff(phix,y);
phixx=diff(phix,x);
phiyy=diff(phiy,y);
disp(['phi=',func2str(phi),';']);
% disp(['phit=@(t,x,y)',char(phit),';']);
disp(['phix=@(x,y)',char(phix),';']);
disp(['phiy=@(x,y)',char(phiy),';']);
disp(['phixy=@(x,y)',char(phixy),';']);
disp(['phixx=@(x,y)',char(phixx),';']);
disp(['phiyy=@(x,y)',char(phiyy),';']);
disp(' ')
disp(' ')


%------ output u1 --------
% u1t = diff(u1,t);
u1x=diff(u1,x);
u1y=diff(u1,y);
u1xy=diff(u1x,y);
u1xx = diff(u1x,x);
u1yy = diff(u1y,y);
disp(['u1=',func2str(u1),';']);
% disp(['u1t=@(t,x,y)',char(u1t),';']);
disp(['u1x=@(x,y)',char(u1x),';']);
disp(['u1y=@(x,y)',char(u1y),';']);
disp(['u1xy=@(x,y)',char(u1xy),';']);
disp(['u1xx=@(x,y)',char(u1xx),';']);
disp(['u1yy=@(x,y)',char(u1yy),';']);
disp(' ')
disp(' ')

%------ output u2 --------
% u2t = diff(u2,t);
u2x=diff(u2,x);
u2y=diff(u2,y);
u2xy=diff(u2x,y);
u2xx = diff(u2x,x);
u2yy = diff(u2y,y);
disp(['u2=',func2str(u2),';']);
% disp(['u2t=@(t,x,y)',char(u2t),';']);
disp(['u2x=@(x,y)',char(u2x),';']);
disp(['u2y=@(x,y)',char(u2y),';']);
disp(['u2xy=@(x,y)',char(u2xy),';']);
disp(['u2xx=@(x,y)',char(u2xx),';']);
disp(['u2yy=@(x,y)',char(u2yy),';']);
disp(' ')
disp(' ')

%------ output p --------
px=diff(p,x);
py=diff(p,y);
disp(['p=',func2str(p),';']);
disp(['px=@(x,y)',char(px),';']);
disp(['py=@(x,y)',char(py),';']);
disp(' ')
disp(' ')


end 