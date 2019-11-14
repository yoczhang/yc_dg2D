%% Tecplot
function exportDGscheme2tecplot(mainK, maxIt, option, setpath_pwd, meshName, ...
    S_DGM, S_DGT, S_uh, S_DGuh1, S_DGuh2, S_DGph, ...
    D_DGM, D_DGT, D_uh, D_DGphix, D_DGphiy, D_DGphi)
    %--- speed
    filenameSpeed=[setpath_pwd,'/outputmat/','tecplotData_',meshName,'_[',option.uTbasestype,...
        ']_k=',num2str(mainK,'%1.e'),'_maxIt',num2str(maxIt),'_Speed.dat'];
    f1=fopen(filenameSpeed,'w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "U" "U1" "U2" \n'); 
    N_Snodes=length(S_uh); N_Dnodes = length(D_uh);
    N_Selems=size(S_DGT,1); N_Delems = size(D_DGT,1);
    SmaxInd = max(max(S_DGT));
    fprintf(f1,'ZONE N=%d,E=%d\n',N_Snodes+N_Dnodes,N_Selems+N_Delems); 
    fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
    for n_Snodes=1:N_Snodes
       fprintf(f1,'%f\t%f\t%f\t%f\t%f\r\n',S_DGM(n_Snodes,1),S_DGM(n_Snodes,2),...
           S_uh(n_Snodes),S_DGuh1(n_Snodes),S_DGuh2(n_Snodes));  
    end
    for n_Dnodes = 1:N_Dnodes
        fprintf(f1,'%f\t%f\t%f\t%f\t%f\r\n',D_DGM(n_Dnodes,1),D_DGM(n_Dnodes,2),...
           D_uh(n_Dnodes),-D_DGphix(n_Dnodes),-D_DGphiy(n_Dnodes));  
    end
    fprintf('\n');
    
    for n_Selems =1:N_Selems
       fprintf(f1,'%d\t%d\t%d\r\n',...
           S_DGT(n_Selems,1),S_DGT(n_Selems,2),S_DGT(n_Selems,3)); 
    end
    for n_Delems =1:N_Delems
       fprintf(f1,'%d\t%d\t%d\r\n',...
           SmaxInd+D_DGT(n_Delems,1),SmaxInd+D_DGT(n_Delems,2),SmaxInd+D_DGT(n_Delems,3)); 
    end
    fclose(f1);
    
    %--- pressure
    filenamePressure=[setpath_pwd,'/outputmat/','tecplotData_',meshName,'_[',option.uTbasestype,...
        ']_k=',num2str(mainK,'%1.e'),'_maxIt',num2str(maxIt),'_Pressure.dat'];
    f1=fopen(filenamePressure,'w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "P"  \n'); 
    N_Snodes=length(S_DGph); N_Dnodes = length(D_DGphi);
    N_Selems=size(S_DGT,1); N_Delems = size(D_DGT,1);
    fprintf(f1,'ZONE N=%d,E=%d\n',N_Snodes+N_Dnodes,N_Selems+N_Delems); 
    fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
    for n_Snodes = 1:N_Snodes
       fprintf(f1,'%f\t%f\t%f\t\r\n',S_DGM(n_Snodes,1),S_DGM(n_Snodes,2),...
           S_DGph(n_Snodes));  
    end
    for n_Dnodes = 1:N_Dnodes
       fprintf(f1,'%f\t%f\t%f\t\r\n',D_DGM(n_Dnodes,1),D_DGM(n_Dnodes,2),...
           D_DGphi(n_Dnodes));  
    end
    fprintf('\n');
    
    for n_Selems =1:N_Selems
       fprintf(f1,'%d\t%d\t%d\r\n',...
           S_DGT(n_Selems,1),S_DGT(n_Selems,2),S_DGT(n_Selems,3)); 
    end
    for n_Delems =1:N_Delems
       fprintf(f1,'%d\t%d\t%d\r\n',...
           SmaxInd+D_DGT(n_Delems,1),SmaxInd+D_DGT(n_Delems,2),SmaxInd+D_DGT(n_Delems,3)); 
    end
    fclose(f1);
    
    %--- only Darcy pressure
    filenamePressure=[setpath_pwd,'/outputmat/','tecplotData_',meshName,'_[',option.uTbasestype,...
        ']_k=',num2str(mainK,'%1.e'),'_maxIt',num2str(maxIt),'_DarcyPressure.dat'];
    f1=fopen(filenamePressure,'w'); fprintf(f1,'TITLE = ""\n'); fprintf(f1,'VARIABLES = "X" "Y" "P"  \n'); 
    N_Dnodes = length(D_DGphi);
    N_Delems = size(D_DGT,1);
    fprintf(f1,'ZONE N=%d,E=%d\n',N_Dnodes,N_Delems); 
    fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
    for n_Dnodes = 1:N_Dnodes
       fprintf(f1,'%f\t%f\t%f\t\r\n',D_DGM(n_Dnodes,1),D_DGM(n_Dnodes,2),...
           D_DGphi(n_Dnodes));
    end
    fprintf('\n');
    
    for n_Delems =1:N_Delems
       fprintf(f1,'%d\t%d\t%d\r\n',...
           D_DGT(n_Delems,1),D_DGT(n_Delems,2),D_DGT(n_Delems,3)); 
    end
    fclose(f1);
    
    
    
end % function