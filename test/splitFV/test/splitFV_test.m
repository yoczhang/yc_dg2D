function splitFV_test

clc;
clearvars;
close all

%---- test 1 --------------------------------
fullpatch.vertices = [2 4; 2 8; 8 4; 8 0; 0 4; 2 6; 2 2; 4 2; 4 0; 5 2; 5 0];
fullpatch.faces = [1 2 3; 1 3 4; 5 6 1; 7 8 9; 11 10 4];
figure, subplot(2,1,1), patch(fullpatch,'facecolor','r'), title('Unsplit mesh'); 

splitpatch = splitFV(fullpatch);
colours = lines(length(splitpatch));
subplot(2,1,2), hold on, title('Split mesh');
for i=1:length(splitpatch)
    patch(splitpatch(i),'facecolor',colours(i,:)); 
end
%----------------------------------------------

%---- test 2 ----------------------------------
n=1;
h_x = 1/2^(n+1); h_y = h_x;
h_partition=[h_x,h_y];
left=0;right=1;bottom=0;top=1;
[M,T]=generate_quad_P_T(left,right,bottom,top,h_partition,1);
fullpatch.vertices = M'; 
fullpatch.faces = T';
splitpatch = splitFV(fullpatch);
%------------------------------------------------

end