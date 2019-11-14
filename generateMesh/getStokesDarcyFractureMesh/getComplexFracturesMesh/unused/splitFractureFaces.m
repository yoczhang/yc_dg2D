function splitedFractureFace = splitFractureFaces(meshInfo, originalFracture)
%
%
%
%
%

patchPlotMesh(meshInfo.node,meshInfo.elem);
hold on

fractureFace = meshInfo.fractureFace;
NfracFace = length(fractureFace);
NorigFracture = length(originalFracture);
baryFracCoord = meshInfo.baryEdge(fractureFace,:);

splitedFractureFace = cell(NorigFracture,1);
Eps = 6e-9;
wholetrue = false(NfracFace,1);
originalSplitedFLength = zeros(NorigFracture,1); % to record the length of splitedFractureFace
for countFrac = 1:NorigFracture
    currFrac = originalFracture{countFrac}; 
	%> each elem in originalFracture is the original [x-coord,y-coord]
    %> of each fracture.
    
    currFracShif = circshift(currFrac,-1);
    finaltrue = false(NfracFace,1);
    
    lengthCurrFrac = length(currFrac) - 1;
    for kk = 1:lengthCurrFrac
        d = DisPtToLine(baryFracCoord,currFrac(kk,:),currFracShif(kk,:));
        whichtrue = (d < Eps);
        finaltrue = finaltrue + whichtrue;
     
%         %--- test -----
%         aa = fractureFace( true(NfracFace,1) & whichtrue );
%         aa
%         baryFracCoord_aa = meshInfo.baryEdge(aa,:);
%         plot(baryFracCoord_aa(:,1),baryFracCoord_aa(:,2),'or','MarkerSize',4);
%         %------------
    end
    
    wholetrue = wholetrue + finaltrue;

    splitedFractureFace{countFrac} = fractureFace( true(NfracFace,1) & finaltrue );
    originalSplitedFLength(countFrac) = length(splitedFractureFace{countFrac});
    
    %--- test -----
    splitedFractureFace_kk = splitedFractureFace{countFrac};
    baryFracCoord_kk = meshInfo.baryEdge(splitedFractureFace_kk,:);
    plot(baryFracCoord_kk(:,1),baryFracCoord_kk(:,2),'sb','MarkerSize',4);
    %------------
end

computedfracture = fractureFace(true(NfracFace,1) & wholetrue);
remainingfracture = setdiff(fractureFace,computedfracture);

% %--- test -----
% remain_baryFracCoord = meshInfo.baryEdge(remainingfracture,:);
% plot(remain_baryFracCoord(:,1),remain_baryFracCoord(:,2),'*r','MarkerSize',7);
% %------------

%--- 
% next we need to compute the distance of remaining fracture to the given
% original fracture, then decide the remaining fracture belong to which original fracture.
length_remain = length(remainingfracture);
dis_mat = zeros(length_remain,NorigFracture);
baryFracCoord_remain = meshInfo.baryEdge(remainingfracture,:);

for countFrac = 1:NorigFracture
    currFrac = originalFracture{countFrac}; 
    currFracShif = circshift(currFrac,-1);
    
    lengthCurrFrac = length(currFrac) - 1;
    d_temp = zeros(length_remain,lengthCurrFrac);
    for kk = 1:lengthCurrFrac
        d = DisPtToLine(baryFracCoord_remain,currFrac(kk,:),currFracShif(kk,:));
        d_temp(:,kk) = d;
    end
    
    dis_mat(:,countFrac) = min(d_temp,[],2); % get the minmum of each row
end

% %--- test -----
for ii = 1:length_remain
    remain_baryFracCoord = meshInfo.baryEdge(remainingfracture(ii),:);
    plot(remain_baryFracCoord(:,1),remain_baryFracCoord(:,2),'*r','MarkerSize',7);
end
% %-------------

%---
[~,minLocation] = min(dis_mat,[],2); % get the location of the minmum value on each row.
for ii = 1:length_remain
    loc = minLocation(ii);
    splitedFractureFace{loc} = [splitedFractureFace{loc}; remainingfracture(ii)];
end

% for ii = 1:NorigFracture
%     curr
%     
% end

% %--- test -----
for countFrac = 1:NorigFracture
    fractures = splitedFractureFace{countFrac};
    baryFracCoord = meshInfo.baryEdge(fractures,:);
    plot(baryFracCoord(:,1),baryFracCoord(:,2),'.r','MarkerSize',7);
end
% %------------

end % function