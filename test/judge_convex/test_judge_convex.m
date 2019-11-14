function test_judge_convex

clearvars;
close all;
clc;

n=20;
p=rand(n,2);

p=createSimplyPoly(p);  %创建简单多边形
p1 = p;

%%
p1_forward = zeros(size(p1));
p1_forward(1,:) = p1(n,:);
p1_forward(2:n,:) = p1(1:n-1,:);

p1_backward = zeros(size(p1));
p1_backward(1:n-1,:) = p1(2:n,:);
p1_backward(n,:) = p1(1,:);

vec1 = p1_forward - p1;
vec2 = p1_backward - p1;

figure
hold on;
for i=1:n    
    if i==1                     %处理第一个点
        v1=p(n,:)-p(1,:);       %当前点到前一点向量
        v2=p(2,:)-p(1,:);       %当前点到后一点向量
    elseif i==n                 %最后一个点
        v1=p(n-1,:)-p(n,:);    
        v2=p(1,:)-p(n,:);        
    else                        %其他点
        v1=p(i-1,:)-p(i,:);     
        v2=p(i+1,:)-p(i,:);
    end
    r=det([v1;v2]);                 %叉乘后第三个向量的方向
    r1 = det([vec1(i,:); vec2(i,:)]);
    disp(num2str(r-r1))
    if r>0
        plot(p(i,1),p(i,2),'*'); % salient points (tu dian)       
    elseif r<0
        plot(p(i,1),p(i,2),'o'); % trough points (ao dian)         
    end
    text(p(i,1)+0.003,p(i,2)-0.01,num2str(i));
end
plot(p(:,1),p(:,2));
p_=circshift(p,1);
plot(p_(:,1),p_(:,2));
hold off;



figure 
hold on;
for i = 1:n
    r2 = det([vec1(i,:); vec2(i,:)]);
    if r2>0
        plot(p(i,1),p(i,2),'x'); % convex-points (tu dian)       
    elseif r2<0
        plot(p(i,1),p(i,2),'s'); % concave-points (ao dian)         
    end
end % for 
plot(p(:,1),p(:,2));
p=circshift(p,1);
plot(p(:,1),p(:,2));

%% 
p2 = p;
[in on] = inpolygon(p2(:,1),p2(:,2),p2(:,1),p2(:,2));


end % function







function p=createSimplyPoly(p)
    cen=mean(p);
    ang=atan2(p(:,1)-cen(1),p(:,2)-cen(2)); %每个点到坐标中心极角

    p=[p,ang];
    p=sortrows(p,3);    %按极角排序

    p=p(:,1:2);
end