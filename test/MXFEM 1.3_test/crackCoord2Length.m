global CRACK

a    = zeros(size(CRACK,1)-1,1);
a(1) = sqrt((CRACK(2,1)-CRACK(1,1))^2+(CRACK(2,2)-CRACK(1,2))^2);
for i = 2:(size(CRACK,1)-1)
    a(i) = a(i-1)+sqrt((CRACK(i+1,1)-CRACK(i,1))^2+(CRACK(i+1,2)-CRACK(i,2))^2);
end