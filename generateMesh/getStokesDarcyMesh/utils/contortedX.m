function mdxx = contortedX(xx, h)
%
%   limit the xx in the domain:[0,1]
%   

mdxx = zeros(length(xx),1);

%--- we let 'begin point' short as 'bp', and 'center point' short as 'cp'
max_h = (abs(max(xx)-min(xx)))/h;
bp = 0;
for ii_h = 1:max_h/2
    cp = h + 2*(ii_h-1)*h;
    theta = 0.8*exp(cp);
    k = theta*1/h;
    
    if mod(ii_h,2) == 1
        for ii_x = 1:length(xx)
            if cp-h < xx(ii_x) && xx(ii_x) <= cp+h
                mdxx(ii_x) = -k*abs(xx(ii_x)-cp) + k*h;
            end
        end % for ii_x
    else 
        for ii_x = 1:length(xx)
            if cp-h < xx(ii_x) && xx(ii_x) <= cp+h
                mdxx(ii_x) = (-k*abs(xx(ii_x)-cp) + k*h);
            end
        end % for ii_x
    end
    
end % for 

% figure
% plot(xx(1:max_h+1),mdxx(1:max_h+1))


end