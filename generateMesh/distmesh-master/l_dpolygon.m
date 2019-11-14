function [ dist ] = l_dpolygon( p, v )

n_p = size(p,1);
n_s = size(v,1)-1;

dist = zeros(n_s,n_p);
for i_s=1:n_s
  v_i = v([i_s,i_s+1],:);
  n_p = size(p,1);
  w  = v_i(2,:)-v_i(1,:);
  ix1 = ones(n_p,1);
  vp = v_i(ix1,:)-p;
  w1 = w(ix1,:);
  s = dot(w1,vp,2);
  u = -s/(w*w.');
  u(u<0) = 0;
  u(u>1) = 1;
  h = w1.*[u,u]+vp;
  dist(i_s,:) = sqrt(dot(h,h,2));
end
dist = (-1).^(inpolygon(p(:,1),p(:,2),v(:,1),v(:,2))).*min(dist).';