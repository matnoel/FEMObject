function w=assembleblock(u)

s=u.s;
sm=u.sm;

ip=[];
jp=[];
vp=[];
for i=1:sm(1)
for j=1:sm(2)       
[ii,jj,v] = find(reshape(u.value(:,(j-1)*sm(1)+i),s));
ip=[ip;(i-1)*s(1)+ii];
jp=[jp;(j-1)*s(2)+jj];
vp=[vp;v];
end
end

w=sparse(ip,jp,vp,s(1)*sm(1),s(2)*sm(2));
