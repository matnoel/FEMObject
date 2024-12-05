function B = Btri3(xnode)
% B = Btri3(xnode)

DN = DNtri3(xnode);

B=zeros(3,6);
B(1,1:2:end)=DN(1,:);
B(2,2:2:end)=DN(2,:);
B(3,1:2:end)=DN(2,:);
B(3,2:2:end)=DN(1,:);

