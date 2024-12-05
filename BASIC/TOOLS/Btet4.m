function B = Btet4(xnode)
% B = Btet4(xnode)

DN = DNtet4(xnode);
%
B=zeros(6,12);
B(1,1:3:end)=DN(1,:);
B(2,2:3:end)=DN(2,:);
B(3,3:3:end)=DN(3,:);
B(4,1:3:end)=DN(2,:);
B(4,2:3:end)=DN(1,:);
B(5,1:3:end)=DN(3,:);
B(5,3:3:end)=DN(1,:);
B(6,2:3:end)=DN(3,:);
B(6,3:3:end)=DN(2,:);


