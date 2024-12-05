function e = getenrichtypetip(c,k)
% function e = getenrichtypetip(c,k)

if nargin==1
    k=1:getnbtip(c);
end
e = zeros(1,length(k));
for i=1:length(k)    
e(i) = getenrichtype(getlstip(c,k(i)));
end



