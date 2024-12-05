function s=size(u,k)
switch nargin
    case 1
        s=size(u.V);
if u.stodim>length(s)
 s = [s,ones(1,u.stodim-length(s))];
end    
s(u.stodim)=[];
    case 2
s=size(u);
if k>length(s)
s=1; 
else
s = s(k);
end
end
