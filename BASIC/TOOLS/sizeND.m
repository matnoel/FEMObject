function s = sizeND(A,k)
% function s = sizeND(A,k)

s = size(A);
if length(s)==2
    s = [s,1] ;
end
s = s(3:end);
if nargin==2
    if length(s)<max(k)
        s = [s,ones(1,max(k)-length(s))];
    end
    s = s(k);
end
