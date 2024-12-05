function c = getlstypes(i)

if nargin==0
    c = {'g','y','w','r','m','b','c','k','y','g'};
elseif ~isa(i,'double')
    error('argument doit etre un double');
else
    c = {'g','y','w','r','m','b','c','k','y','g'};
    if i==0
        c = c(1:7);  
    else
        c = c(i);
        if length(s)==1
            c=c{1};
        end
    end
end


