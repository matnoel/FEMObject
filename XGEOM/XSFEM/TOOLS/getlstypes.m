function s = getlstypes(i)

if nargin==0
    s = {'indomain','in','out','cut','touchcut','bicut','touchbicut'};

elseif ~isa(i,'double')
    error('argument doit etre un double');
else
    s = {'indomain','in','out','cut','touchcut','bicut','touchbicut'};
    if i==0
        s = s(1:7);  
    else
        s = s(i);
        if length(s)==1
            s=s{1};
        end
    end
end

