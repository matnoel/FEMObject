function e = isenrich(c,choix,k)
if nargin==1
e=0;
for k=1:length(c.LEVELSETS)
e = e | isenrich(c.LEVELSETS{k});
end
else
    switch choix
        case 'tip'
if nargin==3            
    e = isenrich(c.LEVELSETS{1+k});
else
    e=0;
for k=2:length(c.LEVELSETS)
e = e | isenrich(c.LEVELSETS{k});
end 
end
        case 'support'
            e = isenrich(c.LEVELSETS{1});
    end
        
    
end

