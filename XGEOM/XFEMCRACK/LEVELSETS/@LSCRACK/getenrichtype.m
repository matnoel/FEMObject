function e = getenrichtype(c,choix,k)

switch choix
    case 'support'
        e = getenrichtype(c.LEVELSETS{1});
    case 'tip'
        if nargin==3
        e = (getenrichtype(c.LEVELSETS{1+k}));        
        else
            e=[];
            for k=2:length(c.LEVELSETS)
        e = [e,(getenrichtype(c.LEVELSETS{k}))];            
            end
        end
        
    otherwise
        error('pas defini')
end
