function u2 = transfer(u1,T2)
% function u2 = transfer(u1,T2)

T1 = gettimemodel(u1);
p1 =    getapproxparam(T1,'p');
p2 =    getapproxparam(T2,'p');
type1 = getapproxparam(T1,'type');
type2 = getapproxparam(T2,'type');


if cmpapprox(T1,T2)
    u2 = u1 ;
elseif (strcmp(type1,'default') ||  strcmp(type1,'CG')) &&  strcmp(type2,'DG')
    u2 = project(u1,T2);
else
    error('le transfert n''est pas programme')
    switch type1
        case {'default','CG'}
            if p2==0 & p1==1
                
            end
            
        case {'DG'}
            
            
    end
    
    
end