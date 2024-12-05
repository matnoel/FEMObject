function rep = cmpapprox(T1,T2,param)
% function rep = cmpapprox(T1,T2,param)

if nargin==2
    rep = cmpapprox(T1,T2,'p') & cmpapprox(T1,T2,'type');
elseif nargin==3
    switch param
        case 'p'
            rep = getapproxparam(T1,'p') == getapproxparam(T2,'p') ;
            
        case 'type'
            type1 = getapproxparam(T1,'type');
            type2 = getapproxparam(T2,'type');
            if strcmp(type1,'CG')
                type1 = 'default';
            end
            if strcmp(type2,'CG')
                type2 = 'default';
            end
            rep = strcmp(type1,type2);
            
        otherwise
            rep = getapproxparam(T1,param) ==  getapproxparam(T2,param);
    end
    
end
