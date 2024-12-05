function tapprox=gettapprox(L)

switch getapproxparam(L,'type')
    case 'default'
        p = getapproxparam(L,'p') ;
        switch p
            case 1
                tapprox = gett(L);
        end
        
    case 'DG'
        
        p = getapproxparam(L,'p') ;
        switch p
            case 0
                tapprox = gettm(L);
            otherwise
                t = gett(L);
                dt = getdt(L);
                tapprox = t(1:end-1);
                for j=1:p
                    tapprox=[tapprox;tapprox(j,:)+dt/(p)];
                end
                tapprox=tapprox(:)';
        end
        
end
