function [tevol,rep] = gettevol(L)
% function [tevol,rep] = gettevol(L)

switch getapproxparam(L,'type')
    case 'default'
        tevol = gett(L);
        rep = 1:length(tevol);
    case 'DG'
        
        tsol = gettapprox(L);
        p = getapproxparam(L,'p') ;
        
        switch p
            case 0
                tevol = tsol;
                rep = 1:length(tevol);
            otherwise
                rep = [1:(p+1):getnt(L)*(p+1)];
                %for j=2:p
                %rep=[rep;j:(p+1):getnt(L)*(p+1)];
                %end
                rep = [rep,getnt(L)*(p+1)];
                rep = rep(:)';
                tevol = tsol(rep);
        end
        
end
