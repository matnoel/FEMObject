function u = meanperinterval(u)

switch getapproxparam(u,'type')
     case {'DG'}
        switch  getapproxparam(u,'p')   
            case 0
                
             case 1
                 
if isa(u.value,'double')
    for i=1:length(u.TIMEMODEL)/2 
     utemp = (u.value(:,i*2-1)+u.value(:,i*2))/2 ;
     u.value(:,i*2-1)=utemp;
     u.value(:,i*2)=utemp;
    end
end
            otherwise 
                error(' pas prevu')
        end
        
end