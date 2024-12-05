function s = getnbtimedof(u)
% function s = getnbtimedof(u)

p = getapproxparam(u,'p');
type = getapproxparam(u,'type');
nt = getnt(u);
switch type
    case 'default'
        if p==0
            s = nt;
        else
            s = nt*p+1 ;
        end
    case 'DG'
        s = nt*(p+1);
        
end
