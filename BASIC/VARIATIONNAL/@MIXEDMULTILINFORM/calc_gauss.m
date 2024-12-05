function gauss = calc_gauss(a,varargin)

order = varargin{end};
if nargin==getn(a)+1
    order=0;
    for i=1:length(a.p)
    o=0;
    for k=1:getn(a)
        o = max(o,getorder(a.p(i),varargin{k}));
    end
    order = order + o;
    end    
    order = order + getorder(a.pk,varargin{1});
end
    gauss = calc_gauss(varargin{1},order);
    
function o = getorder(p,elem)
if isempty(p)
    o=0;
else
switch p
    case 0
       o = orderN(elem);
    case 1
       o = orderDN(elem); 
end
end


return


