function a = TENSORPRODUCTFUN(V,s,varargin)

if nargin==0
    a.m = 0;
    a.phi = cell(1,a.m);
    a.factor = 0;
    a = class(a,'TENSORPRODUCTFUN',PRODUCTSPACE());

else

    V = PRODUCTSPACE(V);    

    if getnbdim(V)~=length(varargin)
        error('le nombre de fonctions ne correspond pas avec le nombre d''espaces ')
    end

    a.m = length(varargin);
    a.phi = cell(1,a.m);
    a.factor = s;
    for i=1:a.m
        a.phi{i} = varargin{i};
    end
    a = class(a,'TENSORPRODUCTFUN',V);

end
