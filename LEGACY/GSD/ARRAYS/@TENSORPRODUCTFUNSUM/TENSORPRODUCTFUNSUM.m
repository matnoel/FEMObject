function a = TENSORPRODUCTFUNSUM(varargin)

if nargin==0
    a.tensorfuns = cell(1,0);
    a = class(a,'TENSORPRODUCTFUNSUM',PRODUCTSPACE()); 
    superiorto('TENSORPRODUCTFUN');

elseif nargin==1 && isa(varargin{1},'PRODUCTSPACE')
    a.tensorfuns = cell(1,0);
    a = class(a,'TENSORPRODUCTFUNSUM',varargin{1}); 
    superiorto('TENSORPRODUCTFUN');
else
    a.tensorfuns = cell(1,length(varargin));
    V = PRODUCTSPACE();
    for i=1:length(a.tensorfuns)
        if ~isa(varargin{i},'TENSORPRODUCTFUN')
            error('RENTRER DES TENSORPRODUCTFUN en arguments')
        end
        a.tensorfuns{i} = varargin{i};
        V = union(V,getproductspace(varargin{i}));
    end

    a = class(a,'TENSORPRODUCTFUNSUM',V);
    superiorto('TENSORPRODUCTFUN');

end
