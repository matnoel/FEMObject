function ke = rigi(mat,elem,xnode,xgauss,varargin)
% function ke = rigi(mat,elem,xnode,xgauss,varargin)

ke = diff(mat,elem,xnode,xgauss,varargin{:});

ke = ke + adv(mat,elem,xnode,xgauss,varargin{:});

if isparam(mat,'stabilize')
    
    ke = ke+stab(mat,elem,xnode,xgauss,varargin{:});
    
end

