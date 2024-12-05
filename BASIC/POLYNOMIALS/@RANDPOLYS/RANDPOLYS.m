function polys = RANDPOLYS(varargin)
% function polys = RANDPOLYS(varargin)

if nargin==0
    polys.h=cell(1,0);    
    polys.M = 0;
    polys = class(polys,'RANDPOLYS');    
elseif nargin==2 && isa(varargin{1},'RANDPOLY') && isa(varargin{2},'double')
    h=cell(1,varargin{2});
    for i=1:varargin{2}
        h{i}=varargin{1};
    end
    polys = RANDPOLYS(h{:});

elseif nargin==1 && isa(varargin{1},'RANDPOLYS')
    polys = varargin{1};
else    
    polys.h=cell(1,0);
    for k=1:nargin
        if isa(varargin{k},'RANDPOLY')
            polys.h = [polys.h , {varargin{k}}];
        elseif isa(varargin{k},'RANDPOLYS')
            polys.h = [polys.h , varargin{k}.h];    
        elseif isa(varargin{k},'RANDVARS') 
            for j=1:length(varargin{k})
                polys.h = [polys.h , {RANDPOLY(varargin{k}{j})}];  
            end
        elseif isa(varargin{k},'RANDVAR')
            polys.h = [polys.h , {RANDPOLY(varargin{k})}];  
        elseif isa(varargin{k},'cell')
            polysk = RANDPOLYS(varargin{k}{:}) ;
            polys.h = [polys.h , polysk.h]; 
        else
            error([' creation d''un randpolys a partir d''un ' class(varargin{k})  ' impossible'])
        end
    end
    polys.M = length(polys.h);
    polys = class(polys,'RANDPOLYS');
    polys = setnumber(polys);
    polys = unique(polys);
end
