function [elem,node] = convert(elem,node,elemtype,varargin)
% function [elem,node] = convert(elem,node,elemtype,varargin)

switch elemtype
    case 'TRI3'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        conneci = getconnec(elem);
        connec = zeros(2*size(conneci,1),3);
        connec(1:2:end,:) = conneci(:,[1,2,3]);
        connec(2:2:end,:) = conneci(:,[1,3,4]);
        elem = TRI3(getnode(node,unique(connec)),[getnumber(elem);getnumber(elem)],...
            connec,'material',mat,'option',option);
        
    case 'COQ4'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = COQ4(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        if ischarin('fullintegration',varargin)
            elem = setparam(elem,'fullintegration',1);
        end
        
    case 'DKQ'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = DKQ(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
    case 'DSQ'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = DSQ(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
    case 'DKT'
        elem = convert(elem,node,'TRI3',varargin{:});
        elem = convert(elem,node,'DKT',varargin{:});
        
    case 'DST'
        elem = convert(elem,node,'TRI3',varargin{:});
        elem = convert(elem,node,'DST',varargin{:});
        
    case 'TRI6'
        elem = convert(elem,node,'TRI3',varargin{:});
        [elem,node] = convert(elem,node,'TRI6',varargin{:});
        
end

