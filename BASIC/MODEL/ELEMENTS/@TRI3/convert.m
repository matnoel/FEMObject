function [elem,node] = convert(elem,node,elemtype,varargin)
% function [elem,node] = convert(elem,node,elemtype,varargin)

switch elemtype
    case 'DKT'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = DKT(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
    case 'DST'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = DST(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
    case 'DKQ'
        warning(['Change elem type ' elemtype ' to DKT'])
        elem = convert(elem,node,'DKT',varargin{:});
        
    case {'DSQ','COQ4'}
        warning(['Change elem type ' elemtype ' to DST'])
        elem = convert(elem,node,'DST',varargin{:});
        
    case 'TRI3AXI'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = TRI3AXI(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
    case 'TRI6'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        xnode = getcoord(node);
        
        xnodeplus{1} = (xnode(connec(:,1),:)+xnode(connec(:,2),:))/2;
        xnodeplus{2} = (xnode(connec(:,2),:)+xnode(connec(:,3),:))/2;
        xnodeplus{3} = (xnode(connec(:,3),:)+xnode(connec(:,1),:))/2;
        for k=1:3
            num = getnumber(node);
            connecplus = max(num) + (1:getnbelem(elem))';
            connec = [connec,connecplus];
            nodeplus = NODE(POINT(xnodeplus{k}),connecplus);
            node = addnode(node,nodeplus);
        end
        elem = TRI6(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
end
