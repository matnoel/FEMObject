function [elem,node] = convert(elem,node,elemtype,varargin)
% function [elem,node] = convert(elem,node,elemtype,varargin)

switch elemtype
    case 'SEG2AXI'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = SEG2AXI(getnode(node,unique(connec)),getnumber(elem),...
            getconnec(elem),'material',mat,'option',option);
        
    case 'SEG3'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        xnode = getcoord(node);
        
        xnodeplus = (xnode(connec(:,1),:)+xnode(connec(:,2),:))/2;
        num = getnumber(node);
        connecplus = max(num) + (1:getnbelem(elem))';
        connec = [connec,connecplus];
        nodeplus = NODE(POINT(xnodeplus),connecplus);
        node = addnode(node,nodeplus);
        elem = SEG3(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
    case 'BEAM'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        if ischarin('param',varargin)
            elem = BEAM(getnode(node,unique(connec)),getnumber(elem),...
                connec,'material',mat,'option',option,'param',getcharin('param',varargin));
        else
            elem = BEAM(getnode(node,unique(connec)),getnumber(elem),...
                connec,'material',mat,'option',option);
        end
        
    case 'BARR'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        elem = BARR(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
        
end

