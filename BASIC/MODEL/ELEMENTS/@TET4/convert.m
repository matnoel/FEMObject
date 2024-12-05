function [elem,node] = convert(elem,node,elemtype,varargin)
% function [elem,node] = convert(elem,node,elemtype,varargin)

switch elemtype
    case 'TET10'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        connec = getconnec(elem);
        xnode = getcoord(node);
        
        xnodeplus{1} = (xnode(connec(:,1),:)+xnode(connec(:,2),:))/2;
        xnodeplus{2} = (xnode(connec(:,2),:)+xnode(connec(:,3),:))/2;
        xnodeplus{3} = (xnode(connec(:,3),:)+xnode(connec(:,1),:))/2;
        xnodeplus{4} = (xnode(connec(:,1),:)+xnode(connec(:,4),:))/2;
        xnodeplus{5} = (xnode(connec(:,4),:)+xnode(connec(:,2),:))/2;
        xnodeplus{6} = (xnode(connec(:,4),:)+xnode(connec(:,3),:))/2;
        for k=1:6
            num = getnumber(node);
            connecplus = max(num) + (1:getnbelem(elem))';
            connec = [connec,connecplus];
            nodeplus = NODE(POINT(xnodeplus{k}),connecplus);
            node = addnode(node,nodeplus);
        end
        elem = TET10(getnode(node,unique(connec)),getnumber(elem),...
            connec,'material',mat,'option',option);
end

