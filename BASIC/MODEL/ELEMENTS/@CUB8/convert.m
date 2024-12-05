function [elem,node] = convert(elem,node,elemtype,varargin)
% function [elem,node] = convert(elem,node,elemtype,varargin)

switch elemtype
    case 'TET4'
        mat = getclassin('MATERIAL',varargin,getmaterial(elem));
        option = getcharin('option',varargin,getoption(elem));
        conneci = getconnec(elem);
        % Subdivision in 5 tetrahedra
%         connec = zeros(5*size(conneci,1),4);
%         connec(1:5:end,:) = conneci(:,[1,2,4,5]);
%         connec(2:5:end,:) = conneci(:,[3,4,2,7]);
%         connec(3:5:end,:) = conneci(:,[6,7,5,2]);
%         connec(4:5:end,:) = conneci(:,[8,7,5,4]);
%         connec(5:5:end,:) = conneci(:,[2,4,5,7]);
%         elem = TET4(getnode(node,unique(connec)),repmat(getnumber(elem),5,1),...
%             connec,'material',mat,'option',option);
        % Subdivision in 6 tetrahedra
        connec = zeros(6*size(conneci,1),4);
        connec(1:6:end,:) = conneci(:,[1,2,4,5]);
        connec(2:6:end,:) = conneci(:,[2,3,4,6]);
        connec(3:6:end,:) = conneci(:,[2,4,5,6]);
        connec(4:6:end,:) = conneci(:,[6,4,8,5]);
        connec(5:6:end,:) = conneci(:,[6,3,4,8]);
        connec(6:6:end,:) = conneci(:,[6,3,7,8]);
        elem = TET4(getnode(node,unique(connec)),repmat(getnumber(elem),6,1),...
            connec,'material',mat,'option',option);
        
    case 'TET10'
        elem = convert(elem,node,'TET4',varargin{:});
        [elem,node] = convert(elem,node,'TET10',varargin{:});
        
end

