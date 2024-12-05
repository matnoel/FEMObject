function u = setcoord(u,coord,listenode)
% function u = setcoord(u,coord,listenode)

switch nargin
    case 1

    case 2
        u.POINT = setcoord(u.POINT,coord);
        
    case 3
        if isa(listenode,'ELEMENTGEOM')
            connec = getconnec(listenode);
            pos = getpos(u,connec');
        else
            if any(size(listenode)==1)
                listenode = listenode(:);
            end
            pos = getpos(u,listenode);
        end
        u.POINT(pos) = setcoord(u.POINT(pos),coord);
        
end
