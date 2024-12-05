function v = getcoord(u,listenode)
% function v = getcoord(u,listenode)

switch nargin
    case 1
        v = getcoord(u.POINT);
        v = double(permute(v,[3 2 4 1]));
        
    case 2
        if isa(listenode,'ELEMENTGEOM')
            connec = getconnec(listenode);
            pos = getpos(u,connec');
        else
            if any(size(listenode)==1)
                listenode = listenode(:);
            end
            pos = getpos(u,listenode);
        end
        v = getcoord(u.POINT(pos));
        v = permute(v,[3 2 4 1]);
        
end
