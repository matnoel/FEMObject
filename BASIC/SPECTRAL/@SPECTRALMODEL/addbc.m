function B=addbc(B,type,varargin)

switch type
    case 'dirichlet'
        setfemobjectoptions('tolerancepoint',1e-10);
        bc.type = type;
        switch varargin{1}
            case 'all'
                if B.dim==1
                    bc.node = [1,getnbnode(B)];
                elseif B.dim==2
                    [X,Y] = meshgrid(1:getnbpoints(B.L),1:getnbpoints(B.L));
                    bc.node = find(X==1 | X==getnbpoints(B.L) | Y==1 | Y==getnbpoints(B.L));
                end
            case 'right'
                if B.dim==1
                    bc.node = getnbnode(B);
                elseif B.dim==2
                    [X,Y] = meshgrid(1:getnbpoints(B.L),1:getnbpoints(B.L));
                    bc.node = find(X==getnbpoints(B.L));
                end
            case 'left'
                if B.dim==1
                    bc.node = 1;
                elseif B.dim==2
                    [X,Y] = meshgrid(1:getnbpoints(B.L),1:getnbpoints(B.L));
                    bc.node = find(X==1);
                end
            case 'up'
                if B.dim==2
                    [X,Y] = meshgrid(1:getnbpoints(B.L),1:getnbpoints(B.L));
                    bc.node = find(Y==getnbpoints(B.L));
                end
            case 'bottom'
                if B.dim==2
                    [X,Y] = meshgrid(1:getnbpoints(B.L),1:getnbpoints(B.L));
                    bc.node = find(Y==1);
                end
            case 'none'
                bc.node=1;
            otherwise
                error('mauvaise entree')
        end
    case 'periodic'
        bc.type  = type;
end

if ~isempty(B.bc)
    B.bc.node = union(B.bc.node,bc.node);
else
    B.bc = bc;
end

