function P = getvertices(C)
% function P = getvertices(C)

P = cell(1,4);
switch C.indim
    case 2
        P{1} = [C.cx-C.a,C.cy];
        P{2} = [C.cx,C.cy-C.b];
        P{3} = [C.cx+C.a,C.cy];
        P{4} = [C.cx,C.cy+C.b];
        
        v = [C.vx,C.vy];
        v = v/norm(v);
        R = [v(1) v(2);
            -v(2) v(1)];
    case 3
        P{1} = [C.cx-C.a,C.cy,C.cz];
        P{2} = [C.cx,C.cy-C.b,C.cz];
        P{3} = [C.cx+C.a,C.cy,C.cz];
        P{4} = [C.cx,C.cy+C.b,C.cz];
        
        n = [C.nx,C.ny,C.nz];
        n = n/norm(n);
        Q = [0 -n(3) n(2);
            n(3) 0 -n(1);
            -n(2) n(1) 0];
        v = [C.vx,C.vy];
        v = v/norm(v);
        R = eye(3)-v(2)*Q+(1-v(1))*Q^2;
end

for i=1:4
    P{i} = P{i}*R;
end
