function u = freematrix(u,S)
% function u = freematrix(u,S)

if iscell(u.V)
    V = getvalue(u.V);
    for i=1:length(V)
        V{i} = freematrix(S,V{i});
    end
    u.V = setvalue(u.V,V);
    u.V = setsize(u.V,size(V{1}));
else
    u.V = freematrix(S,u.V);
end
