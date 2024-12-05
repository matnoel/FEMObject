function value = getelementfield(M,field,g)

if nargin==2
    g= 1:M.nbgroupelem;
end
value=cell(1,length(g));
for i=1:length(g)
    p = g(i);
    value{i} = getelementfield(M.groupelem{p},field);
end

