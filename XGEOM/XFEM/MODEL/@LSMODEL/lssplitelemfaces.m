function S=lssplitelemfaces(S,varargin)

for i=1:getnbfacets(S)
S = setfacet(S,i,lssplitelem(getfacet(S,i),varargin{:}));
end
for i=1:getnbridges(S)
S = setridge(S,i,lssplitelem(getridge(S,i),varargin{:})); 
end
