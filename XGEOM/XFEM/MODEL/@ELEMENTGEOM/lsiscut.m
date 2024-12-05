function iscut = lsiscut(elem,ls,varargin)

iscut = ~(lsisin(elem,ls,varargin{:}) | lsisout(elem,ls,varargin{:}));
