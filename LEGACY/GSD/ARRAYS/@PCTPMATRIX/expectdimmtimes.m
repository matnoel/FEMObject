function a = expectdimmtimes(dim,a,varargin)

nodim = setdiff(1:getnbdim(a),dim);
a = expectnodimmtimes(nodim,a,varargin{:});

