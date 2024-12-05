function rv = eval(ma,x)
% function rv = eval(ma,x)

n=nbparam(ma.RV);
paramnames = fieldnames(ma.param);
param = struct2cell(ma.param);
rv=ma.RV;

for k=1:n
    if isa(param{k},'inline') | isa(param{k},'function_handle')
    val=param{k}(varargin{:});  
    else
    val=param{k};
    end
rv=setparam(rv,paramnames{k},val);
end
