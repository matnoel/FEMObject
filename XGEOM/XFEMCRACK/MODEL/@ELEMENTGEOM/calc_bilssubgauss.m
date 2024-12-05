function varargout = calc_bilssubgauss(elem,ls1,ls2,order,varargin)
% function varargout = calc_bilssubgauss(elem,ls1,ls2,order,varargin)

p = calc_gaussorder(elem,order);

varargout=cell(1,4);
[varargout{:}] = elem_bilssubgauss(elem,ls1,ls2,p,varargin{:});
for i=1:length(varargout)
varargout{i} = permutegaussND(varargout{i});
end


