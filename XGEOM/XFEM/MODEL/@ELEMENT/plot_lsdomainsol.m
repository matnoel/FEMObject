function varargout = plot_lsdomainsol(elem,node,q,ls,varargin)

if ~isenrich(elem)
Helem = plot_sol(elem,node,q,varargin{:});
else
if getenrichtype(ls)>1 
    error('pas programme')
else
    Helem = plot_sol(elem,node,q,varargin{:});
end

end

  if nargout==1
      varargout{1}=Helem;
  end
  
