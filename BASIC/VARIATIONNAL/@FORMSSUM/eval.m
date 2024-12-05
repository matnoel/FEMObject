function ae = eval(a,S,varargin)

ae = eval(a.forms{1},S,varargin{:});
for i=2:length(a.forms)
 ae = ae + eval(a.forms{i},S,varargin{:});  
end

