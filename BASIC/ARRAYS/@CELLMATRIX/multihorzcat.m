function A=multicat(varargin)

for k=1:length(varargin)
  s=  varargin{k}.s;
  A{k}=varargin{k}.value;  
end

A = [A{:}];

A = MULTIMATRIX(A,s);

