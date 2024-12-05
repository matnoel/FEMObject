function A=multicat(varargin)

isadouble = [] ;
for  k=1:length(varargin)
    isadouble=[isadouble, isa(varargin{k}.value,'double')];
end

if all(isadouble)

for k=1:length(varargin)
  s=  varargin{k}.s;
  A{k}=varargin{k}.value;  
end

A = [A{:}];

A = MULTIMATRIX(A,s);


elseif all(isadouble==0)

    
A = varargin{1};
A.value = A.value(:);
for k=2:length(varargin)
A.value = [A.value;varargin{k}.value(:)];    
end
A.sm = size(A.value);    
    
else

for k=find(isadouble)
    varargin{k}=mat2cell(varargin{k});
end


A = multicat(varargin{:});

end


