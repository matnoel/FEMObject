function A=multivertcat(varargin)

isadouble = [] ;
for  k=1:length(varargin)
    isadouble=[isadouble, isa(varargin{k}.value,'double')];
end

if all(isadouble)
alllines = 1;    
for k=1:length(varargin)
alllines = alllines & varargin{k}.sm(1)==1;
end

if alllines
sm=zeros(1,2);
sm(2)=1;
for k=1:length(varargin)
  s=  varargin{k}.s;
  A{k}=varargin{k}.value;  
  sm(1) = sm(1)+varargin{k}.sm(1);
end
A = [A{:}];
A = MULTIMATRIX(A,s);
else
    
for k=1:length(varargin)
varargin{k}=multitranspose(varargin{k});
end  
A = multihorzcat(varargin{:});
A = multitranspose(A);

end

elseif  all(isadouble==0)

    
A = varargin{1};
for k=2:length(varargin)
A.value = [A.value;varargin{k}.value];    
end
A.sm = size(A.value);    
    
else

for k=find(isadouble)
    varargin{k}=mat2cell(varargin{k});
end
A = multivertcat(varargin{:});

end


