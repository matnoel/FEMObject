function A = calc_ximasse(A,varargin)

for i=1:length(A.L)
   A.L{i} = calc_ximasse(A.L{i},varargin{:}); 
end


