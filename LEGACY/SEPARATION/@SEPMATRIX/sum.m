function u = sum(u,varargin)

if nargin==1
   options={};
else
   options = varargin(1);
end
if nargin<=2
   scandim = 1:u.dim;
else
   scandim = varargin{2}; 
end
for i=1:u.m
    for k=scandim
u.F{i,k}=sum(u.F{i,k},options{:});        
    end
end
