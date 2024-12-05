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

for k=scandim
    for i=1:u.m(k)
        u.F{k}{i}=sum(u.F{k}{i},options{:});        
    end
end
