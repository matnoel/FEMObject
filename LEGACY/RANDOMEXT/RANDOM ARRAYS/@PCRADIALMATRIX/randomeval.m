function a = randomeval(pcr,x,varargin)
% function a = randomeval(pcr,x,varargin)
% see also PCMATRIX/randomeval

s = size(pcr.V);


L = pcr.D*double(randomeval(pcr.L(:),x,varargin{:}));


if iscell(pcr.V)
    
    L = MULTIMATRIX(L,[pcr.m,1]);
    if length(L)==1
        L = double(L);
    end
    a = pcr.V{1}*L(1);
    for k=2:pcr.m
        a = a + pcr.V{k}*L(k);
    end
    
    if prod(sizem(L))==1
        a = double(a);
    end
else
    
    a = double(pcr.V)* L;
    if size(a,2)==1
        a = reshape(a,s);
    elseif all(s==1)
        a = reshape(a,1,size(a,2));
    else
        a = MULTIMATRIX(a,s);
    end
    
end

