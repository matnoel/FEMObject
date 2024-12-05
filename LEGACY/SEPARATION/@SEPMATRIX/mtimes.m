function w = mtimes(u,v,dim,varargin)

if isa(u,'double')
    if all(size(u)==1)
        w=v;
        w.alpha = v.alpha*u;
    elseif nargin==2
        error('pas prevu, inverser la multiplication')
        
    else
        w=v;
        for i=1:v.m
            for j=1:length(dim)
                k=dim(j);
                w.F{i,k} = u*v.F{i,k};
            end
        end
    end
    
elseif isa(v,'double')
    if all(size(v)==1)
        w=u;
        w.alpha = u.alpha*v;
    elseif nargin==2
        if length(size(v))~=getdim(u)
            error('')
        end
        %s = size(v);
        %w=zeros(size(v));
        s=cellfun('size',u.F(1,:),1);
        w = zeros(s);
        for i=1:u.m
            w = w + u.alpha(i)*double(ttm(tensor(v),u.F(i,:)));
%             for k=1:getdim(u)
%                 r = [k,1:k-1,k+1:getdim(u)];
%                 wtemp = permute(wtemp,r);
%                 wtemp(:,:) = u.F{i,k}*wtemp(:,:);
%                 wtemp = permute(wtemp,[2:k,1,k+1:getdim(u)]);
%             end
%            w=w+wtemp;
        end
    else
        w=u;
        for i=1:u.m
            for j=1:length(dim)
                k=dim(j);
                w.F{i,k} = u.F{i,k}*v;
            end
        end
    end
    
    
else
    w=u;
    w.dim = u.dim;
    w.m = u.m*v.m;
    w.F = cell(w.m,w.dim);
    w.alpha = zeros(1,w.m);
    if nargin==2
        dim = 1:u.dim;
    end
    global mexcompiled;
    if ~mexcompiled
     [I,J] = ind2sub([u.m,v.m],1:w.m);  
     w.F(:,dim) = cellfun(@mtimes,u.F(I,dim),v.F(J,dim),'UniformOutput',0);
     w.alpha = u.alpha(I).*v.alpha(J);
%         for i=1:u.m
%             for j=1:v.m
%                 I=(i-1)*v.m+j;
%                 for kkk=1:length(dim)
%                     k=dim(kkk);
%                     w.F{I,k} = u.F{i,k}*v.F{j,k};
%                 end
%                 w.alpha(I) = u.alpha(i)*v.alpha(j);
%             end
%         end
    else
    w.alpha=kron(u.alpha(1:u.m),v.alpha(1:v.m));
    w.F(:,dim)=multiplyF(u.F(:,dim),v.F(:,dim));
    end
    
    
end
