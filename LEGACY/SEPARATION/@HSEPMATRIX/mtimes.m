function w = mtimes(u,v,dim,varargin)

if isa(u,'double')
    if all(size(u)==1)  % Multiplication de v par un scalaire u
        w=v;
        w.alpha = v.alpha*u;

    else % On multiplie v par la matrice u sur la (ou les) dimensions dim 
        error('non programme')
        
    end
elseif isa(v,'double') 
    if all(size(v)==1)  % Multiplication de u par un scalaire v
        w=u;
        w.alpha = u.alpha*v;
        

    else % On multiplie v par la matrice u sur la (ou les) dimensions dim 
        error('non programme')
    end

else %     (H)SM * (H)SM
    w=HSEPMATRIX(u.tree);
    w.m = u.m*v.m; 
    w.F=cell(w.m,u.dim);
    w.alpha = zeros(1,w.m);
    if nargin==2
        dim = 1:u.dim;
    end
%     for i=1:u.m
%         for j=1:v.m
%             I=(i-1)*v.m+j;
%             for kkk=1:length(dim)
%                 k=dim(kkk);
%                 w.F{I,k} = u.F{i,k}*v.F{j,k};
%             end
%             w.alpha(I) = u.alpha(i)*v.alpha(j);
%         end
%     end
    [I,J] = ind2sub([u.m,v.m],1:w.m);  
    w.F(:,dim) = cellfun(@mtimes,u.F(I,dim),v.F(J,dim),'UniformOutput',0);
    w.alpha = u.alpha(I).*v.alpha(J);
    
    
end