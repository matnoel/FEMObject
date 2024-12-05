function M = reduce(x,mu)
% function M = reduce(A,mu)

% t = A.dim2ind(mu);
% A = change_root(A,t);
% t = A.dim2ind(mu);
% U = incomplete_full(A,t);
% U = U{1}(:);
%
% M = U(1)*A.U{t}(:,1);
% for m=2:size(U,1);
%     M = M + U(m) * A.U{t}(:,m);
% end

p = getparents(x,x.dim2ind(mu));

U = cell(1, x.nr_nodes);

x_is_leaf = x.is_leaf;

% Loop through all nodes, starting with the leaf nodes
% except on parents
for ii=x.nr_nodes:-1:1
    if all(p~=ii)
        if(x_is_leaf(ii))
            % Matrix U already known
            U{ii} = x.U{ii};
        else
            % Find child nodes
            ii_left  = x.children(ii, 1);
            ii_right = x.children(ii, 2);
            
%             BUU = ttm(x.B{ii}, {U{ii_left}, U{ii_right}}, [1 2]);
%             U{ii} = matricize(BUU, [1 2], 3, false);
           
            U{ii} = U{ii_left}*reshape(x.B{ii},size(x.B{ii},1),size(x.B{ii},2)*size(x.B{ii},3));
            U{ii} = U{ii_right}*reshape(U{ii},size(x.B{ii},2),size(x.B{ii},3));

            % Clear variables to save memory
%            clear BUU;
            U{ii_left} = [];
            U{ii_right} = [];
        end
    end
end

% % Contract the parents
% p = fliplr(p);
% 
% for i=1:numel(p)-1
%     ii = p(i);
%     
%     ii_left = x.children(ii,1);
%     ii_right = x.children(ii,2);
%     
%     if p(i+1) == ii_left
%         U{ii} = ttm(x.B{ii},U{ii_right},2);
%         U{ii_right} = [];
%     else
%         U{ii} = ttm(x.B{ii},U{ii_left},1);
%         U{ii_left} = [];
%     end
% end
% 
% % Get back down to the leaf
% U = U(p);
% 
% for i=1:numel(p)-2
%     U{i+1} = ttm(U{i+1},U{i}(:)',3);
%     U{i} = [];
% end
% 
% M = x.U{x.dim2ind(mu)}*U{i+1}(:);


p = fliplr(p);
p(end) = [];
np = numel(p);

V = U;
V(p) = x.B(p);

reshape_for_prod = @(Y) reshape(Y,size(Y,1)*size(Y,2),size(Y,3));
reshape_again = @(Y,t) reshape(Y,size(x.B{t},1),size(x.B{t},2));

for i = 1:np-1
    ii = p(i);
    
    ii_left = x.children(ii,1);
    ii_right = x.children(ii,2);
    
    if any(p==ii_left)
        V{ii} = V{ii}*V{ii_right}(:);
        V{ii_left} = reshape_for_prod(V{ii_left})*V{ii}(:);
        V{ii_left} = reshape_again(V{ii_left},ii_left);
        V{ii_right} = [];
    else
        V{ii} = V{ii_left}(:)'*V{ii};
        V{ii_right} = reshape_for_prod(V{ii_right})*V{ii}(:);
        V{ii_right} = reshape_again(V{ii_right},ii_right);
        V{ii_left} = [];
    end
    V{ii} = [];
end

p = p(end);
p_left = x.children(p,1);
p_right = x.children(p,2);

p_sp = x.dim2ind(mu);

if p_sp == p_left
    V{p} = V{p}*V{p_right}(:);
else
    V{p} = V{p_left}(:)'*V{p};
end
M = x.U{p_sp}*V{p}(:);




