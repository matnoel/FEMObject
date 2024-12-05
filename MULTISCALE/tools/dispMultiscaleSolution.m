function dispMultiscaleSolution(U,w,lambda,varargin)
% function dispMultiscaleSolution(U,w,lambda,varargin)
% Display multiscale solution u=(U,w) and Lagrange multiplier lambda
% U: FunctionalBasisArray of global solution U
% w: FunctionalBasisArray of local solution w
% lambda: FunctionalBasisArray of Lagrange multiplier lambda

n = numel(w);

fprintf('spatial dimension = %d for U\n',U.sz(1))
for k=1:n
    fprintf('                  = %d for w{%u}\n',w{k}.sz(1),k)
    fprintf('                  = %d for lambda{%u}\n',lambda{k}.sz(1),k)
end
if length(U.sz)==2
    fprintf('time dimension = %d for U\n',U.sz(2))
    for k=1:n
        fprintf('               = %d for w{%u}\n',w{k}.sz(2),k)
        fprintf('               = %d for lambda{%u}\n',lambda{k}.sz(2),k)
    end
end
fprintf('parametric dimension = %d\n',ndims(U.basis))
fprintf('basis dimension = %d for U\n',cardinal(U.basis))
for k=1:n
    fprintf('                = %d for w{%u}\n',cardinal(w{k}.basis),k)
    fprintf('                = %d for lambda{%u}\n',cardinal(lambda{k}.basis),k)
end
fprintf('order = [ %s ] for U\n',num2str(max(U.basis.indices.array)))
for k=1:n
    fprintf('      = [ %s ] for w{%u}\n',num2str(max(w{k}.basis.indices.array)),k)
    fprintf('      = [ %s ] for lambda{%u}\n',num2str(max(lambda{k}.basis.indices.array)),k)
end
% fprintf('multi-index set for U = \n')
% disp(num2str(U.basis.indices.array))
% for k=1:n
%     fprintf('multi-index set for w{%u} = \n',k)
%     disp(num2str(w{k}.basis.indices.array))
%     fprintf('multi-index set for lambda{%u} = \n',k)
%     disp(num2str(lambda{k}.basis.indices.array))
% end

end
