function [x,err_min,delta] = select_opt_path(y,A,p,varargin)
% function [x,err_min,delta] = select_opt_path(y,A,p)
% function [x,err_min,delta] = select_opt_path(y,A,p,'leaveout')
% Select optimal solution x from path p using relative leave-one-out cross-validation error
% 
% function [x,err_min,delta] = select_opt_path(y,A,p,'kfold',k)
% Select optimal solution x from path p using relative k-fold cross-validation error
% 
% function [x,err_min,delta] = select_opt_path(y,A,p,'leaveout','correction')
% Select optimal solution x from path p using corrected relative leave-one-out cross-validation error
% 
% function [x,err_min,delta] = select_opt_path(y,A,p,'kfold',k,'correction')
% Select optimal solution x from path p using corrected relative k-fold cross-validation error

if ischarin('parallel',varargin)
    parallel = 1;
    varargin = delonlycharin('parallel',varargin);
else
    parallel = 0;
end

pattern = (p~=0);
rep = all(~pattern,1);
pattern(:,rep) = [];
p(:,rep) = [];

[~,I,~] = unique(pattern','rows','first');
I = sort(I);
p = p(:,I);
pattern = pattern(:,I);

n = size(pattern,2);
x = cell(1,n);
x(:) = {zeros(size(pattern,1),1)};
delta = cell(1,n);
err = zeros(1,n);

if parallel
    parfor i=1:n
        % indices of non-zero coefficents in p(:,i)
        % ind = find(p(:,i));
        ind = pattern(:,i);
        % compute reduced matrix
        A_red = A(:,ind);
        % recompute non-zero coefficients of p(:,i)
        % using ordinary least-squares minimization and store them in x_red
        x_red = solve(A_red'*A_red,A_red'*y);
        % x_red = regress(y,A_red);
        % compute relative cross-validation error
        [err(i),delta{i}] = calc_fast_cv_error(y,A_red,x_red,varargin{:});
        % fun = @(y,A_red) solve(A_red'*A_red,A_red'*y);
        % [err(i),delta{i}] = calc_cv_error(y,A_red,fun,varargin{:});
        x{i}(ind) = x_red;
    end
else
    for i=1:n
        % find indices of non-zero coefficents in p(:,i)
        % ind = find(p(:,i));
        ind = pattern(:,i);
        % compute reduced matrix
        A_red = A(:,ind);
        % recompute non-zero coefficients of p(:,i)
        % using ordinary least-squares minimization and store them in p_red
        x_red = solve(A_red'*A_red,A_red'*y);
        % x_red = regress(y,A_red);
        % compute relative cross-validation error
        [err(i),delta{i}] = calc_fast_cv_error(y,A_red,x_red,varargin{:});
        % fun = @(y,A_red) solve(A_red'*A_red,A_red'*y);
        % [err(i),delta{i}] = calc_cv_error(y,A_red,fun,varargin{:});
        x{i}(ind) = x_red;
    end
end
x_path = [x{:}];
delta_path = [delta{:}];

[err_min,ind] = min(err);
% [err_max, ind_max] = max(err);
% complexity = sum(pattern,1);
% indclose = find(log10(err/err_min)<1);
% [~,ind2] = min(complexity(indclose));
% ind = indclose(ind2);

x = x_path(:,ind);
delta = delta_path(:,ind);

end

