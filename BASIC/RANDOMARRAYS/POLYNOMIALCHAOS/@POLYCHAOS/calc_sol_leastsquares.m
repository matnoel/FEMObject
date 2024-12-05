function [u,err] = calc_sol_leastsquares(D,y,varargin)
% function [u,err] = calc_sol_leastsquares(D,y,varargin)
% Solve (regularized) least-squares minimization problem
% and compute cross-validation error
% D : matrix containing the evaluations of basis functions
% y : response vector or matrix

n = length(y);
u = cell(1,n);
err = cell(1,n);
for i=1:n
    if isa(y{i},'cell')
        for l=1:length(y{i})
            [u{i}{l},err{i}{l}] = solve_leastsquares(D,y{i}{l},varargin{:});
        end
    else
        [u{i},err{i}] = solve_leastsquares(D,y{i},varargin{:});
    end
end

end
