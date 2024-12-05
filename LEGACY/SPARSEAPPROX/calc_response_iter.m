function [y,s] = calc_response_iter(fun,x,n,varargin)
% function [y,s] = calc_response_iter(fun,x,n,varargin)
% Compute response matrix y returned by function fun applied to each sample vector x(k,:)
% fun : function_handle of random variables
% x : sample matrix
% n : number of output arguments for fun
% varargin : arguments of function fun
% call yk = fun(x(k,:),varargin{:}) for each sample vector x(k,:)
% x(:,i) is the sample vector containing all the evaluations of random variable i
% x(k,:) is the sample vector containing the k-th evaluation of all random variables
% yk is the response vector containing the values returned by function fun applied to sample vector x(k,:)

if ischarin('display',varargin)
    varargin = delonlycharin('display',varargin);
end

if nargin<3 || isempty(n)
    n = 1;
end

N = size(x,1);
ytmp = cellfun(@(x) cell(1,n),cell(1,N),'UniformOutput',false);
parfor k=1:N
    [ytmp{k}{:}] = fun(x(k,:),varargin{:});
end

y = cell(1,n);
s = cell(1,n);
for k=1:N
    yk = ytmp{k};
    for i=1:n
        if isa(yk{i},'cell')
            for l=1:length(yk{i})
                if k==1
                    s{i}{l} = size(yk{i}{l});
                end
                y{i}{l}(:,k) = yk{i}{l}(:);
            end
        else
            if k==1
                s{i} = size(yk{i});
            end
             y{i}(:,k) = yk{i}(:);
        end
    end
end

end
