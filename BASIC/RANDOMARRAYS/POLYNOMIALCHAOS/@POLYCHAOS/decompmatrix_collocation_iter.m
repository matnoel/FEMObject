function varargout = decompmatrix_collocation_iter(apc,fun,varargin)
% function u = decompmatrix_collocation_iter(PC,fun,varargin)
% Decomposition of a function of random variables on polynomial chaos basis
% by Lagrange interpolation (or collocation)
% PC : POLYCHOAS or PCMODEL or FESTOMODEL or PCMATRIX or PCARRAY
% fun : function_handle of random variables RANDVARS(PC) associated to PC
% varargin : arguments of function fun
% call yk = fun(x(k,:),varargin{:}) for each sample vector x(k,:)
% x(:,i) is the sample vector containing all the evaluations of random variable i
% x(k,:) is the sample vector containing the evaluation k of all random variables
% yk is the response vector contaning the values returned by function fun applied to sample vector x(k,:)
% if PC is a POLYCHAOS, x is the sample matrix containing the evaluations of random variables associated to POLYCHAOS
%          a FESTOMODEL or a PCMODEL, x is the sample matrix containing the evaluations of random variables associated to PCMODEL
%          a PCMATRIX or a PCARRAY, x is the sample matrix containing the evaluations of random variables associated to PCMATRIX

if ischarin('display',varargin)
    display_ = 1;
    varargin = delonlycharin('display',varargin);
else
    display_ = 0;
end

n = nargout;
fun = fcnchk(fun);
PC = getPC(apc);

xpc = calc_collocationpoints(PC);

N = size(x,1);
if display_
    fprintf('\nNumber of collocation points = %d\n',N);
end

x = transfer(apc,xpc);

% t0 = clock;
% for k=1:N
%     if display_
%         pourcentage(k,N)
%     end
%     yk = cell(1,n);
%     [yk{:}] = fun(x(k,:),varargin{:});
%     for i=1:n
%         if isa(yk{i},'cell')
%             for l=1:length(yk{i})
%                 if k==1
%                     s{i}{l} = size(yk{i}{l});
%                 end
%                 u{i}{l}(:,k) = yk{i}{l}(:);
%             end
%         else
%             if k==1
%                 s{i} = size(yk{i});
%             end
%             u{i}(:,k)= yk{i}(:);
%         end
%         
%     end
%     if display_ && k==3
%         fprintf('\nNumber of collocation points = %d -- Estimated time : %4.1f s\n',N,etime(clock,t0)*N/k);
%     end
% end

[u,s] = calc_response_iter(fun,x,n,varargin{:});

for i=1:n
    if isa(u{i},'cell')
        for l=1:length(u{i})
            u{i}{l} = PCMATRIX(u{i}{l},s{i}{l},PC);
        end
    else
        u{i} = PCMATRIX(u{i},s{i},PC);
    end
end

varargout = u;
