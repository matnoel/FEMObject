function [u,result,W,Lambda] = solve_gsd_arnoldi(A,b,solver,solverlambda,solverlambdaupdate)

dim = getdim(A);
param = getparam(solver);

n = size(b);
bu = b;

u = SEPMATRIX(dim);
u0=u;
W = zeros(n(1),0);
Lambda = SEPMATRIX(dim);

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);
result.resultlambda = {};
switch param.errorindicator
    case 'reference'
normref = normND(param.reference);                
    case 'residual'
normref = norm(b);                
end

for iii=0:param.restart

if size(W,2)==param.maxorder 
    break
end
    
fprintf('Restart #%d \n',iii)

if ~param.fullupdate
    Wr = zeros(n(1),0);
else
    Wr = W;
end

    lambda = seprand([1,n(2:end)]);
    lambda.F{1,1}=1;
    alpha0=0;
    
for i=1:param.maxorder-size(W,2)
fprintf('arnoldi iterate #%d - ',i);

result.stolambda{i+(~param.fullupdate)*size(W,2)}=lambda;

lAl = lambda'*A*lambda;
lbu = lambda'*bu;
M = reduce(lAl,1); 
S = reduce(lbu,1);
U = M\S;
alpha = norm(U);
U = U/alpha;

for jjj=1:size(Wr,2)
   U = U - Wr(:,jjj)*(Wr(:,jjj)'*U); 
end

H = norm(U);
fprintf('residu = %.3e',H)
if norm(U)<param.orthocrit
    disp(' -> break')
break
end
 fprintf('\n')
U = U/H;
Wr = [Wr,full(U)];

AU = mtimes(A,U,1);
UAU = mtimes(U',AU,1);
Ubu = mtimes(U',bu,1);
if dim>2
[lambda,ressolver] = solve(UAU,Ubu,solverlambda);
result.resultlambda = [result.resultlambda , {ressolver} ];
else
M = reduce(UAU,2);
S = reduce(Ubu,2);
lambda.F{2} = M\S;
end


end  %%%%  END ITERATION ARNOLDI

if param.fullupdate
W=Wr;
bu=b;
else
W = [W,Wr];
end

AWr = mtimes(A,Wr,1);
WrAWr = mtimes(Wr',AWr,1);
Wrbu = mtimes(Wr',bu,1);
if dim>2
Lambda = solve(WrAWr,Wrbu,solverlambdaupdate);
else
M = assembleblock(WrAWr,2);
S = assembleblock(Wrbu,2);
Lambda = SEPMATRIX({eye(size(Wr,2)),reshape(M\S,n(2),size(Wr,2))});
Lambda = splitvectors(Lambda);
Lambda.alpha=ones(1,size(Wr,2));
end

if param.fullupdate
u = mtimes(Wr,Lambda,1);
bu = b - A*u;
else
ur = mtimes(Wr,Lambda,1);
u = u + ur;
bu = bu - A*ur;
end

mu = size(W,2);
switch param.errorindicator
    case 'reference'
errorder(mu)=normND(expand(u)-param.reference)/normref;                
    case 'residual'
errorder(mu)=norm(bu)/normref;        
    otherwise
errorder(mu)= norm(u-u0)/norm(u);
end
u0=u;

fprintf('order #%d - error = %d \n',size(W,2),errorder(mu))

if errorder(mu)<param.tol
    hasconv = 1;
 break
end

end

 
result.error = errorder;


function M = reduce(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
a = A.alpha.*prod(cell2mat(A.F),2)';
M = a(1)*Z{1};
for i=2:A.m
M = M + a(i)*Z{i};
end

return

% function M = assembleblock(A,j)
% 
% Z = A.F(:,j);
% A.F(:,j) = [];
% n = size(Z{1},1);
% nb = size(A.F{1},1);
% B = cell(A.m,1);
% for i=1:A.m
% B{i} = A.F{i,1};
% for k=2:size(A.F,2)   
% B{i}=B{i}.*A.F{i,k};    
% end
% end
% 
% if size(Z{1},2)>1
% M = sparse([],[],[],n,n,nnz(Z{1})*nb^2);
% for ii=1:nb
% I = (ii-1)*n+(1:n);
% for jj=1:nb
% J = (jj-1)*n+(1:n);
% MIJ = (A.alpha(1)*B{1}(ii,jj))*Z{1};
% for i=2:A.m
% MIJ = MIJ + (A.alpha(i)*B{i}(ii,jj))*Z{i};
% end
% M(I,J)=MIJ;
% end
% end
% else 
% M = zeros(nb*n,1);
% for ii=1:nb
% I = (ii-1)*n+(1:n);
% MI = (A.alpha(1)*B{1}(ii))*Z{1};
% for i=2:A.m
% MI = MI + (A.alpha(i)*B{i}(ii))*Z{i};
% end
% M(I)=MI;
% end
% end
%   
% 
% return


function M = timesblock(A)

i=1;
Mi = A.F{i,1};
for k=2:A.dim
Mi = Mi.*A.F{i,k};    
end
M = Mi*A.alpha(i);

for i=2:A.m
Mi = A.F{i,1}; 
for k=2:A.dim
Mi = Mi.*A.F{i,k};    
end
M = M + Mi*A.alpha(i);
end

return

