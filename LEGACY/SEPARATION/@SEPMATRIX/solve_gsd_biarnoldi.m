function [u,result,W1,W2] = solve_gsd_biarnoldi(A,b,solver)

dim = getdim(A);
param = getparam(solver);

n = size(b);
bu = b;

u = SEPMATRIX(dim);
u0=u;
W1 = zeros(n(1),0);
W2 = zeros(n(2),0);

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

if size(W1,2)==param.maxorder 
    break
end
    
fprintf('Restart #%d \n',iii)

if ~param.fullupdate
    W1r = zeros(n(1),0);
    W2r = zeros(n(2),0);
else
    W1r = W1;
    W2r = W2;
end

    w2 = rand(n(2),1);
    alpha0=0;
    
for i=1:param.maxorder-size(W1,2)
fprintf('arnoldi iterate #%d - ',i);

Aw2 = mtimes(A,w2,2);
w2Aw2 = mtimes(w2',Aw2,2);
w2bu = mtimes(w2',bu,2);
M = reduce(w2Aw2,1);
S = reduce(w2bu,1);
w1 = M\S;
alpha1 = norm(w1);
w1=w1/alpha1;

for jjj=1:size(W1r,2)
   w1 = w1 - W1r(:,jjj)*(W1r(:,jjj)'*w1); 
end

H1 = norm(w1);
fprintf('residu 1= %.3e',H1)
if H1<param.orthocrit
    disp(' -> break')
break
end
w1 = w1/H1;

Aw1 = mtimes(A,w1,1);
w1Aw1 = mtimes(w1',Aw1,1);
w1bu = mtimes(w1',bu,1);
M = reduce(w1Aw1,2);
S = reduce(w1bu,2);
w2 = M\S;
alpha2 = norm(w2);
w2=w2/alpha2;

for jjj=1:size(W2r,2)
   w2 = w2 - W2r(:,jjj)*(W2r(:,jjj)'*w2); 
end

H2 = norm(w2);
fprintf(' , residu 2= %.3e',H2)
if H2<param.orthocrit
    disp(' -> break')
break
end
w2 = w2/H2;

fprintf('\n')
W1r = [W1r,full(w1)];
W2r = [W2r,full(w2)];


end  %%%%  END ITERATION ARNOLDI

if param.fullupdate
W1=W1r;
W2=W2r;
else
W1 = [W1,W1r];
W2 = [W2,W2r];
end
u=SEPMATRIX(dim);
for kkk=1:size(W1,2)
    u = u + SEPMATRIX({W1(:,kkk),W2(:,kkk)},1);
end
u = solve_updatetucker(A,b,u,u);
bu = b - A*u;


mu = size(W1,2);
switch param.errorindicator
    case 'reference'
errorder(mu)=normND(expand(u)-param.reference)/normref;                
    case 'residual'
errorder(mu)=norm(bu)/normref;        
    otherwise
errorder(mu)= norm(u-u0)/norm(u);
end
u0=u;

fprintf('order #%d - error = %d \n',size(W1,2),errorder(mu))

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

