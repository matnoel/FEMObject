function [ut,result,vt] = dsolve(L,b,A,B,u0,solver)
% function [ut,result,vt] = dsolve(L,b,A,B,u0,solver)
% Discontinuous Galerkin solver: solves Au'+Bu=b
% L: DGTIMESOLVER
% b: right-hand side
% A,B: double or random matrices
% u0: vector of initial conditions
% solver: function for solving linear systems
%         solver(K,f) solves Ku=f
% ut: TIMEMATRIX of solution u
% vt: TIMEMATRIX of velocity v=u'
% result: struct containing outputs
% result.totaltime : total CPU time

tStart = tic;

display_ = getparam(L,'display');
if display_
    fprintf('\n ---------------------------------------');
    fprintf('\n ------------ DGtime solver ------------')
    fprintf('\n ---------------------------------------\n');
end

T = gettimemodel(L);
% t = gett(T);
nt = getnt(T);
dt = getdt(T);
p  = getp(L);

if p>2
    error('DG pas programme pour un ordre superieur a 2')
end

[b,A,B] = init_resolution(T,b,A,B);

if isuniform(T) && ~istime(A) && ~istime(B)
    switch p
        case 0
            C = A+B*dt(1);
        case 1
            C = [1/2*A+B*(dt(1)/3) , 1/2*A+B*(dt(1)/6) ; ...
                -1/2*A+B*(dt(1)/6) , 1/2*A+B*(dt(1)/3)];
        case 2
            C = [1/2*A+B*(dt(1)*2/15) , 2/3*A+B*(dt(1)/15), -1/6*A+B*(-dt(1)/30) ; ...
                -2/3*A+B*(dt(1)/15) , B*(dt(1)*8/15), 2/3*A+B*(dt(1)*1/15) ; ...
                1/6*A+B*(-dt(1)/30) , -2/3*A+B*(dt(1)/15), 1/2*A+B*(dt(1)*2/15) ];
    end
else
    C = [];
end

if nargin>=6 && ~isempty(solver)
    solver =  fcnchk(solver);
else
    if ~isempty(C) && isa(C,'double') && getparam(L,'lu')
        if size(C,1)>100 && issparse(C)
            perm = symrcm(C);
            [LC,UC] = lu(C(perm,perm));
            solver = @(C,b) unperm(solve(UC,(LC\b(perm))),perm);
        else
            [LC,UC] = lu(C);
            solver = @(C,b) solve(UC,(LC\b));
        end
    else
        solver = @solve;
    end
end

n = size(A,1);
if nargin<5 || isempty(u0)
    u0 = zeros(n,1);
end

if display_
    fprintf('DG order %d : ',p);
end
ut = cell(1,nt);
switch p
    case 0
        for i=1:nt
            if display_
                pourcentage(i,nt)
            end
            
            if i>1
                bi = b{i}*dt(i)+A*ut{i-1};
            else
                bi = b{i}*dt(i)+A*u0;
            end
            if ~isuniform(T) || istime(A) || istime(B)
                if istime(A)
                    Ai = getmatrixatstep(A,i);
                else
                    Ai = A;
                end
                if istime(B)
                    Bi = getmatrixatstep(B,i);
                else
                    Bi = B;
                end
                C = Ai+Bi*dt(i);
            end
            
            if nargin(solver)==3 && i>1
                [ut{i},result.solver{i}] = solver(C,bi,ut{i-1});
            elseif nargin(solver)==3
                [ut{i},result.solver{i}] = solver(C,bi,[]);
            else
                [ut{i},result.solver{i}] = solver(C,bi);
            end
        end
        
    case 1
        for i=1:nt
            if display_
                pourcentage(i,nt)
            end
            
            if i>1
                U0 = ut{i-1}(p*n+1:end,:) ;
            else
                U0 = u0;
            end
            
            bi = [b{i*(p+1)-1}*(dt(i)/3) + b{i*(p+1)}*(dt(i)/6) ; ...
                b{i*(p+1)-1}*(dt(i)/6) + b{i*(p+1)}*(dt(i)/3) ] + [A*U0;zeros(n,1)] ;
            
            if ~isuniform(T) || istime(A) || istime(B)
                if istime(A)
                    Ai = (getmatrixatstep(A,i*(p+1)-1)+getmatrixatstep(A,i*(p+1)))/2;
                else
                    Ai = A;
                end
                if istime(B)
                    Bi = (getmatrixatstep(B,i*(p+1)-1)+getmatrixatstep(B,i*(p+1)))/2;
                else
                    Bi = B;
                end
                C = [1/2*Ai+Bi*(dt(i)/3) , 1/2*Ai+Bi*(dt(i)/6) ; ...
                    -1/2*Ai+Bi*(dt(i)/6) , 1/2*Ai+Bi*(dt(i)/3)];
            end
            
            if nargin(solver)==3 && i>1
                [ut{i},result.solver{i}] = solver(C,bi,ut{i-1});
            elseif nargin(solver)==3
                [ut{i},result.solver{i}] = solver(C,bi,[]);
            else
                [ut{i},result.solver{i}] = solver(C,bi);
            end
        end
        
    case 2
        for i=1:nt
            if display_
                pourcentage(i,nt)
            end
            
            if i>1
                U0 = ut{i-1}(p*n+1:end,:);
            else
                U0 = u0;
            end
            
            bi = [b{i*(p+1)-2}*(dt(i)*2/15) + b{i*(p+1)-1}*(dt(i)*1/15) + b{i*(p+1)}*(-dt(i)/30) ; ...
                b{i*(p+1)-2}*(dt(i)*1/15) + b{i*(p+1)-1}*(dt(i)*8/15) + b{i*(p+1)}*(dt(i)/15) ; ...
                b{i*(p+1)-2}*(-dt(i)/30) + b{i*(p+1)-1}*(dt(i)*1/15) + b{i*(p+1)}*(dt(i)*2/15) ] ...
                + [A*U0;zeros(2*size(A,1),1)];
            
            if ~isuniform(T) || istime(A) || istime(B)
                if istime(A)
                    Ai = getmatrixatstep(A,i*(p+1)-1);
                else
                    Ai = A;
                end
                if istime(B)
                    Bi = getmatrixatstep(B,i*(p+1)-1);
                else
                    Bi = B;
                end
                C = [1/2*Ai+Bi*(dt(i)*2/15) , 2/3*Ai+Bi*(dt(i)/15), -1/6*Ai+Bi*(-dt(i)/30) ; ...
                    -2/3*Ai+Bi*(dt(i)/15) , Bi*(dt(i)*8/15), 2/3*Ai+Bi*(dt(i)*1/15) ; ...
                    1/6*Ai+Bi*(-dt(i)/30) , -2/3*Ai+Bi*(dt(i)/15), 1/2*Ai+Bi*(dt(i)*2/15) ];
            end
            
            if nargin(solver)==3 && i>1
                [ut{i},result.solver{i}] = solver(C,bi,ut{i-1});
            elseif nargin(solver)==3
                [ut{i},result.solver{i}] = solver(C,bi,[]);
            else
                [ut{i},result.solver{i}] = solver(C,bi);
            end
        end
end

result.totaltime = toc(tStart);
if display_
    fprintf('Elapsed time = %f\n',result.totaltime)
end

Ut = cell(1,length(T));
for i=1:nt
    ut{i} = reshape(ut{i},n,p+1);
    for k=0:p
        Ut{(p+1)*(i-1)+k+1} = ut{i}(:,k+1);
    end
end

ut = TIMEMATRIX(Ut,T,[n,1]);

if nargout>=3
    vt = derivative(L,ut,b,A,B);
end


function [u,flag] = unperm(u,perm)
u(perm,:)=u;
flag=0;
return

