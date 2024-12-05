%%
RV=RANDVARS(RVUNIFORM(1,2),2);
X = PCMODEL(RV,'pcg','order',4);
PC = getPC(X);
A = eye(3)*X{1}+[2,1,0;1,2,1;0,1,2]*X{2};
b = one(PC)*[1;0;0];
%A = eye(4)*X{1}+[2,1,0,0;1,2,1,0;0,1,2,1;0,0,1,2]*X{2};
%b = one(PC)*[1;0;0;0];

%
u = solve(A,b);
%%

A = calc_masse(A,PC);

%%

f = @(U) solve(U'*A*U,U'*b);
fm = @(U,um) solve(U'*A*U,U'*(b-A*um));

F = @(l) solve(assembleblock(expect(A,l,l)),assembleblock(expect(b,l)));
Fm = @(l,um) solve(expect(A,l,l),expect(b,l)-expectmtimes(A*l,um));

Ray = @(U) trace(full(expectmtimes(U'*b,f(U)')));
Raym = @(U,um) full(expectmtimes(U'*b,fm(U,um))-expectmtimes(U'*A*fm(U,um),um));

T = @(U) reshape(F(f(U)),size(U));
Tm = @(U,um) Fm(fm(U,um),um);

%%
um = PCRADIALMATRIX(size(b),PC);
%
m=2;
for i=1:m
    l=rand(PC);l=normalize(l);

    if getm(um)>0
        U0=Fm(l,um);    
        R0 = Raym(U,um);
    else
        U0=F(l);    
        R0 = Ray(U0);    
    end
    kmax=20;
    tolpfix = 1e-13;
    update = 0;
    for k=1:kmax
        if getm(um)>0
            l = fm(U0,um);
        else
            l = f(U0);    
        end
        l = normalize(l);
        if getm(um)>0
            U = Fm(l,um);
            alpha = norm(U);
            U = U/alpha;
            R = Raym(U,um);
        else
            U = F(l);
            alpha = norm(U);
            U = U/alpha;
            R = Ray(U);    
        end
%fprintf('m = %d , iteration %d , error Rayleigh = %.4d\n',i,k,abs((R-R0)/R));

        if sqrt(abs((R-R0)/R))<tolpfix
            fprintf('m = %d , iteration %d , error Rayleigh = %.4d\n',i,k,abs((R-R0)/R)); 
            break 
        elseif k==kmax
            fprintf('m = %d , iteration %d , error Rayleigh = %.4d , non convergence\n',i,k,abs((R-R0)/R));   
        end
        U0=U;
        R0=R;
    end
    um = um + U*(l*alpha);
    if update 
        W = double(getV(um));
        L = f(W);
        um = PCRADIALMATRIX(W,size(b),L);
    end
    fprintf('m = %d , error L2 = %.4d\n',i,norm(um-u)/norm(u));
    fprintf('m = %d , Raym = %.2d\n',i,R);
    fprintf('m = %d , Ray = %.2d\n',i,Ray(U));

end

%% test points fixes deflates

Wm = double(full(getV(um)));
for i=1:getm(um)
    if i==1
        e = norm(T(Wm(:,1))-Wm(:,1));
    else
        ui = truncate(um,1:i-1);
        e = norm(Tm(Wm(:,i),ui)-Wm(:,i));
    end
    fprintf('norm(Tm(Um)-Um) = %.3d\n',e)
end

%%
m=2
W0 = orth(rand(size(b,1),m));
R0 = Ray(W0);
for k=1:20
    W = T(W0);
    W = orth(full(W));
    R = Ray(W);
    fprintf('iteration %d , error Rayleigh = %.4d\n',k,abs((R-R0)/R)); 
    if sqrt(abs((R-R0)/R))<tolpfix
        break 
    end
    R0=R;
    W0=W;

end
e=norm(T(W)-W);
fprintf('norm(T(W)-W) = %.3d\n',e)

%%

a0=rand(3,1);
alpha = norm(W*a0);
a0 = a0/alpha;

[a,err,ierr] = nsol(a0,@(a) (T(W*a)-W*a),[1e-14,1e-14],[50,-1,1])
norm(T(W*a)-W*a)/norm(W*a)

%%
U0=rand(size(b));U0=U0/norm(U0);

for k=1:10
    U = T(U0);
    U=U/norm(U);
    U0=U;
end
U1=U;
l1=f(U1);
%%
Ur = getV(um,1);
options = optimset('Maxiter',100,'Display','iter','TolX',1e-12,'TolCon',1e-12,'MaxFunEvals',600);

U = fmincon(@(U) Ur'*U,rand(4,1),[],[],[],[],[],[],@(U) constTUU(U,T),options)

