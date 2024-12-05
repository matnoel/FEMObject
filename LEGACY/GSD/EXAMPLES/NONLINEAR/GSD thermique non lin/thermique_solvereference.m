function [qpc,result] = thermique_solvereference(tol,S,PC,X,cas)
if nargin==4
    cas=1
end
thermique_nonlin_forms

N = NEWTONSOLVER('type','full','tol',tol,...
    'tolreact',1e-2,'increment',true);

if nbsec==1
    f = X{3}*l{S}(:);
else
    f = X{3}*l{S}(:)+X{4}*l2{S}(:);
end


if cas==1
    fsolve = @(A,b) cgs(A,b,getparam(N,'tol')/100);
else
    fsolve = @(A,b) solve(A,b);
end

clock0=clock;
[qpc,result] = solve(N,f,@(u) calc_Bu(S,u,PC,a,n,X),...
    @(u) calc_tang(S,u,PC,a,n,dn,X,cas),[],[],fsolve);
etime(clock,clock0)



function Bt=calc_tang(S,u,PC,a,n,dn,X,cas)

%Bt = X{1}*a{S}(:,:);
%Bt = X{1}*a{S}(:,:)+X{2}*();
if cas==1
    Bt = X{1}*a{S}(:,:)+X{2}*dn{S}(:,:,expect(u),expect(u));
else
    Bt = expect(X{1})*a{S}(:,:)+expect(X{2})*dn{S}(:,:,expect(u),expect(u));
end
return


function b=calc_Bu(S,u,PC,a,n,X)
s=size(u);
u=expand(u);
b = X{1}*(a{S}(:,:)*u);
%for k=1:length(PC)
b = b+X{2}*n{S}(:,u);
%end

return
