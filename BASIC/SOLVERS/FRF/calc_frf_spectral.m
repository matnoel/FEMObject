function [ur,ui,u]=calc_frf_spectral(wband,p,K,M,f,alphaw,betaw,solveur)
% function [ur,ui,u]=calc_frf_spectral(wband,p,K,M,f,alphaw,betaw,solveur)

wmax = max(max(wband));
wmin = min(min(wband));

xi = (wband-wmin)/(wmax-wmin);
if size(xi,1)==1
    h = POLYFE(xi);
else
    h = POLYFE(xi,[],'elem');    
end

PC = POLYCHAOS(h,p);
PC = calc_masse(PC);
w = decompfun(PC,[],[],@(xi) wmin+xi*wmax);

% RV = RANDVARS(RVUNIFORM(wmin,wmax));
% PC = PCMODEL(RV,'order',p,'fedim',1,'femesh',{4});
% w = PC{1};

if nargin<=5 || (isempty(alphaw) && isempty(betaw))
    amort=0;
else
    amort=1;
end


if amort==1
    Z = sparse(size(K,1),size(K,2));
    Aw = [K,Z;Z,K] + w*alphaw*[Z,-M;M,Z] + w*betaw*[Z,-K;K,Z] -(w*w)*[M,Z;Z,M];
    bw = one(PC)*[f;0*f];
else
    Aw = K - w.^2*M;
    bw = one(PC)*f;
end


if nargin<8
    solveur = @(A,b) solve(A,b);
end
u = solveur(Aw,bw);

ur = u(1:size(K,1));
ui = u(1+size(K,1):end);



