function varargout=KL(RF,m,S,varargin)
% function [V,L]=KL(RF,m,S)
% RF : RANDGAUSSIAN
% S : MODEL de type MAILLAGE
% m : ordre de la decomposition
%
% function [V,L]=KL(RF,m,S,'rescale')
% les vecteurs de la decomposition de Karhunen-Loeve tronqu�e sont
% modifies de tel maniere a obtenir un champ de distribution marginale
% gaussienne centree reduite
%

disp('---- Decomposition de Karhunen-Loeve , champ gaussien ----')
if ~isa(S,'MODEL')
    error('troisieme argument doit etre un MODEL')
end

    
S=createddlnode(S,DDL('U'));
D=calc_massgeom(S);
x = getcoord(S.node);
disp('Compute correlation matrix ...')

R=correl(RF,x,x);
s=std(RF,x);
R = R.*(s*s');

opt.disp=0;
disp('Solve eigenvalue problem ...')

C = D*R*D ; 
C = (C+C')/2;

[V,L]=eigs(C,D,m+1,'LM',opt);
smax = sign(V(1,:));


if any(smax)==0
error('mauvais point, le mode est nul')
end
V= V*diag(smax);

V = V/diag(sqrt(sum(V.*(D*V),1)));


lam=diag(L);
[lam,ia]=sort(lam,'descend');
err = sqrt(abs(lam(end))/sum(lam));

lam = lam(1:m);
V=V(:,ia(1:m));

if ischarin('rescale',varargin)
V = V.*repmat(s./sqrt(sum(V.^2*diag(lam),2)),1,m);   
end



disp(['order ' num2str(m)])
fprintf('error estimate = %3d\n',err);

L=diag([1;sqrt(abs(lam))]);
v{1}=mean(RF,x);    
for k=1:m
v{k+1}=(V(:,k));
l{k+1}=RVNORMAL();
end


if ischarin('modes',varargin) || nargout==2

V = [v{1},V];
if nargout==2
varargout{1} = V;
varargout{2} = L;
else
varargout{1} = V*L;    
end

else
    
   
p = getcharin('order',varargin,1);

PCg = PCMODEL(RANDVARS(l(2:end)),'order',p);
l = vertcat(one(PCg),PCg{:});
varargout{1} = FENODEFIELD(PCRADIALMATRIX(v,size(v{1}),L*l));
end

