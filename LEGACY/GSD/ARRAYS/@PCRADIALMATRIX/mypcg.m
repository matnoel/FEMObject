function qpc=mypcg(K,f,varargin)
%function q=mypcg(K,f)
% 
% Gradient conjugue pour un systeme stochastique   K*q=f
% K PCRADIAL
% f PCRADIAL ou PCCELL ou PCARRAY
% q : PCARRAY
% 
% function q=pcg(K,f,option,valeur)
% option = 'tol'   valeur = precision souhaitee (1e-10 par defaut) 
% option = 'noreortho' pas de reorthogonalisation
% option = 'noprecond' pas de preconditionnement
% option = 'doublereortho' double reorthogonalisation
% option = 'nstomax' valeur = nombre maxi de vecteurs stockes
% option = 'display' affichage des informations 

PC=getPC(K);      
P=getP(PC);
n=size(K,1);


paramgc.prec = getcharin('tol',varargin,1e-10);
paramgc.reortho = ~ischarin('noreortho',varargin);
paramgc.precond = ~ischarin('noprecond',varargin);
paramgc.deuxortho = getcharin('doublereortho',varargin);
paramgc.nstomax=getcharin('nstomax',varargin,n*(P+1)) ;
dislay_ = ischarin('display',varargin);

if dislay_==1
Precondtype{1}='Auncun';Precondtype{2}='Bloc Diagonal (esperance de K)';
Reorthotype{1}='non';Reorthotype{2}='oui';
fprintf('\n ------- Gradient Conjugue -------\n')
fprintf(' Precision   arret    : %3d\n', paramgc.prec )
disp([' Preconditionnement  : ' Precondtype{paramgc.precond+1}])
disp([' Reorthogonalisation : ' Reorthotype{paramgc.reortho+1}])
if paramgc.reortho==1
fprintf(' Nb max de fonctions : %3d\n' , paramgc.nstomax)
end
end

switch paramgc.precond
case 0
ML=speye(n);
case 1
ML = expect(K);
if ischarin('cholinc',varargin)
 ML = cholinc(ML,'0');
else
 ML = chol(ML,'lower')';
end
end
global masse

masse=getmasse(K);
K = K.V;

f=double(expand(f));

q  = zeros(n,P+1) ;
r  = f  ; 

z  = ML\(ML'\(r)) ;
normf= sqrt(prod_q1q2(f,f)) ;


for i=1:n*(P+1)
if dislay_==1
    fprintf('iteration %3d : ',i)
end
r0=r;z0=z;
if paramgc.reortho==1
isto=mod(i-1,paramgc.nstomax)+1; 
else
isto=1;
end


W{isto} = z0 ;

KW{isto} = prod_Kq(z0,K) ;

WKW{isto} = prod_q1q2(z0,KW{isto}) ;

a=-prod_q1q2(r0,z0)/WKW{isto};
q=q - a*W{isto} ;
r = r0 + a * KW{isto} ;
z = ML\(ML'\(r)) ; normz = sqrt(prod_q1q2(z,z));

scanortho=[1:isto];

if (paramgc.reortho==1 && i>paramgc.nstomax) 
scanortho=[isto+1:paramgc.nstomax,scanortho];
end

for j=scanortho
z = z - (prod_q1q2(z,KW{j})/WKW{j})* W{j} ;
end
if paramgc.deuxortho==1
for j=scanortho
z = z - (prod_q1q2(z,KW{j})/WKW{j})* W{j} ;
end
end

err=sqrt(prod_q1q2(r,r))/normf;

if dislay_==1    
    fprintf('erreur = %3d \r',full(err));
end
if err<paramgc.prec
if dislay_==1    
fprintf('\n ---->  converge en %3d iterations',i);
end
break
end

end

if dislay_ ==1
    fprintf('\n')
end

qpc=PCMATRIX(q,[n,1],PC);

return

function a=prod_q1q2(q1,q2)

a=q1(:)'*q2(:);

return

function x=prod_Kq(q,Kglob)

global masse

x=zeros(size(q));
for k=1:length(Kglob)
x=x+Kglob{k}*q*masse{k}';
end

return


