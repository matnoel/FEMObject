function masse=calc_masse(h,p,p2,varargin)
% masse=calc_masse(h,p,p2)
%
% calcul des matrice de masse 
% masse : MULTIMATRIX (p2+1)-by-(p2+1) et de longueur (p+1) 
% E(hk * a * b) = a'*masse{k+1}*b  ou a et b sont les vecteurs
% representant les coeff de a et b sur la base (h0, ... , hp2) 

if nargin==2
p2=p;    
end
if ischarin('display',varargin)
fprintf(['Calcul des masse pour ' get(h,'type') ' ... '])
end
param=get(h,'param');
n=param.n;
I=param.I;
dx = param.dx ;
subpoints = param.subpoints;
reppoints = param.reppoints;

gauss=calc_gausslobattopoints(POLYLEGENDRE(),p2+1);
M = calc_masse(POLYLAGRANGE(gauss.coord),p,p2);


masse=cell((p+1)*n,1);
for i=1:p*n+1
masse{i}=sparse((p2)*n+1,(p2)*n+1);
end

for i=1:n
repI = reppoints{i} ;
for k=0:p
masse{reppoints{i}(k+1)}(repI,repI) =masse{reppoints{i}(k+1)}(repI,repI)+ M{k+1}*(dx(i))/(I(end)-I(1));
end
end

masse = MULTIMATRIX(masse);

if ischarin('double',varargin) || ischarin('cell2mat',varargin)
   masse = cell2mat(masse); 
end

if ischarin('double',varargin)
    warning('remplacer argument double par cell2mat')
end
if ischarin('display',varargin)
fprintf('\n')
end
