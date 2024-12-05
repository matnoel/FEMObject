function masse=calc_masse(h,p,p2,varargin)
% masse=calc_masse(h,p,p2)
% calcul des matrice de masse 
% masse : MULTIMATRIX (p2+1)-by-(p2+1) et de longueur (p+1) 
% E(hk * a * b) = a'*masse{k+1}*b  ou a et b sont les vecteurs
% representant les coeff de a et b sur la base (h0, ... , hp2) 
if nargin==2 || isempty(p2)
   p2=p ;
end
try
%load(['M' h.type])
end

if exist('masse') && length(masse)>=p+1 && size(masse,1)>=p2+1 && polycmp(h,hsave)
   masse = masse{1:p+1}(1:p2+1,1:p2+1);
   masse = MULTIMATRIX(masse,[p2+1,p2+1],[p+1,1]);
else
    if ischarin('display',varargin)
 fprintf(['Calcul des masse pour ' h.type ' ... '])
    end
    
intxn=calc_intxn(h,[0:p+2*p2]);
masse=cell(p+1,1);
for k=0:p
masse{k+1}=sparse(p2+1,p2+1);
for i=0:p2
for j=i:p2
masse{k+1}(i+1,j+1)=moment(h,[i,j,k],intxn);
masse{k+1}(j+1,i+1)=masse{k+1}(i+1,j+1);
end
end
if ischarin('display',varargin)
pourcentage(k+1,p+1);
end
end

masse=MULTIMATRIX(masse,[p2+1,p2+1],[p+1,1]);
hsave=h;
save(['M' h.type],'masse','hsave')
if ischarin('display',varargin)
    fprintf('\n')
end
end

if ischarin('double',varargin) || ischarin('cell2mat',varargin)
   masse = cell2mat(masse); 
end

if ischarin('double',varargin)
    warning('remplacer argument double par cell2mat')
end