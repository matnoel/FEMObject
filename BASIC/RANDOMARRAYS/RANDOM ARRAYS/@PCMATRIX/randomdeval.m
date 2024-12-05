function a = randomdeval(apc,j,x,RV)
% function a = randomeval(apc,x)
% evaluation de la PCMATRIX pour les realisations x
% M variables aleatoires
% x : double n-by-M (n : nombre de realisations)
% ou x : 1-by-M cell de n-by-1 double
%
% a : MULTIMATRIX de taille size(apc) et de longueur n si n>1 et numel(apc)>1
% a : double sinon
%
% function a = randomeval(apc,x,RV)
% RV indique les dimensions stochastiques de x (pour la correspondance avec apc)
% -> utilisation de [ok,rep] = ismember(apc,RV) et x = x(:,rep)
% RV peut etre un RANDPOLYS ou une RANDVARS ou encore un double

if nargin == 3
    RV = RANDVARS(apc.POLYCHAOS);
    if (isa(x,'cell') && length(x)~=length(RV)) || (isa(x,'double') && size(x,2)~=length(RV))
        error(['on attend ' num2str(length(RV)) ' jeux de valeurs' ])
    end
    if isa(x,'cell')
        x = [x{:}];
    end
end

if nargin==4
    x = transfer(RANDVARS(RV),RANDVARS(apc.POLYCHAOS),x);
end

H = dpolyval(apc.POLYCHAOS,j,x);

if iscell(apc)
    a = getvalue(multimtimes(H,apc.MULTIMATRIX));
else
    a = double(apc) * H';
end

% keyboard
if size(H,1)==1
    if iscell(apc)
        a = a{1};
    end
    a = reshape(a,size(apc));
elseif all(size(apc)==1)
    if iscell(apc)
        a = [a{:}];
    end
    a = reshape(a,1,size(H,1));
else
    a = MULTIMATRIX(a,size(apc),[size(H,1),1]);
end
