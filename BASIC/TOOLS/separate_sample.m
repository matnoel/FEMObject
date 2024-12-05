function [X,varargout] = separate_sample(Xs,type,varargin)
% function [Xss,varargout] = separate_sample(Xs,varargin)
% Xs : echantillons d'une variable aeatoire
% separation des echantillons � probabilit� fix�e ou � intervalle fix�
% Xss : blocks d'�chantillons (sous forme de cellule)
%
% function [Xss,x,P] = separate_sample(Xs,'proba',P)
% P : vecteur de taille n 
% separation en n block d'echantillons de probabilit� Pi = P(i), i=1...n
% (si P scalaire Pi=P (equiprobable) P doit pouvoir s'�crire P=1/n, n entier)
% x : points de s�paration des blocks 
% P : reactualisation des proba des intervalles (car decomposition non exacte)
%
% function [Xss,x,P] = separate_sample(Xs,'interval',x)
% x : vecteur de taille N+1 (points de s�paration des blocks) (tri� par ordre croissant)
% si des r�alisations existent en dehors de l'intervalle [x(1),x(end)],
% on redefinit x(1) et x(end)
% (si x est un scalaire : on separe [min(Xs),max(Xs)] en x intervalles de meme longueur )

N = length(Xs);
Xs = sort(Xs);

switch type
case 'proba'
    P = varargin{1};
    if length(P)==1
        n = floor(1/P);
        if n~=1/P
            error('la probabilit� doit s''ecrire 1/n avec n entier')
        end
        P = repmat(P,1,n);
    else
        n=length(P);
    end
    if any(P<=0 | P>=1) || sum(P)~=1
        error('les probabilit� doivent etre comprises entre 0 et 1 strictement et la somme doit etre egame a 1')
    end

    x=min(Xs);
    for i=1:n
        rep = [1:ceil(N*P(i))];
        X{i} = Xs(rep);
        x = [x,Xs(rep(end))];
        if i>1
            x(end-1) = 1/2*(x(end-1)+Xs(rep(1)));    
        end
        Xs(rep)=[];
    end

    for i=1:n
        P(i) = numel(X{i})/N;
    end
case 'interval'
    x= varargin{1};
    if length(x)==1        
        x = linspace(min(Xs),max(Xs),x+1);
    else
        x = sort(x);
        if x(1)>min(Xs)
            x(1)=min(Xs);
            warning('on redefinit la borne inf des intervalles pour coller aux echantillons')
        end
        if x(end)<max(Xs)
            x(end)=max(Xs);
            warning('on redefinit la borne sup des intervalles pour coller aux echantillons')
        end
    end

    rep = find(Xs>=x(1) && Xs<=x(2));
    X{1} = Xs(rep);
    P(1) = numel(rep)/N;
    for i=2:length(x)-1
        rep = find(Xs>x(i) && Xs<=x(i+1)); 
        X{i} = Xs(rep);
        P(i) = numel(rep)/N;
    end


end


if nargout>1
    varargout{1}=x;
end

if nargout>2
    varargout{2}=P;
end
