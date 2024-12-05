function [P,x,y] = pdfsample(us,x,y)
% function [P,xP] = pdfsample(us,x)
% us : vecteur de realisations
% x  : vecteur x=[x0,x1,...,xN]
% sorties : 
% P : vecteur donnant la densite de proba associee � us aux points xP
% xP : points milieux de l'ensemble de points x (pour l'affichage de P)
%
% function [P,xP,yP] = pdfsample(us,x,y)
% meme chose en 2D
% us doit contenir 2 colonnes, us(:,i) contient des r�alisations d'une
% variable aleatoire 
% x et y sont un ensemble de points correspondant � des r�alisations des
% deux variables
% P est une matrice dont les composantes sont les probabilites (renormalisee)
% d'etre dans less patch [xi,xi+1]*[y_j,y_j+1]
% xP, yP : points milieux des patchs (-> pour l'affichage par surf ou mesh)

if any(size(us)==1)
    us=us(:);
end

if size(us,2)==1
    P=zeros(1,length(x)-1);

    for k=1:length(x)-1
        P(k)=length(find(us>= x(k) & us<x(k+1)))/(x(k+1)-x(k))/length(us) ;  
    end

    x=[x(1:end-1)+x(2:end)]/2;
else
    P=zeros(length(y)-1,length(x)-1);
    x=x(:);
    y=y(:);
    dx = x(2:end)-x(1:end-1);
    dy = y(2:end)-y(1:end-1);
    DX = dy*dx';

    U1 = repmat(us(:,1),1,length(x)-1);
    U2 = repmat(us(:,2),1,length(y)-1);
    Y = repmat(y',size(us,1),1);
    X = repmat(x',size(us,1),1);
    rep1 = U1>=X(:,1:end-1) & U1<X(:,2:end);
    rep2 = U2>=Y(:,1:end-1) & U2<Y(:,2:end);


    for k=1:length(x)-1
        P(:,k) = sum(repmat(rep1(:,k),1,length(y)-1) & rep2,1)'./DX(:,k)/size(us,1); 
    end



%for k=1:length(x)-1
%for l=1:length(y)-1
%P(l,k)=  length(find(rep1(:,k) & rep2(:,l)))/DX(k,l)/size(us,1); 
%end
%end

%tic
    for k=1:length(x)-1
        for l=1:length(y)-1
            P(l,k)=length(find(us(:,1)>= x(k) & us(:,1)<x(k+1) & us(:,2)>= y(l) & us(:,2)<y(l+1)))/(x(k+1)-x(k))/(y(l+1)-y(l))/size(us,1) ;  
        end
    end


    x=[x(1:end-1)+x(2:end)]/2;
    y=[y(1:end-1)+y(2:end)]/2;
    [x,y]=meshgrid(x,y);

end
