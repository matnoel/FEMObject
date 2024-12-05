function A=calc_trimatrix(SM,a,varargin)
% function A=calc_trimatrix(SM,MULTILINFORM,varargin)
%   SM : SEPMODEL
%   A  : SEPMATRIX
% Attention : A est une HSEPMATRIX du format le plus compact.
% Ne pas oublier de changer de format ((H)SEPMATRIX) si necessaire.


p=getp(a);

if isa(a,'TRILINFORM')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%            (:,:,:)            %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Trois cas sont a conciderer :
    % 1 - Sortie de type scalaire :
    %     \int{ u v w }
    %     \int{ grad(u).grad(v) w }  || \int{ u grad(v).grad(w) }
    %     Les cas 2 avec une dimension spatiale=1 sont aussi a prendre en
    %     compte.
    %     La sortie est alors une (H)SEP
    % 2 - Sortie de type vecteur :
    %     \int{ grad(u) v w } || \int{ u grad(v) w } || \int{ u v grad(w) }
    %     La sortie est alors un veteur-cell {(H)SEP}
    % 3 - Sortie de type matrice :
    %     \int{ grad(u).v.grad(w) }
    %     La sortie est alors une matrice-cell {(H)SEP} (v est une matrice)
    
    cas1=all(p==[0 0 0]) || all(p==[1 1 0]) || all(p==[0 1 1]);
    cas2=all(p==[1 0 0]) || all(p==[0 1 0]) || all(p==[0 0 1]);
    cas3=all(p==[1 0 1]);
    if cas1
        A=calc_trimatrix_SCALAR(SM,a,varargin{:});
    elseif cas2
        A=calc_trimatrix_VECTOR(SM,a,varargin{:});
    elseif cas3
        A=calc_trimatrix_MATRIX(SM,a,varargin{:});
    else
        error('Cas non repertorie')
    end
    
else
    error('Ne calcule (pour l''intant) que des TRILINFORM(:,:,:)');
end


end







