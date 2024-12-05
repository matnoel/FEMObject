function A = calc_matrix(SM,a,varargin)
% function A = calc_matrix(SM,MULTILINFORM,varargin)
%   SM : SEPMODEL
%   A  : SEPMATRIX


p = getp(a);


if isa(a,'TRILINFORM')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% (V,:,:) || (:,V,:) || (:,:,V) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    V = varargin(1:3);
    varargin = varargin(4:end);
    locV = find(cellfun(@isempty,V(:))==0);
    V = V{locV};
    
    A = calc_trimatrix(SM,a,[],[],[],varargin{:});
    
    % V impose son style :
    A = impose_style(A,V);
    
    
elseif isa(a,'BILINFORM')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%           (:,:)           %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Astuce de la mort qui tue :
    %   1-Demander la T=tlf( [p 0] )
    %   2-Evaluter T{SM}(:,:,1),
    %     ou bien  T{SM}(:,:,k) si k specifie dans BILINFORM
    
    locV=3;
    if ~isempty(getk(a)) && ~isempty(getpk(a))
        % Un champs k est specifie dans BILINFORM :
        %   b = BILINFORM(~,~,k,pk)
        V = getk(a);
        
        % Detecter si k est un CHamps par Element :
        if isa(V,'cell'), sV = size(V{1});
        else              sV = size(V);
        end
        if any(sV~=getnbddl(SM))
            varargin = [varargin 'CHE' 3];
        end
        
        % La forme trilineaire
        ta = TRILINFORM(p(1),p(2),getpk(a));
        
        % avec un cas problematique qu'il ne faut pas oublier:
        if isa(V,'SEP')
            % cas problematique 1 :
            if [p(1),p(2),getpk(a)]==[1 0 1];
                locV = 2;
                ta = TRILINFORM(p(1),getpk(a),p(2));
                if ischarin('CHE',varargin)
                    [~,nia] = ischarin('CHE',varargin);
                    varargin{nia+1} = 2;
                end
            end
        elseif isa(V,'cell') && min(size(V))~=1 && all(p==[1 1])
            % cas problematique 2 :
            % V est une matrice : forme du type TLF(1,0,1)
            locV = 2;
            ta = TRILINFORM(p(1),0,p(2));
            if ischarin('CHE',varargin)
                [~,nia] = ischarin('CHE',varargin);
                varargin{nia+1} = 2;
            end
        end
        
        % Assemblage de la TLF :
        A = calc_trimatrix(SM,ta,[],[],[],varargin{:});
        % V impose son style :
        if isa(V,'cell')
            A = impose_style(A,V{1});
        else
            A = impose_style(A,V);
        end
    else
        % Evaluter T{SM}(:,:,1)
        V = one(SM);
        ta = TRILINFORM(p(1),p(2),0);
        A = calc_trimatrix(SM,ta,[],[],[],varargin{:});
        % A impose son style :
        V = impose_style(V,A);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Realiser la multiplication  %%%%
%%%%     A.V sur la dim locV      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(V,'SEP') % Cas scalaire
    if isa(A,'SEP')
        % Finalisation
        A = finalize(A,V,locV,varargin{:});
    elseif isa(A,'cell')
        % Finalisation
        A = cellfun(@(A) finalize(A,V,locV,varargin{:}),...
            A,'UniformOutput',false);
        if numel(A)==1
            A = A{1};
        end
    end
elseif isa(V,'cell')
    if isa(A,'SEP')
        error('Reflechis a ce que tu viens d''ecrire')
    elseif isa(A,'cell') && all(size(A)==size(V))
        % Deux possibilite :
        %  - V est un vecteur qui vient completer une tlf du type :
        %    [1 0 0]  [0 1 0]  [0 0 1]
        %  - V est une matrice qui vient completer le seul cas :
        %    [1 0 1]
        % Dans les deux cas, la sotie est une HSEP, ie sortie scalaire !
        A = cellfun(@(A,V) finalize(A,V,locV,varargin{:}),...
            A,V,'UniformOutput',false);
        % Sommation :
        for i=2:numel(A)
            A{1} = A{1}+A{i};
        end
        A = A{1};
    end
end



end



function A = finalize(A,V,locV,varargin)
% function A = finalize(A,V,locV,varargin)

switch order(A)
    case 3
        A = mtimeslike(A,V,@(u,v) SEPtimesSEP_FH_MULTM(u,v,locV),'MATRIX');
    otherwise
        A = mtimeslike(A,V,@(u,v) SEPtimesSEP_FH_MULTM(u,v,locV),'');
end
end


function A=impose_style(A,V)
% function A=impose_style(A,V)

if isa(A,'cell')
    if isa(V,'SEP')
        A = cellfun(@(A)   change_tree(A,V),A,'UniformOutput',false);
    else
        A = cellfun(@(A,V) change_tree(A,V),A,V,'UniformOutput',false);
    end
else
    if isa(V,'SEP')
        A = change_tree(A,V);
    else
        A = change_tree(A,V{1});
    end
end

end


