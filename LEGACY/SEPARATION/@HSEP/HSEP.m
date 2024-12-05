classdef HSEP < SEP
    properties
        tree=TREE();
    end
    
    methods
        function H = HSEP(F,alpha,varargin)
            % function A = HSEP()                   -> Initialisation
            % function A = HSEP(TREE)               -> Initialisation
            % function A = HSEP(struct/HSEP)        -> Copie
            % function A = HSEP( {(H)SEP} ,alpha)   -> Construction
            % function A = HSEP(HSEP,tree)          -> Morphisme (Yeah!)
            H=H@SEP();
            
            if nargin==0
                % Initialisation
                H.tree = TREE();
                
            elseif isa(F,'TREE')
                H.tree = F;
                H.dim = length(find(getconnect(F)==find(getconnect(F)==0)));
                
            elseif nargin==1 && (isa(F,'struct')||isa(F,'HSEP'))
                % Copie
                H.dim   = F.dim;
                H.m     = F.m;
                H.alpha = F.alpha;
                H.F     = F.F;
                H.tree  = F.tree;
                
            elseif isa(F,'cell') && isa(F{1},'SEP')
                % Construction
                H.F     = F;
                H.m     = size(F,1);
                H.dim   = size(F,2);
                if nargin==1, alpha = ones(1,H.m); end
                H.alpha = alpha;
                % Creation de l'arbre :
                T=cell(H.dim,1);
                for i=1:H.dim % dimension
                    if isa(F{1,i},'HSEP')
                        % On recupere l'arbre de la HSEP
                        T{i}=F{1,i}.tree;
                    else
                        % Creation de l'arbre de la SEP (ou autre...)
                        connect = [0, ones(1,F{1,i}.dim) ];
                        T{i}=TREE(connect);
                    end
                end
                H.tree = TREE(T);
                
            elseif nargin>=2 && isa(F,'SEP') && isa(alpha,'TREE')
                % Morphisme
                T=alpha;
                if isa(F,'HSEP')
                    % Le but ici est de concerver au maximum la reprentation compacte.
                    % On ne change l'arbre que lorsque c'est vraiment necessaire.
                    HSEPfils=str2func(class(F));
                    clF= class(F);
                    SEPfils=str2func(clF(2:end));
                    FT = F.tree;
                    if treecmp(FT,T)
                        % Premier niveau de decomposition identiques :
                        % On fait suivre la demande de transformation.
                        H=F;
                        H.tree=T;
                        Cvar = getCvar(T);
                        subnode =find(getconnect(T)==1);
                        for d=1:F.dim
                            extrdim = Cvar{subnode(d)};
                            ST      = subtree(T,subnode(d));
                            if length( find(getconnect(ST)==1)) < length(extrdim)
                                H.F(:,d)=cellfun(@(HF) HSEPfils(HF,ST),H.F(:,d),'UniformOutput',0);
                            end
                        end
                    else
                        % Recuperer la localisation premier etage des
                        % dimensions
                        dim2varFT=getdim2var(FT);
                        dim2varT =getdim2var(T);
                        for dFT=1:max(dim2varFT)
                            dimFT(dFT)=isallsame(dim2varT(dim2varFT==dFT));
                            varFT(dim2varFT==dFT)=dimFT(dFT);
                        end
                        for dT=1:max(dim2varT)
                            dimT(dT)=isallsame(dim2varFT(dim2varT==dT));
                            varT(dim2varT==dT)=dimT(dT);
                        end
                        d2stay=all([varT;varFT]);
                        dimFT2stay=unique(dim2varFT(d2stay));
                        dimFT2change=unique(dim2varFT(~d2stay));
                        dimT2stay=unique(dim2varT(d2stay));
                        dimT2change=unique(dim2varT(~d2stay));
                        
                        % Si tout est different :
                        if length(dimFT2change)==F.dim
                            H=HSEPfils(SEPfils(F),T);
                            return
                        end
                        
                        % Creation de l'arbre Tstay et Tchange :
                        Cvar=getCvar(T);
                        Nc=[];Ns=[];
                        search=find(getconnect(T)==1);
                        for d=1:max(dim2varT)
                            locD=find(getvar2dim(T,1:totaldim(F),1)==d)';
                            for n=search
                                if length(Cvar{n})==length(locD) &&...
                                        all(Cvar{n}==locD)
                                    if ~isempty(find(dimT2change==d, 1))
                                        Nc=[Nc,n];
                                    elseif ~isempty(find(dimT2stay==d, 1))
                                        Ns=[Ns,n];
                                    end
                                    break
                                end
                            end
                        end
                        STc=cell(1,length(Nc));
                        STs=cell(1,length(Ns));
                        for n=1:length(Nc), STc{n}=subtree(T,Nc(n));end
                        for n=1:length(Ns), STs{n}=subtree(T,Ns(n));end
                        STc=TREE(STc);
                        STs=TREE(STs);
                        
                        % Changer tous les rangs :
                        newF=cell(0,max(dim2varT));
                        newA=[];
                        for m=1:F.m
                            % F2change
                            F2change=HSEPfils(F.F(m,dimFT2change));
                            F2change=HSEPfils(F2change,STc);
                            newm    =F2change.m;
                            F2change=F2change.F;
                            
                            % F2stay
                            F2stay  =F.F(m,dimFT2stay);
                            F2stay  =HSEPfils(HSEPfils(F2stay),STs);
                            F2stay  =repmat(F2stay.F,[newm 1]);
                            
                            % Finalisation :
                            newAm   =repmat(F.alpha(m),[1 newm]);
                            newFm   =cell(newm,max(dim2varT));
                            newFm(:,dimT2stay) = F2stay;
                            newFm(:,dimT2change) = F2change;
                            newF=[newF ; newFm ];
                            newA=[newA newAm];
                        end
                        
                        % Remplissage de H :
                        H=HSEPfils();
                        H.tree=T;
                        H.F=newF;
                        H.alpha=newA;
                        H.m=length(newA);
                        H.dim=max(dim2varT);
                    end
                    
                else
                    % Transformer chaque rang de F en HSM
                    if F.m==0
                        H=HSEP();
                        H.tree=alpha;
                        return
                    end
                    connect = getconnect(T);
                    subnode = find(connect==1);
                    Cvar    = getCvar(T);
                    D=length(subnode);
                    SEPfils=str2func(class(F));
                    HSEPfils=str2func([ 'H' class(F) ]);
                    SH=cell(F.m,D);
                    for d=1:D
                        extrdim = Cvar{subnode(d)};
                        ST      = subtree(T,subnode(d));
                        if length( find(getconnect(ST)==1)) < length(extrdim)
                            %                             SH(:,d)=cellfun(@(FF) HSEP(SEP(FF),ST),F.F(:,extrdim),'UniformOutput',0);
                            for r=1:F.m
                                SH{r,d} = HSEPfils(SEPfils(F.F(r,extrdim)),ST);
                            end
                        else
                            %                             SH(:,d)=cellfun(@(FF) SEP(FF),F.F(:,extrdim),'UniformOutput',0);
                            for r=1:F.m
                                SH{r,d} = SEPfils( F.F(r,extrdim) );
                            end
                        end
                    end
                    H=HSEPfils(SH);
                    H.alpha = F.alpha;
                end
                
                
            else
                error('Non prevu')
            end
            
        end
    end
end

function s=isallsame(a)
s=min(a)==max(a);
end



