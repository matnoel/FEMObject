function subplot(M,a,b,c,varargin)
% function subplot(M,a,b,c,'propertyname',propertyvalue,...)
% M : MODEL
% propertyname       propertyvalue        effet
% 'selgroup'         liste                affiche les groupes d'elements de liste
% 'selelem'          liste                affiche les elements liste
% 'node'             ---                  affiche les noeuds avec des petits points
% 'numnode'          ---                  affichage des numeros des noeuds
% 'numelem'          ---                  affichage des numeros des elements
% 'mat'              ---                  affichage des materiaux
% 'group'            ---                  affichage des groupes d'elements
% 'color'            ---                  definit la couleur du maillage
% -> voir patchoptions pour la suite des parametres associes a la commande
% patch

issurf = ischarin('surface',varargin);
varargin = delonlycharin('surface',varargin);

nodetext = ischarin('numnode',varargin);
elemtext = ischarin('numelem',varargin);
matplot = ischarin('mat',varargin);
lsenrichplot = ischarin('lsenrich',varargin);
groupplot = ischarin('group',varargin);
listegroup = getcharin('selgroup',varargin,1:M.nbgroupelem);
listeelem = getcharin('selelem',varargin,getnumelem(M));

dim = getindim(M);
options = patchoptions(dim,varargin{:});

color = getcharin('color',varargin);
if ~isempty(color) && isa(color,'char')
    options = setcharin('edgecolor',options,color);
elseif ~isempty(color) && isa(color,'double')
    options = setcharin('facecolor',options,'flat');
    options = setcharin('facevertexcdata',options,color);
end
node = M.node;

Helem = [];
Hnode = [];
for p=listegroup
    if ischarin('selelem',varargin)
        if ischarin('local',varargin)
            elem = getelem(M.groupelem{p},listeelem,'local');
        else
            elem = getelem(M.groupelem{p},listeelem,'global');
        end
    else
        elem = M.groupelem{p};
    end
    optionselem = options;
    if getdim(elem)==1
        optionselem = delcharin('facelighting',optionselem);
        optionselem = delcharin('edgelighting',optionselem);
    end
    if getnbelem(elem)>0 && (~strcmp(getlstype(elem),'out') || ischarin('out',varargin))
        
        if matplot || groupplot || lsenrichplot
            if matplot
                colelem = getmaterialnumber(elem);
            elseif groupplot
                colelem = p;
            elseif lsenrichplot
                colelem = getlsenrich(elem);
            end
            optionselem = setcharin('facecolor',optionselem,'flat');
            optionselem = setcharin('facevertexcdata',optionselem,colelem);
        end
        subplot(a,b,c);
        plot(elem,node,optionselem{:});
        
        if elemtext==1
            plotnumber(elem,node,'color','r',varargin{:});
        end
        
        if  matplot || groupplot || lsenrichplot
            connec = getconnec(elem);
            xplot = double(mean(getcoord(node,connec),3));
            xplot = POINT(mean(xplot,1));
            plottext(xplot,num2str(colelem),'Color','w','FontSize',14);
        end
        
    end
end

if nodetext==1
    plotnumber(node,varargin{:});
end

axis image
axis off

if dim==3 || issurf
    view(3)
end

if nargout>=1
    varargout{1} = Helem;
end
