function elem=actualise_ddl(elem,varargin)
% function elem=actualise_ddl(elem,ddlnode,ddlgauss,ddlnodedual,ddlgaussdual)


if nargin>=2 && isa(varargin{1},'DDL')
    ddlnode = varargin{1};
    if nargin>=3
        ddlnodedual = varargin{2};
    else
        ddlnodedual = ddlnode ;
    end
    if nargin>4
        ddlgauss = varargin{3};
        ddlgaussdual = varargin{4};
    else
        ddlgauss = ddlnode;
        ddlgaussdual = ddlnodedual ;
    end
else
    
    try
        [ddlnode,ddlgauss,ddlnodedual,ddlgaussdual]=calc_ddl(getmaterial(elem),elem);
%         if getlsenrich(elem)>0
%             ddlnode = enrich(ddlnode);
%             ddlnodedual = enrich(ddlnodedual);
%         end
    catch
        ddlnode = DDL('u');
        ddlnodedual = DDL('q');
        ddlgauss = ddlnode;
        ddlgaussdual = ddlnodedual;
%         warning('primal dof set to ''u'', dual dof set to ''q''');
    end
    
    
end

elem.ddlnode = DDL(ddlnode);
elem.ddlgauss = DDL(ddlgauss);
elem.ddlnodedual = DDL(ddlnodedual);
elem.ddlgaussdual = DDL(ddlgaussdual);

elem.nbddlpernode=length(ddlnode);
elem.nbddlpergauss=length(ddlgauss);
elem.nbddl=elem.nbddlpernode*getnbnode(elem);
