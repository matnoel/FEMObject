function elem=setddl(elem,varargin)
% function elem=setddl(elem,'ddlnode',ddlnode,'ddlnodedual',ddlnodedual,...)

elem.ddlnode=DDL(getcharin('ddlnode',varargin,DDL()));
elem.nbddlpernode=length(elem.ddlnode);
elem.ddlgauss=DDL(getcharin('ddlgauss',varargin,elem.ddlnode));
elem.nbddlpergauss=length(elem.ddlgauss);
elem.ddlnodedual=DDL(getcharin('ddlnodedual',varargin,elem.ddlnode));
elem.ddlgaussdual=DDL(getcharin('ddlgaussdual',varargin,elem.ddlgauss));

elem.nbddl=elem.nbddlpernode*getnbnode(elem);
elem.numddl=[];
