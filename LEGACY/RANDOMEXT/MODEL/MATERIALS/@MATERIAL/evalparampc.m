function varargout = evalparampc(mat,name,PC,varargin)
% function [val,L,x] = evalparampc(mat,name,PC,elem,xnode,xgauss)
% function [val,L,x] = evalparampc(mat,name,PC,x)
% mat : MATERIAL
% PC : POLYCHAOS
% name : nom du parametre à évaluer
% x : points en lesquels on evalue le parametre
% Plusieurs cas possibles
% - Si le parametre est deterministe, val est l'evaluation du parametre
%   et L est vide
% - si le parametre est une PCMATRIX de taille [1,1]
%    val=1 et L est la variable aléatoire (ou la va constante exprimee sur PC)
%    exprimée sur PC
% - si le parametre est un FENODEFIELD
%   1) si le FENODEFIELD est une PCRADIALMATRIX   sum_i(U_i lambda_i)
%      val est un MYDOUBLEND contenant l'estimation des U_i en x. La 5eme
%      dimension aux differents U_i
%      L est la PCMATRIX des lambda_i
%   2) si le FENODEFIELD est une PCMATRIX  sum_i(U_i H_i) (H_i fonctions de base de PC)
%      val est un MYDOUBLEND contenant l'estimation des des U_i en x. La 5eme
%      dimension aux differents U_i
%      L est le POLYCHAOS PC

param=getparam(mat,name);

if nargin ==4
    x=varargin{1};
else
    elem=varargin{1};
    xnode=varargin{2};
    xgauss=varargin{3};
    [x,Nuni]=calc_x(elem,xnode,xgauss);    
end

    if isa(param,'FENODEFIELD')
    param = getvalue(param);
    if isa(param,'PCRADIALMATRIX')
        L = getL(param);
        param = getV(param);
        if iscell(param)
        param=cell2mat(param);
        end
        param = double(param);        
        
    elseif isa(param,'PCMATRIX')
        L = getPC(param);
        param = double(param);
        
    elseif isa(param,'double')
        param = double(param);
        L = [];
    else
       error('pas programme') 
    end
    
    param=full(param);
    con = getconnec(elem)';
    param=param(con(:),:);
    param=reshape(param,[size(con,1) 1 size(con,2) 1 size(param,2)]);
    
    val = Nuni*MYDOUBLEND(param);
    if nargout<=1 && ~isempty(L) %%%%%%
        val = PCMYDOUBLEND(val,L,5);
    end
    
    elseif isa(param,'double') || isa(param,'MYDOUBLE')
    val = double(param) ;
    L = [];
    elseif isa(param,'PCMATRIX')
    
    if nargout<=1
    val = convertradial(PCMYDOUBLEND(param,5));
    
    else
    if numel(param)==1
    val = 1;
    L = param; 
    else
    val = zerosND(numel(param),1,1,1,numel(param));
    for k=1:numel(param)
    val(k,1,1,1,k)=1;
    end
    L = param;
    end
        
    end
    
    
    elseif isa(param,'RANDVAR')
    error('decomposer le MATERIAL sur le chaos')    
       
    end
    
    varargout{1}=val;
    if nargout>1
    varargout{2}=L ; 
    end
    if nargout>2
      varargout{3}=x;  
    end
 