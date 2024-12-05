function [sto,PC] = PCMODEL(a,varargin)
% function sto = PCMODEL(a)
% a: PCMATRIX
% cree un model stochastique
%
% See also RANDVARS/PCMODEL

if nargin==0
    sto.X = PCMATRIX();
    PC = POLYCHAOS();
    sto=class(sto,'PCMODEL',PC);
elseif nargin==1  && isa(a,'PCMODEL')
    sto=a;
    PC = getPC(a);
elseif isa(a,'PCMATRIX')
    PC = getPC(a);

    sto.X = a(:);
    if ~ischarin('nomasse',varargin)
        sto.X = calc_ximasse(sto.X);
    end

    sto=class(sto,'PCMODEL',PC);
    superiorto('POLYCHAOS')

elseif isa(a,'RANDVAR')
    [sto,PC] = PCMODEL(RANDVARS(a),varargin{:});
elseif isa(a,'RANDPOLYS') 
    [sto,PC] = PCMODEL(RANDVARS(a),varargin{:});    
else
    error('mauvais arguments')
end
