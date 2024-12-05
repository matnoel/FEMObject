function T = TREE(connect,varargin)

% Arbre vide
if nargin==0
    T.connect = [];
    T.var     = [];
    T.Cvar    = {};
    T.dim     = [];
    T.var2dim = {};
    T.leaf    = [];
    T.param   = struct();
    T = class(T,'TREE');

% Copie d'arbre
elseif nargin==1 && isa(connect,'TREE') %% Un classique sans intéret...
    T = connect;

% Creation a partir de connect
elseif isa(connect,'double')
    T         = TREE();
    T.connect = connect;
    [T.dim,T.var,T.var2dim,T.leaf,T.Cvar] = init_tree(T);

% Assemblage d'arbre
elseif nargin>1 && isa(connect,'TREE') && isa(varargin{1},'TREE')
    trees = cell(1+length(varargin),1);
    trees{1} = connect;
    for i=1:length(varargin)
        trees{1+i} = varargin{i};
    end
    T = TREE();
    T = assembled_tree(T,trees);

% Assemblage d'arbre, a partir d'une cellule d'arbre 
elseif nargin==1 && isa(connect,'cell') && isa(connect{1},'TREE')
    T = TREE();
    T = assembled_tree(T,connect);
end






