function O=one(SM,varargin)
% Fonction constante :
O=cell(1,SM.dim);
for d=1:SM.dim
    if nargin==2 && isa(varargin{1},'char') && strcmp(SM.F{d}.type,'SPACE')
        o=double(ones(1,getnbelem(SM.F{d}.model)));
    elseif isa(SM.F{d}.model,'double')
        o=ones(SM.F{d}.model,1);
    else
       o=double(one(SM.F{d}.model));
    end
    O{1,d}=o(:);
end
O=SEPMATRIX(O);