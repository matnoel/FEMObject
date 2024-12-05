function x = PCTPMATRIX(PC,phi0,varargin)
% function x = PCTPMATRIX(PC,phi0,varargin)

if nargin==0
    x.phi0 = [];
    x.phi=cell(1,0);
    x.isranddim = zeros(1,0);
    x.ximasse = {};
    PC = POLYCHAOSTP();
    x = class(x,'PCTPMATRIX',PC);
elseif nargin==2
    x.phi0 = phi0;
    x.phi  = getphi(one(PC));%repmat({1},1,getnbgroups(PC));
    x.isranddim = ones(1,getnbgroups(PC));
    x.ximasse = {};
    x = class(x,'PCTPMATRIX',PC);
elseif ischarin('dim',varargin)
    dim = getcharin('dim',varargin);  
%phi = getphi(one(PC));
    phi = cell(1,getnbgroups(PC));
    phi(:)={1};
    phi(dim) = varargin(1:length(dim)); 

    phi = checkphi(phi,PC);
    x.isranddim = ones(1,getnbdim(PC));
    x.isranddim(dim) = getisranddim(phi(dim));
    x = PCTPMATRIX(PC,phi0,phi{:});
else
    x.phi0 = phi0;

    x.phi = cell(1,getnbgroups(PC));
    phi = checkphi(varargin(1:getnbgroups(PC)),PC);
    x.phi = phi;
    x.isranddim = getisranddim(phi);
    x.ximasse = {};
    x = class(x,'PCTPMATRIX',PC);

end


function phi = checkphi(phi,PC)

for i=1:length(phi)
    if ~isa(phi{i},'double') || (size(phi{i},1)~= length(PC,i) && size(phi{i},1)>1)
        error('les arguments doivent etre des double')
    end
end

return

function isr = getisranddim(phi)

isr = zeros(1,length(phi));
for k=1:length(phi)
    isr(k)= ~(numel(phi{k})==1);    
end

return
