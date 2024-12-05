function u = setphi(u,phi,i)


if nargin==2 && isa(phi,'cell')
u.phi=phi;
u.isranddim = getisranddim(phi);
u.ximasse={};
elseif nargin==2 && isa(phi,'PCTPMATRIX') && getnbdim(u)==getnbdim(phi)
u.phi=phi.phi;
u.isranddim = phi.isranddim;
u.ximasse = phi.ximasse;
elseif nargin==3 && isa(phi,'cell') && length(i)==length(phi)
u.phi(i) = phi;
u.isranddim = getisranddim(phi);
elseif nargin==3 && isa(phi,'double') && length(i)==1
u.phi{i} = phi;
u.isranddim(i) = (numel(phi)>1);
elseif nargin==3 && isa(phi,'PCTPMATRIX')
u.phi(i) = phi.phi(i);
u.isranddim(i) = phi.isranddim(i);

else
    error('pas prevu')
end

    
function isr = getisranddim(phi)

isr = zeros(1,length(phi));
for k=1:length(phi)
    isr(k)= ~(numel(phi{k})==1);    
end
    
return
