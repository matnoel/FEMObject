function rep = isenrich(u,num,choix)

if nargin==1
  rep = zeros(1,u.nbnode);
  rep(u.lsenrichnode)=1;
else
    if nargin==2
        choix='global';
    end
    if strcmp(choix,'global')
     pos = getpos(u,num);
    elseif strcmp(choix,'local')
     pos = num;
    else
        error('pas prevu')
    end
    rep = ismember(pos,u.lsenrichnode);
    rep = reshape(rep,size(num));
end