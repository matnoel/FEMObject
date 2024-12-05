function H = calc_historyfield(S,q,H_old,varargin)
% function H = calc_historyfield(S,q,H_old,varargin)
% S : MODEL
% q : solution
% H_old : previous field

if ~isa(S,'MODEL')
    H = calc_historyfield(q,S,H_old,varargin{:});
    return
end

if isa(H_old,'FENODEFIELD')
    h_old = double(H_old);
    H = FENODEFIELD(calc_energyint(S,q,'node','positive','local'));
    h = double(H);
    rep = find(h <= h_old);
    h(rep) = h_old(rep);
else
    h_old = getvalue(H_old);
    H = calc_energyint(S,q,'intorder','mass','positive','local');
    h = getvalue(H);
    for p=1:S.nbgroupelem
        he = double(h{p});
        he_old = double(h_old{p});
        rep = find(he <= he_old);
        he(rep) = he_old(rep);
        h{p} = MYDOUBLEND(he);
    end
end
H = setvalue(H,h);
