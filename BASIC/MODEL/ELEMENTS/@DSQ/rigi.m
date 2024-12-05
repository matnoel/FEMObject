function ke = rigi(elem,node,varargin)
% function ke = rigi(elem,node,varargin)

xnode = node(elem);

kem = integrate(elem,xnode,2,@eval_kem);
keb = integrate(elem,xnode,4,@eval_keb);
[Bsb,Bsa,Ab] = calc_Bs(elem,xnode);
kesb = integrate(elem,xnode,0,@eval_kesb,Bsb);
kesba = integrate(elem,xnode,0,@eval_kesba,Bsb,Bsa);
kesa = integrate(elem,xnode,0,@eval_kesa,Bsa);
kebs = keb(1:12,1:12) + kesb + Ab'*(keb(13:16,13:16)+kesa)*Ab + ...
    (keb(1:12,13:16)+kesba)*Ab + Ab'*(keb(1:12,13:16)'+kesba'); 
ked = integrate(elem,xnode,2,@eval_ked);

nbddl = getnbddl(elem);
nbelem = getnbelem(elem);
P = calc_P(elem);

repm = [1,2];repm=[repm,repm+6,repm+12,repm+18];
repbs = [3,4,5];repbs=[repbs,repbs+6,repbs+12,repbs+18];
repd = [6:6:24];

ke = zerosND(nbddl,nbddl,nbelem);
ke(repm,repm) = kem;
ke(repbs,repbs) = kebs;
ke(repd,repd) = ked;

ke = P'*ke*P;


function ke = eval_kem(xi,elem,xnode)
mat = getmaterial(elem);
D = calc_opmatmembrane(mat,elem,xnode,xi);
B = calc_Bm(elem,xnode,xi);
ke = B'*D*B;
return

function ke = eval_keb(xi,elem,xnode)
mat = getmaterial(elem);
D = calc_opmatbending(mat,elem,xnode,xi);
B = calc_Bb(elem,xnode,xi);
ke = B'*D*B;
return

function ke = eval_kesb(xi,elem,xnode,Bsb)
mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnode,xi);
ke = Bsb'*(Hs\Bsb);
return

function ke = eval_kesba(xi,elem,xnode,Bsb,Bsa)
mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnode,xi);
ke = Bsb'*(Hs\Bsa);
return

function ke = eval_kesa(xi,elem,xnode,Bsa)
mat = getmaterial(elem);
Hs = calc_opmatshear(mat,elem,xnode,xi);
ke = Bsa'*(Hs\Bsa);
return

function ke = eval_ked(xi,elem,xnode)
mat = getmaterial(elem);
D = calc_opmatdrilling(mat,elem,xnode,xi);
B = calc_Bd(elem,xnode,xi);
ke = B'*D*B;
return
