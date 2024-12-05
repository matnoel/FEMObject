function ke = rigi(elem,node,varargin)
% function ke = rigi(elem,node,varargin)

xnode = node(elem);

kem = integrate(elem,xnode,2,@eval_kem);
keb = integrate(elem,xnode,2,@eval_keb);

if getparam(elem,'fullintegration')==1
    kes = integrate(elem,xnode,2,@eval_kes);
else
    kes = integrate(elem,xnode,1,@eval_kes);
end
ked = integrate(elem,xnode,2,@eval_ked);

nbddl = getnbddl(elem);
nbelem = getnbelem(elem);
P = calc_P(elem);

repm = [1,2];repm=[repm,repm+6,repm+12,repm+18];
repb = [4,5];repb=[repb,repb+6,repb+12,repb+18];
reps = [3,4,5];reps=[reps,reps+6,reps+12,reps+18];
repd = [6:6:24];

ke = zerosND(nbddl,nbddl,nbelem);
ke(repm,repm) = kem;
ke(reps,reps) = kes;
ke(repb,repb) = ke(repb,repb)+keb;
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

function ke = eval_kes(xi,elem,xnode)
mat = getmaterial(elem);
D = calc_opmatshear(mat,elem,xnode,xi);
B = calc_Bs(elem,xnode,xi);
ke = B'*D*B;
return

function ke = eval_ked(xi,elem,xnode)
mat = getmaterial(elem);
D = calc_opmatdrilling(mat,elem,xnode,xi);
B = calc_Bd(elem,xnode,xi);
ke = B'*D*B;
return
