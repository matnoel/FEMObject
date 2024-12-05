function xix = funtransfer(a,ax)
% function xix = funtransfer(a,ax)
% a : RANDVAR    
% xix=Phi(F_a(ax))
xi=RVNORMAL();
xix=transfer(a,xi,ax);
