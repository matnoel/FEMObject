function hval = polyval(h,liste,x)
% hval = polyval(h,liste,x)
%
% Evaluate a RANDOM POLYNOMIAL h_n(x)
%

x = transfer(h.rv,RANDVAR(h.POLYFE),x);
hval = polyval(h.POLYFE,liste,x);
