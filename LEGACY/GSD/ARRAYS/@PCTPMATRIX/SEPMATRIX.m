function x = SEPMATRIX(x)

try
phi = get_ximasse(x);
catch
phi = get_ximasse(calc_ximasse(x));    
end

x = SEPMATRIX([{x.phi0},phi]);

