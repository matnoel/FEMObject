function x = SEPMATRIX(x)

try
B = getximasse(x);
catch
B = getximasse(calc_ximasse(x));    
end

x = mat2cell(x);
A = getvalue(getV(x));
B = getvalue(mat2cell(B));

x = SEPMATRIX([A,B]);
