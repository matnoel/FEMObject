function s = expand(s)

if getnbdim(s)==1
    s = s.factor*s.phi{1};
end