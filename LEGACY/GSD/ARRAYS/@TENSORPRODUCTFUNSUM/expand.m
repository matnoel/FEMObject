function z = expand(s)

if getnbdim(s)<=1
z = expand(s.tensorfuns{1});    
for i=2:length(s.tensorfuns)
z = z + expand(s.tensorfuns{i});
end
end


