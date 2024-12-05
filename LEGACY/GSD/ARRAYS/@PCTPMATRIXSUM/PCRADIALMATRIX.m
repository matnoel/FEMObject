function y = PCRADIALMATRIX(x)
% function y = PCRADIALMATRIX(x)

y = PCRADIALMATRIX(x.funs{1});
for i=2:length(x.funs)
y = y + PCRADIALMATRIX(x.funs{i});    
end