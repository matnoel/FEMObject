function x = full(x)

x.phi0=full(x.phi0);
for i=1:length(x.phi)
    x.phi{i} = full(x.phi{i});
end

