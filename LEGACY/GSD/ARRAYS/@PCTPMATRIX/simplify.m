function x = simplify(x)

isr = isranddim(x);
rdim = find(isr);
ddim = find(~isr);

if isempty(rdim)
phi = prod(vertcat(x.phi{:}),1);
x = x.phi0 * phi;
elseif length(rdim)==1 && numel(x)==1
phi = prod(vertcat(x.phi{ddim}),1);
x = x.phi0 * phi * x.phi{rdim};
end
  

