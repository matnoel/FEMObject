function gauss=calc_subgausspoints(h,n,m)
% function gauss=calc_subgausspoints(h,n,m)
% calcul d'une quadrature par morceaux

domain = getdomain(h);
hfe = POLYFE(linspace(domain(1),domain(2),m+1));
gauss = calc_gausspoints(hfe,n);
gauss.w = gauss.w*1/(domain(2)-domain(1));
if ~isa(h,'POLYLEGENDRE')
    error('calcul pas correct')
end
