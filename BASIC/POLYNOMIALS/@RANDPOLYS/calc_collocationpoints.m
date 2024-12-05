function [points,points1D]=calc_collocationpoints(H,multiindex)

% function points=calc_collocationpoints(H,multiindex)
% calcul des points de collocation pour des interpolants de lagrange
% multiindex permet de déterminer comment trier les points de Gauss

M=H.M;
points1D = cell(1,M);
for k=1:M
    
    if ~(isa(H.h{k},'POLYLAGRANGE') || isa(H.h{k},'POLYFELAGRANGE'))
        error('chaos must contain POLYFELAGRANGE or POLYLAGRANGE')
    end
   points1D{k}=getpoints(H.h{k});
end

points=zeros(size(multiindex));
for k=1:M
    points(:,k) = points1D{k}(multiindex(:,k));
end

