function plottext(u,texto,varargin)
% function plottext(u,texto,varargin)

if isa(texto,'double') && numel(u)==numel(texto)
    texto = texto(:);
end

if ~isa(texto,'char')
    texto = num2str(texto);
end

coord = double(getcoord(u));
coord = permute(coord,[1,3,2]);
coord = reshape(coord,[size(coord,1)*size(coord,2),size(coord,3)]);

switch size(coord,2)
    case 1
        coord = [coord,zeros(size(coord,1),1)];
        text(coord(:,1),coord(:,2),texto,varargin{:});
    case 2
        text(coord(:,1),coord(:,2),texto,varargin{:});
    case 3
        text(coord(:,1),coord(:,2),coord(:,3),texto,varargin{:});
end
